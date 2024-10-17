#include <autotuner.h>

#include <iostream>

size_t calc_first_level_in_sd(const rocksdb::Options &options) {
  uint64_t fd_size = options.db_paths[0].target_size;
  size_t level = 0;
  uint64_t level_size = options.max_bytes_for_level_base;
  while (level_size <= fd_size) {
    fd_size -= level_size;
    if (level > 0) {
      level_size *= options.max_bytes_for_level_multiplier;
    }
    level += 1;
  }
  return level;
}

// Solve the equation: x^a + x^{a-1} + ... + x = b
double calc_size_ratio(size_t a, double b) {
  assert(a > 0);
  double inv_a = 1 / (double)a;
  // Let f(x) = x^a + x^{a-1} + ... + x - b = 0
  // x^a + x^{a-1} + ... + 1 = (x^{a+1} - 1) / (x - 1)
  // x^a + x^{a-1} + ... + x = (x^{a+1} - 1) / (x - 1) - 1
  // (x^{a+1} - 1) / (x - 1) - 1 = b
  // x^{a+1} - 1 = (x - 1) * (b + 1)
  // x^{a+1} - 1 = (b + 1) x - b - 1
  // x^{a+1} - (b + 1) x + b = 0
  // Let g(u) = u^{a+1} - (b + 1) u + maxb
  // g'(u) = (a + 1) u^a - b - 1
  // Let g'(u_min) = 0, then u_min = ((b + 1) / (a + 1)) ^ (1 / a)
  // So x >= ((b + 1) / (a + 1)) ^ (1 / a)
  // x^a + x^{a-1} + ... + x = b > x^a
  // So x < b ^ (1 / a)
  // In conclusion, ((b + 1) / (a + 1)) ^ (1 / a) <= x < b ^ (1 / a)
  double min = pow((b + 1) / (a + 1), inv_a);
  double max = pow(b, inv_a);
  auto f = [a, b](double x) {
    double sum = 0;
    double xa = x;
    for (size_t i = 1; i <= a; ++i) {
      sum += xa;
      xa *= x;
    }
    return sum - b;
  };
  while (max - min >= 0.001) {
    double x = (max + min) / 2;
    if (f(x) > 0) {
      max = x;
    } else {
      min = x;
    }
  }
  return (max + min) / 2;
}
void calc_fd_size_ratio(rocksdb::Options &options, size_t first_level_in_sd,
                        uint64_t max_ralt_size) {
  options.max_bytes_for_level_multiplier_additional.clear();
  // It seems that L0 and L1 are not affected by
  // options.max_bytes_for_level_multiplier_additional
  if (first_level_in_sd <= 2) return;

  assert(options.db_paths[0].target_size > max_ralt_size);
  uint64_t fd_size = options.db_paths[0].target_size - max_ralt_size;
  double ratio =
      calc_size_ratio(first_level_in_sd - 2,
                      (double)(fd_size - 2 * options.max_bytes_for_level_base) /
                          options.max_bytes_for_level_base);
  // Multiply 0.999 to make room for floating point error
  ratio *= 0.999;
  assert(options.max_bytes_for_level_multiplier_additional.empty());
  options.max_bytes_for_level_multiplier_additional.push_back(1.0);
  for (size_t i = 2; i < first_level_in_sd; ++i) {
    options.max_bytes_for_level_multiplier_additional.push_back(
        ratio / options.max_bytes_for_level_multiplier);
  }
}
void calc_sd_size_ratio(rocksdb::Options &options, rocksdb::DB &db,
                        size_t last_level_in_fd, uint64_t last_level_in_fd_size,
                        uint64_t hot_set_size) {
  std::string str;
  db.GetProperty(rocksdb::DB::Properties::kLevelStats, &str);
  std::istringstream in(str);
  // The first two lines are headers.
  size_t lines_to_skip = 2 + last_level_in_fd + 1;
  while (in && lines_to_skip) {
    --lines_to_skip;
    while (in && in.get() != '\n')
      ;
  }
  if (!in) return;

  size_t last_level = last_level_in_fd;
  uint64_t sd_level_size = 0;
  while (in) {
    size_t level;
    size_t num_files;
    size_t size;
    in >> level >> num_files >> size;
    if (size == 0) break;
    last_level = level;
    sd_level_size += size;
  }
  sd_level_size *= 1048576;
  if (last_level <= last_level_in_fd + 1) return;
  // unlikely
  if (sd_level_size <= last_level_in_fd_size) return;

  assert(last_level_in_fd > 1);  // Not implemented yet
  uint64_t last_level_in_fd_effective_size =
      last_level_in_fd_size - hot_set_size;
  size_t a = last_level - last_level_in_fd;
  double b = (double)sd_level_size / last_level_in_fd_effective_size;
  double sd_ratio = calc_size_ratio(a, b);
  size_t level = last_level_in_fd + 1;
  if (sd_ratio > options.max_bytes_for_level_multiplier) {
    do {
      last_level += 1;
      a += 1;
      sd_ratio = calc_size_ratio(a, b);
    } while (sd_ratio > options.max_bytes_for_level_multiplier);
  }
  // The first level in SD shouldn't be smaller than the last level in FD.
  if (last_level_in_fd_effective_size * sd_ratio <
      last_level_in_fd_size * 1.01) {
    double ratio_additional = 1.01 / options.max_bytes_for_level_multiplier;
    options.max_bytes_for_level_multiplier_additional.push_back(
        ratio_additional);
    uint64_t first_level_in_sd_size = last_level_in_fd_size *
                                      options.max_bytes_for_level_multiplier *
                                      ratio_additional;
    a -= 1;
    b = (double)(sd_level_size - first_level_in_sd_size) /
        first_level_in_sd_size;
    sd_ratio = calc_size_ratio(a, b);
    level += 1;

    if (sd_ratio > options.max_bytes_for_level_multiplier) {
      do {
        last_level += 1;
        a += 1;
        sd_ratio = calc_size_ratio(a, b);
      } while (sd_ratio > options.max_bytes_for_level_multiplier);
    }
  } else {
    options.max_bytes_for_level_multiplier_additional.push_back(
        last_level_in_fd_effective_size * sd_ratio / last_level_in_fd_size /
        options.max_bytes_for_level_multiplier);
    level += 1;
  }
  sd_ratio /= options.max_bytes_for_level_multiplier;
  for (; level < last_level; ++level) {
    options.max_bytes_for_level_multiplier_additional.push_back(sd_ratio);
  }
  options.max_bytes_for_level_multiplier_additional.push_back(100.0);
}
bool should_update_max_bytes_for_level_multiplier_additional(
    const std::vector<double> &ori, const std::vector<double> &cur) {
  if (ori.size() != cur.size()) {
    assert(ori.size() < cur.size());
    return true;
  }
  for (size_t i = 0; i < ori.size(); ++i) {
    if (cur[i] < ori[i] * 0.99) return true;
    if (cur[i] > ori[i] * 1.01) return true;
  }
  return false;
}

double MaxBytesMultiplerAdditional(const rocksdb::Options &options, int level) {
  if (level >= static_cast<int>(
                   options.max_bytes_for_level_multiplier_additional.size())) {
    return 1;
  }
  return options.max_bytes_for_level_multiplier_additional[level];
}

std::vector<std::pair<uint64_t, uint32_t>> predict_level_assignment(
    const rocksdb::Options &options) {
  std::vector<std::pair<uint64_t, uint32_t>> ret;
  uint32_t p = 0;
  size_t level = 0;
  assert(!options.db_paths.empty());

  // size remaining in the most recent path
  uint64_t current_path_size = options.db_paths[0].target_size;

  uint64_t level_size;
  size_t cur_level = 0;

  // max_bytes_for_level_base denotes L1 size.
  // We estimate L0 size to be the same as L1.
  level_size = options.max_bytes_for_level_base;

  // Last path is the fallback
  while (p < options.db_paths.size() - 1) {
    if (current_path_size < level_size) {
      p++;
      current_path_size = options.db_paths[p].target_size;
      continue;
    }
    if (cur_level == level) {
      // Does desired level fit in this path?
      assert(ret.size() == level);
      ret.emplace_back(level_size, p);
      ++level;
    }
    current_path_size -= level_size;
    if (cur_level > 0) {
      if (options.level_compaction_dynamic_level_bytes) {
        // Currently, level_compaction_dynamic_level_bytes is ignored when
        // multiple db paths are specified. https://github.com/facebook/
        // rocksdb/blob/main/db/column_family.cc.
        // Still, adding this check to avoid accidentally using
        // max_bytes_for_level_multiplier_additional
        level_size =
            static_cast<uint64_t>(static_cast<double>(level_size) *
                                  options.max_bytes_for_level_multiplier);
      } else {
        level_size = static_cast<uint64_t>(
            static_cast<double>(level_size) *
            options.max_bytes_for_level_multiplier *
            MaxBytesMultiplerAdditional(options, cur_level));
      }
    }
    cur_level++;
  }
  assert(ret.size() == level);
  ret.emplace_back(level_size, p);
  return ret;
}

void AutoTuner::update_thread() {
  rocksdb::Options options = db_.GetOptions();
  RALT &ralt = *static_cast<RALT *>(options.ralt);

  const uint64_t initial_max_hot_set_size = ralt.GetMaxHotSetSizeLimit();
  std::cerr << "Initial max hot set size: " << initial_max_hot_set_size
            << std::endl;

  const uint64_t initial_hot_set_size_limit = ralt.GetHotSetSizeLimit();
  std::cerr << "Initial hot set size limit: " << initial_hot_set_size_limit
            << std::endl;

  uint64_t phy_size_limit = ralt.GetPhySizeLimit();
  std::cerr << "Initial physical size limit: " << phy_size_limit << std::endl;

  const size_t last_level_in_fd = first_level_in_sd_ - 1;
  std::vector<double> ori_multiplier_additional =
      options.max_bytes_for_level_multiplier_additional;
  bool warming_up = true;
  bool first = true;
  while (!stop_signal_) {
    std::this_thread::sleep_for(std::chrono::nanoseconds(wait_time_ns_));
    if (stop_signal_) {
      break;
    }
    uint64_t real_hot_set_size = ralt.GetRealHotSetSize();
    if (ralt.DecayCount() > 10) {
      if (first) {
        first = false;
        warming_up = false;
        ralt.SetMinHotSetSizeLimit(min_hot_set_size_);
      }
      double hs_step = max_hot_set_size_ / 20.0;
      uint64_t real_phy_size = ralt.GetRealPhySize();
      std::cerr << "real_phy_size " << real_phy_size << '\n';
      auto rate = real_phy_size / (double)real_hot_set_size;
      auto delta =
          rate * hs_step;  // std::max<size_t>(rate * hs_step, (64 << 20));
      phy_size_limit = real_phy_size + delta;
      std::cerr << "rate " << rate << std::endl;
      ralt.SetPhysicalSizeLimit(phy_size_limit);
      std::cerr << "Update physical size limit: " << phy_size_limit
                << std::endl;
    }
    calc_fd_size_ratio(options, first_level_in_sd_, phy_size_limit);
    assert(first_level_in_sd_ > 0);
    uint64_t last_level_in_fd_size =
        predict_level_assignment(options)[last_level_in_fd].first;
    uint64_t min_effective_size_of_last_level_in_fd =
        last_level_in_fd_size / options.max_bytes_for_level_multiplier;
    // to avoid making the size of the first level in the slow disk too small
    uint64_t max_hot_set_size =
        last_level_in_fd_size - min_effective_size_of_last_level_in_fd;
    if (warming_up) {
      max_hot_set_size = std::min(max_hot_set_size, initial_max_hot_set_size);
    } else {
      max_hot_set_size = std::min(max_hot_set_size, max_hot_set_size_);
    }
    if (ralt.GetMaxHotSetSizeLimit() != max_hot_set_size) {
      std::cerr << "Update max hot set size limit: " << max_hot_set_size
                << std::endl;
      ralt.SetMaxHotSetSizeLimit(max_hot_set_size);
    }

    uint64_t hot_set_size;
    if (warming_up) {
      hot_set_size = real_hot_set_size;
    } else {
      hot_set_size = ralt.GetHotSetSizeLimit();
      std::cerr << "hot set size limit: " << hot_set_size << std::endl;
    }
    calc_sd_size_ratio(options, db_, last_level_in_fd, last_level_in_fd_size,
                       hot_set_size);
    if (should_update_max_bytes_for_level_multiplier_additional(
            ori_multiplier_additional,
            options.max_bytes_for_level_multiplier_additional)) {
      ori_multiplier_additional.assign(
          options.max_bytes_for_level_multiplier_additional.begin(),
          options.max_bytes_for_level_multiplier_additional.end());
      std::ostringstream out;
      for (size_t i = 0; i < ori_multiplier_additional.size(); ++i) {
        out << ori_multiplier_additional[i];
        if (i != ori_multiplier_additional.size()) {
          out << ':';
        }
      }
      std::string str = out.str();
      std::cerr << "Update max_bytes_for_level_multiplier_additional: " << str
                << std::endl;
      db_.SetOptions({{"max_bytes_for_level_multiplier_additional", str}});
    }
  }
}
