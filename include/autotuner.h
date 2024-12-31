#ifndef AUTOTUNER_H_
#define AUTOTUNER_H_

#include <rocksdb/db.h>

#include <cmath>
#include <sstream>
#include <thread>

size_t calc_first_level_in_last_tier(const rocksdb::Options &options);

void calc_fd_size_ratio(rocksdb::Options &options, size_t first_level_in_sd,
                        uint64_t max_ralt_size);

std::vector<std::pair<uint64_t, uint32_t>> predict_level_assignment(
    const rocksdb::Options &options);

class AutoTuner {
 public:
  AutoTuner(rocksdb::DB &db, size_t first_level_in_sd,
            uint64_t min_hot_set_size, double max_hot_ratio_in_last_level_in_fd,
            uint64_t max_unstable_record_size, size_t wait_time_ns = 20e9)
      : db_(db),
        first_level_in_sd_(first_level_in_sd),
        wait_time_ns_(wait_time_ns),
        min_hot_set_size_(min_hot_set_size),
        max_hot_ratio_in_last_level_in_fd_(max_hot_ratio_in_last_level_in_fd),
        max_unstable_record_size_(max_unstable_record_size),
        ralt_(db.GetOptions().ralt) {
    th_ = std::thread([&]() { update_thread(); });
  }

  ~AutoTuner() { Stop(); }

  void Stop() {
    stop_signal_ = true;
    th_.join();
  }

 private:
  void update_thread();

  rocksdb::DB &db_;
  const size_t first_level_in_sd_;

  ssize_t wait_time_ns_;
  uint64_t min_hot_set_size_;
  double max_hot_ratio_in_last_level_in_fd_;
  uint64_t max_unstable_record_size_;

  // To make sure that RALT is not deleted while the autotuner is running.
  std::shared_ptr<rocksdb::RALT> ralt_;

  bool stop_signal_{false};
  std::thread th_;
};

#endif  // AUTOTUNER_H_
