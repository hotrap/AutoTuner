cmake_minimum_required(VERSION 3.15)
project(autotuner CXX)

aux_source_directory(src SRCS)
add_library(${PROJECT_NAME} ${SRCS})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_include_directories(${PROJECT_NAME} PUBLIC include)

include_directories(
	${ROCKSDB_INCLUDE}
    ${VISCNTS_INCLUDE}
)

file(GLOB_RECURSE HEADERS "include/*")
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

install(TARGETS ${PROJECT_NAME} DESTINATION "include/")
