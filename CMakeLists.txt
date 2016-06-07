cmake_minimum_required(VERSION 2.8)

project(sfm-min)

add_subdirectory(sfm-min-lib)
include_directories(include)


add_executable(sfm-min-run main.cpp)
target_link_libraries(sfm-min-run sfm-min-lib)
