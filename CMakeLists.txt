cmake_minimum_required(VERSION 2.8)

project(sfm-min)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library code
include_directories(include)
add_subdirectory(sfm-min-lib)


# Submodular function minimisation application
add_subdirectory(min-cut-pb)
include_directories(min-cut-pb)

# Executable
add_executable(sfm-min-run main.cpp)
target_link_libraries(sfm-min-run sfm-min-lib min-cut linalg)
