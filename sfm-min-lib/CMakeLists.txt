add_library(sfm-min-lib optim.cpp)
add_library(graph graph.cpp)


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_library(linalg linalg.cpp)
target_link_libraries(linalg ${ARMADILLO_LIBRARIES})
