add_library(obb STATIC obb.cpp obb.h)
target_link_libraries(obb PUBLIC Eigen3::Eigen)
target_include_directories(obb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})