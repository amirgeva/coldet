add_library(obb STATIC 
	obb.cpp
	obb.h
	obbtree.cpp
	obbtree.h
)
target_link_libraries(obb PUBLIC coldet Eigen3::Eigen)
target_include_directories(obb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})