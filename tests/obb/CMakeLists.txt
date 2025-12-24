add_executable(test.obb main.cpp)
target_link_libraries(test.obb PRIVATE coldet obb GTest::gtest_main)
