add_executable(test.obb
	main.cpp
	tritri_test.cpp
	scene_test.cpp
	coldet_api_test.cpp
	obb_box_test.cpp
	tritri_equivalence_test.cpp
	../reference/tritri_reference.c
	test_helpers.h
)
target_link_libraries(test.obb PRIVATE coldet obb GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test.obb DISCOVERY_TIMEOUT 60 DISCOVERY_MODE PRE_TEST)
