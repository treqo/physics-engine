enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(vectors_test vectors_test.cpp)
target_link_libraries(vectors_test PRIVATE physics_engine GTest::GTest GTest::Main)
gtest_discover_tests(vectors_test)

add_executable(matrix_test matrix_test.cpp)
target_link_libraries(matrix_test PRIVATE physics_engine GTest::GTest GTest::Main)
gtest_discover_tests(matrix_test)