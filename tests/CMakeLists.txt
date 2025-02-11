enable_testing()
find_package(GTest REQUIRED)

add_executable(vectors_test vectors_test.cpp)
target_link_libraries(vectors_test PRIVATE physics_engine GTest::GTest GTest::Main)

add_test(NAME Vector2D_Test COMMAND vectors_test)

