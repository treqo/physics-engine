cmake_minimum_required(VERSION 3.10)
project(PhysicsEngine)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Option to build only tests (useful for CI)
option(BUILD_TESTS_ONLY "Build only the tests without SFML dependencies" OFF)

# Set SFML_DIR before calling find_package only if we're not in tests-only mode
if(NOT BUILD_TESTS_ONLY)
  set(SFML_DIR "/opt/homebrew/lib/cmake/SFML")
  find_package(SFML 3.0 COMPONENTS Graphics Window System REQUIRED)
endif()

# Output directory setup
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directory setup
include_directories(include)

# Subdirectory setup
add_subdirectory(src)
add_subdirectory(tests)

# Add simulations only if not in tests-only mode
if(NOT BUILD_TESTS_ONLY)
  add_subdirectory(simulations)
endif()