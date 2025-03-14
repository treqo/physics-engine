cmake_minimum_required(VERSION 3.10)
project(PhysicsEngine)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set SFML_DIR before calling find_package
set(SFML_DIR "/opt/homebrew/lib/cmake/SFML")

find_package(SFML 3.0 COMPONENTS Graphics Window System REQUIRED)

# Output directory setup
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directory setup
include_directories(include)

# Subdirectory setup
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(simulations)

