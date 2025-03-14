# simulations/CMakeLists.txt
if(NOT BUILD_TESTS_ONLY)
  add_executable(SFMLSim main.cpp)
  target_link_libraries(SFMLSim PRIVATE SFML::Graphics SFML::Window SFML::System physics_engine)
endif()