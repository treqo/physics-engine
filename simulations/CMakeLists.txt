# simulations/CMakeLists.txt
add_executable(SFMLSim main.cpp)

target_link_libraries(SFMLSim PRIVATE SFML::Graphics SFML::Window SFML::System physics_engine)

