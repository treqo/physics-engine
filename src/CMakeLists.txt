add_library(physics_engine STATIC
    math/Vector2D.cpp
)

target_include_directories(physics_engine PUBLIC ${PROJECT_SOURCE_DIR}/include)

