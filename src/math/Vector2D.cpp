#include "math/Vector2D.hpp"
#include <stdexcept>

Vector2D::Vector2D() : x(0.0), y(0.0) {}

Vector2D::Vector2D(double x_, double y_) : x(x_), y(y_) {}

double Vector2D::getX() const noexcept { return x; }
double Vector2D::getY() const noexcept { return y; }
