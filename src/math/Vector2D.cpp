#include "math/Vector2D.hpp"
#include <stdexcept>

Vector2D::Vector2D() : x(0.0), y(0.0) {}

Vector2D::Vector2D(double x_, double y_) : x(x_), y(y_) {}

double Vector2D::getX() const noexcept { return x; }
double Vector2D::getY() const noexcept { return y; }

Vector2D Vector2D::operator+(Vector2D const& v) const {
	return Vector2D(this->x + v.x, this->y + v.y);
}

Vector2D Vector2D::operator-(Vector2D const& v) const {
	return Vector2D(this->x - v.x, this->y - v.y);
}

Vector2D Vector2D::operator*(double a) const {
	return Vector2D(this->x * a, this->y * a);
}

Vector2D Vector2D::operator/(double a) const {
	return Vector2D(this->x / a, this->y / a);
}

Vector2D& Vector2D::operator+=(Vector2D const& v) {
	this->x += v.x;
	this->y += v.y;
	return *this;
}

Vector2D& Vector2D::operator+=(double a) {
	this->x += a;
	this->y += a;
	return *this;
}

Vector2D& Vector2D::operator-=(Vector2D const& v) {
	this->x -= v.x;
	this->y -= v.y;
	return *this;
}


Vector2D& Vector2D::operator-=(double a) {
	this->x -= a;
	this->y -= a;
	return *this;
}

Vector2D& Vector2D::operator*=(double a) {
	this->x *= a;
	this->y *= a;
	return *this;
}

Vector2D& Vector2D::operator/=(double a) {
	this->x /= a;
	this->y /= a;
	return *this;
}

bool Vector2D::operator==(const Vector2D& v) const {
	return (this->x == v.x) && (this->y == v.y);
}

bool Vector2D::operator!=(const Vector2D& v) const {
	return (this->x != v.x) || (this->y != v.y);
}

std::ostream& operator<<(std::ostream& os, const Vector2D& v) {
	os << "(" << v.x << "," << v.y << ")";
	return os;
}

double Vector2D::magnitude() const noexcept {
	return sqrt(pow(this->x, 2) + pow(this->y, 2));
}

Vector2D Vector2D::unit(Vector2D& v) const {
	return v / v.magnitude();
}

double Vector2D::dot(Vector2D& v) const {
	return v.x * this->x + v.y * this->y;
}


double dot(const Vector2D& v1, const Vector2D& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}




