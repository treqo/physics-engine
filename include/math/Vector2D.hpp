#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cmath>

class Vector2D {
private:
	double x, y;
public:
	/* Constructors */
	Vector2D();
	Vector2D(double x_, double y_);

	/* Getters */
	double getX() const noexcept;
	double getY() const noexcept;

	/* Vector Operators */
};

#endif // VECTOR2D_HPP
