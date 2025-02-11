#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cmath>
#include <iostream>

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

	/* --- Vector Operators --- */

	/* Arithmetic Operators */
	Vector2D operator+(Vector2D const& v) const;
	Vector2D operator-(Vector2D const& v) const;
	Vector2D operator*(double a) const;
	Vector2D operator/(double a) const;

	/* Compound Assignment Operators */
	Vector2D& operator+=(Vector2D const& v); 
	Vector2D& operator+=(double a);
	Vector2D& operator-=(Vector2D const& v);
	Vector2D& operator-=(double a);
	Vector2D& operator*=(double a);
	Vector2D& operator/=(double a);

	/* Comparison Operators */
	bool operator==(const Vector2D& v) const;
	bool operator!=(const Vector2D& v) const;

	/* Stream Operator */
	friend std::ostream& operator<<(std::ostream& os, const Vector2D&v);
};

#endif // VECTOR2D_HPP
