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
	Vector2D(const Vector2D& v);

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

	/* Vector Operations */
	double magnitude() const noexcept;

	void normalize();

	Vector2D unit() const;
	friend Vector2D unit(const Vector2D& v);

	double dot(const Vector2D& v) const;
	friend double dot(const Vector2D& v1, const Vector2D& v2);
	
	double angle(const Vector2D &v) const;

	double distanceTo(const Vector2D& v) const;
	friend double distanceBetween(const Vector2D& v1, const Vector2D& v2);

	Vector2D projectionOnto(const Vector2D& v) const;

	// Vector2D perpendicular() const;

	/* Additional Vector Operations */
	double cross(const Vector2D& v) const;
	friend double cross(const Vector2D& v1, const Vector2D& v2);

	Vector2D reflect(const Vector2D& normal) const;
	
	Vector2D rotate(double angle) const;
	
	static Vector2D lerp(const Vector2D& start, const Vector2D& end, double t);
	
	bool isZero(double epsilon = 1e-10) const noexcept;
	
	double squaredMagnitude() const noexcept;
};

#endif // VECTOR2D_HPP
