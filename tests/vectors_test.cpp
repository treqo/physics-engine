#include <gtest/gtest.h>
#include "math/Vector2D.hpp"
#include <cmath>

// Utility function to compare floating point values
bool isEqual(double a, double b, double epsilon = 1e-10) {
    return std::abs(a - b) < epsilon;
}

// Test fixture for Vector2D tests
class Vector2DTest : public ::testing::Test {
protected:
    // Test vectors that will be used across multiple tests
    Vector2D zero;           // (0, 0)
    Vector2D unit_x;         // (1, 0)
    Vector2D unit_y;         // (0, 1)
    Vector2D one;            // (1, 1)
    Vector2D negative;       // (-1, -2)
    Vector2D arbitrary;      // (3.5, -2.7)

    void SetUp() override {
        // Initialize test vectors
        unit_x = Vector2D(1.0, 0.0);
        unit_y = Vector2D(0.0, 1.0);
        one = Vector2D(1.0, 1.0);
        negative = Vector2D(-1.0, -2.0);
        arbitrary = Vector2D(3.5, -2.7);
    }
};

// Test constructors and accessors
TEST_F(Vector2DTest, Construction) {
    // Default constructor
    EXPECT_DOUBLE_EQ(zero.getX(), 0.0);
    EXPECT_DOUBLE_EQ(zero.getY(), 0.0);
    
    // Value constructor
    Vector2D v(2.5, -3.7);
    EXPECT_DOUBLE_EQ(v.getX(), 2.5);
    EXPECT_DOUBLE_EQ(v.getY(), -3.7);
    
    // Copy constructor
    Vector2D v_copy(v);
    EXPECT_DOUBLE_EQ(v_copy.getX(), 2.5);
    EXPECT_DOUBLE_EQ(v_copy.getY(), -3.7);
    
    // Setters
    v.setX(1.1);
    v.setY(2.2);
    EXPECT_DOUBLE_EQ(v.getX(), 1.1);
    EXPECT_DOUBLE_EQ(v.getY(), 2.2);
}

// Test basic arithmetic operators
TEST_F(Vector2DTest, ArithmeticOperators) {
    // Addition
    Vector2D sum = unit_x + unit_y;
    EXPECT_DOUBLE_EQ(sum.getX(), 1.0);
    EXPECT_DOUBLE_EQ(sum.getY(), 1.0);
    
    // Subtraction
    Vector2D diff = arbitrary - one;
    EXPECT_DOUBLE_EQ(diff.getX(), 2.5);
    EXPECT_DOUBLE_EQ(diff.getY(), -3.7);
    
    // Scalar multiplication
    Vector2D scaled = arbitrary * 2.0;
    EXPECT_DOUBLE_EQ(scaled.getX(), 7.0);
    EXPECT_DOUBLE_EQ(scaled.getY(), -5.4);
    
    // Scalar division
    Vector2D divided = arbitrary / 2.0;
    EXPECT_DOUBLE_EQ(divided.getX(), 1.75);
    EXPECT_DOUBLE_EQ(divided.getY(), -1.35);
}

// Test compound assignment operators
TEST_F(Vector2DTest, CompoundAssignmentOperators) {
    // +=
    Vector2D v1 = arbitrary;
    v1 += one;
    EXPECT_DOUBLE_EQ(v1.getX(), 4.5);
    EXPECT_DOUBLE_EQ(v1.getY(), -1.7);
    
    Vector2D v1_scalar = arbitrary;
    v1_scalar += 1.0;
    EXPECT_DOUBLE_EQ(v1_scalar.getX(), 4.5);
    EXPECT_DOUBLE_EQ(v1_scalar.getY(), -1.7);
    
    // -=
    Vector2D v2 = arbitrary;
    v2 -= one;
    EXPECT_DOUBLE_EQ(v2.getX(), 2.5);
    EXPECT_DOUBLE_EQ(v2.getY(), -3.7);
    
    Vector2D v2_scalar = arbitrary;
    v2_scalar -= 1.0;
    EXPECT_DOUBLE_EQ(v2_scalar.getX(), 2.5);
    EXPECT_DOUBLE_EQ(v2_scalar.getY(), -3.7);
    
    // *=
    Vector2D v3 = arbitrary;
    v3 *= 2.0;
    EXPECT_DOUBLE_EQ(v3.getX(), 7.0);
    EXPECT_DOUBLE_EQ(v3.getY(), -5.4);
    
    // /=
    Vector2D v4 = arbitrary;
    v4 /= 2.0;
    EXPECT_DOUBLE_EQ(v4.getX(), 1.75);
    EXPECT_DOUBLE_EQ(v4.getY(), -1.35);
}

// Test comparison operators
TEST_F(Vector2DTest, ComparisonOperators) {
    // ==
    EXPECT_TRUE(zero == Vector2D(0.0, 0.0));
    EXPECT_TRUE(arbitrary == Vector2D(3.5, -2.7));
    EXPECT_FALSE(zero == one);
    
    // !=
    EXPECT_TRUE(zero != one);
    EXPECT_FALSE(arbitrary != Vector2D(3.5, -2.7));
}

// Test magnitude calculations
TEST_F(Vector2DTest, Magnitude) {
    // magnitude()
    EXPECT_DOUBLE_EQ(zero.magnitude(), 0.0);
    EXPECT_DOUBLE_EQ(unit_x.magnitude(), 1.0);
    EXPECT_DOUBLE_EQ(unit_y.magnitude(), 1.0);
    EXPECT_NEAR(one.magnitude(), std::sqrt(2.0), 1e-10);
    
    // squaredMagnitude()
    EXPECT_DOUBLE_EQ(zero.squaredMagnitude(), 0.0);
    EXPECT_DOUBLE_EQ(unit_x.squaredMagnitude(), 1.0);
    EXPECT_DOUBLE_EQ(one.squaredMagnitude(), 2.0);
    EXPECT_DOUBLE_EQ(arbitrary.squaredMagnitude(), 3.5*3.5 + (-2.7)*(-2.7));
}

// Test normalization
TEST_F(Vector2DTest, Normalization) {
    // normalize() - Note: There's a bug in your implementation!
    Vector2D v = Vector2D(3.0, 4.0); // 3-4-5 triangle
    Vector2D v_copy = v;
    v_copy.normalize();
    EXPECT_NEAR(v_copy.getX(), 0.6, 1e-10);
    EXPECT_NEAR(v_copy.getY(), 0.8, 1e-10);
    EXPECT_NEAR(v_copy.magnitude(), 1.0, 1e-10);
    
    // unit()
    Vector2D v_unit = v.unit();
    EXPECT_NEAR(v_unit.getX(), 0.6, 1e-10);
    EXPECT_NEAR(v_unit.getY(), 0.8, 1e-10);
    EXPECT_NEAR(v_unit.magnitude(), 1.0, 1e-10);
    
    // unit() friend function
    Vector2D v_unit_friend = unit(v);
    EXPECT_NEAR(v_unit_friend.getX(), 0.6, 1e-10);
    EXPECT_NEAR(v_unit_friend.getY(), 0.8, 1e-10);
    EXPECT_NEAR(v_unit_friend.magnitude(), 1.0, 1e-10);
}

// Test dot product
TEST_F(Vector2DTest, DotProduct) {
    // dot() method
    EXPECT_DOUBLE_EQ(unit_x.dot(unit_y), 0.0);
    EXPECT_DOUBLE_EQ(unit_x.dot(unit_x), 1.0);
    EXPECT_DOUBLE_EQ(one.dot(one), 2.0);
    
    // dot() friend function
    EXPECT_DOUBLE_EQ(dot(unit_x, unit_y), 0.0);
    EXPECT_DOUBLE_EQ(dot(arbitrary, one), 3.5 - 2.7);
}

// Test angle calculation
TEST_F(Vector2DTest, Angle) {
    // angle() method - returns angle in radians
    EXPECT_NEAR(unit_x.angle(unit_y), M_PI/2, 1e-10);
    EXPECT_NEAR(unit_x.angle(unit_x), 0.0, 1e-10);
    
    // 45-degree test
    Vector2D v45(1.0, 1.0);
    EXPECT_NEAR(unit_x.angle(v45), M_PI/4, 1e-10);
}

// Test distance calculations
TEST_F(Vector2DTest, Distance) {
    // distanceTo() method
    EXPECT_DOUBLE_EQ(zero.distanceTo(unit_x), 1.0);
    EXPECT_DOUBLE_EQ(unit_x.distanceTo(unit_y), std::sqrt(2.0));
    
    // distanceBetween() friend function
    EXPECT_DOUBLE_EQ(distanceBetween(zero, unit_x), 1.0);
    EXPECT_DOUBLE_EQ(distanceBetween(arbitrary, one), std::sqrt(std::pow(3.5-1.0, 2) + std::pow(-2.7-1.0, 2)));
}

// Test projection
TEST_F(Vector2DTest, Projection) {
    // projectionOnto() method
    Vector2D proj = one.projectionOnto(unit_x);
    EXPECT_NEAR(proj.getX(), 1.0, 1e-10);
    EXPECT_NEAR(proj.getY(), 0.0, 1e-10);
    
    // Project (3,4) onto (1,0) should give (3,0)
    Vector2D v(3.0, 4.0);
    Vector2D proj2 = v.projectionOnto(unit_x);
    EXPECT_NEAR(proj2.getX(), 3.0, 1e-10);
    EXPECT_NEAR(proj2.getY(), 0.0, 1e-10);
}

// Test cross product (in 2D, this is a scalar)
TEST_F(Vector2DTest, CrossProduct) {
    // cross() method
    EXPECT_DOUBLE_EQ(unit_x.cross(unit_y), 1.0);
    EXPECT_DOUBLE_EQ(unit_y.cross(unit_x), -1.0);
    EXPECT_DOUBLE_EQ(unit_x.cross(unit_x), 0.0);
    
    // cross() friend function
    EXPECT_DOUBLE_EQ(cross(unit_x, unit_y), 1.0);
    EXPECT_DOUBLE_EQ(cross(arbitrary, one), 3.5*1.0 - (-2.7)*1.0);
}

// Test reflection
TEST_F(Vector2DTest, Reflection) {
    // reflect() method
    Vector2D v(1.0, -1.0);
    Vector2D normal(0.0, 1.0); // Unit vector pointing up
    
    Vector2D reflected = v.reflect(normal);
    EXPECT_NEAR(reflected.getX(), 1.0, 1e-10);
    EXPECT_NEAR(reflected.getY(), 1.0, 1e-10);
    
    // Test against a horizontal normal
    Vector2D horizontal_normal(1.0, 0.0);
    Vector2D v2(1.0, 1.0);
    Vector2D reflected2 = v2.reflect(horizontal_normal);
    EXPECT_NEAR(reflected2.getX(), -1.0, 1e-10);
    EXPECT_NEAR(reflected2.getY(), 1.0, 1e-10);
}

// Test rotation
TEST_F(Vector2DTest, Rotation) {
    // rotate() method
    Vector2D rotated = unit_x.rotate(M_PI/2); // 90 degrees
    EXPECT_NEAR(rotated.getX(), 0.0, 1e-10);
    EXPECT_NEAR(rotated.getY(), 1.0, 1e-10);
    
    // Full circle should return to original
    Vector2D full_circle = arbitrary.rotate(2 * M_PI);
    EXPECT_NEAR(full_circle.getX(), arbitrary.getX(), 1e-10);
    EXPECT_NEAR(full_circle.getY(), arbitrary.getY(), 1e-10);
}

// Test linear interpolation
TEST_F(Vector2DTest, LinearInterpolation) {
    // lerp() static method
    Vector2D start(0.0, 0.0);
    Vector2D end(10.0, 10.0);
    
    // t = 0 should be start
    Vector2D lerp0 = Vector2D::lerp(start, end, 0.0);
    EXPECT_DOUBLE_EQ(lerp0.getX(), 0.0);
    EXPECT_DOUBLE_EQ(lerp0.getY(), 0.0);
    
    // t = 1 should be end
    Vector2D lerp1 = Vector2D::lerp(start, end, 1.0);
    EXPECT_DOUBLE_EQ(lerp1.getX(), 10.0);
    EXPECT_DOUBLE_EQ(lerp1.getY(), 10.0);
    
    // t = 0.5 should be midpoint
    Vector2D lerp_mid = Vector2D::lerp(start, end, 0.5);
    EXPECT_DOUBLE_EQ(lerp_mid.getX(), 5.0);
    EXPECT_DOUBLE_EQ(lerp_mid.getY(), 5.0);
    
    // t > 1 should be clamped to 1
    Vector2D lerp_over = Vector2D::lerp(start, end, 2.0);
    EXPECT_DOUBLE_EQ(lerp_over.getX(), 10.0);
    EXPECT_DOUBLE_EQ(lerp_over.getY(), 10.0);
    
    // t < 0 should be clamped to 0
    Vector2D lerp_under = Vector2D::lerp(start, end, -1.0);
    EXPECT_DOUBLE_EQ(lerp_under.getX(), 0.0);
    EXPECT_DOUBLE_EQ(lerp_under.getY(), 0.0);
}

// Test isZero method
TEST_F(Vector2DTest, IsZero) {
    // isZero() method
    EXPECT_TRUE(zero.isZero());
    EXPECT_FALSE(one.isZero());
    
    // Test with custom epsilon
    Vector2D almost_zero(1e-6, 1e-6);
    EXPECT_FALSE(almost_zero.isZero(1e-10)); // Smaller epsilon, should fail
    EXPECT_TRUE(almost_zero.isZero(1e-5));   // Larger epsilon, should pass
}

// Test stream operator
TEST_F(Vector2DTest, StreamOperator) {
    // This test checks if the ostream operator produces expected output
    std::stringstream ss;
    ss << arbitrary;
    EXPECT_EQ(ss.str(), "(3.5,-2.7)");
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}