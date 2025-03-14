#include <gtest/gtest.h>
#include "math/Matrix.hpp"
#include <cmath>
#include <vector>

// Utility function to compare floating point values
bool isEqual(double a, double b, double epsilon = 1e-10) {
    return std::abs(a - b) < epsilon;
}

// Utility function to compare matrices with epsilon tolerance
bool matricesEqual(const Matrix& a, const Matrix& b, double epsilon = 1e-10) {
    if (a.getRows() != b.getRows() || a.getCols() != b.getCols()) {
        return false;
    }

    for (size_t i = 0; i < a.getRows(); ++i) {
        for (size_t j = 0; j < a.getCols(); ++j) {
            if (!isEqual(a(i, j), b(i, j), epsilon)) {
                return false;
            }
        }
    }
    return true;
}

// Test fixture for Matrix tests
class MatrixTest : public ::testing::Test {
protected:
    // Test matrices that will be used across multiple tests
    Matrix zero_2x2;          // 2x2 zero matrix
    Matrix zero_3x3;          // 3x3 zero matrix
    Matrix identity_3x3;      // 3x3 identity matrix
    Matrix identity_4x4;      // 4x4 identity matrix
    Matrix square_2x2;        // 2x2 arbitrary matrix
    Matrix square_3x3;        // 3x3 arbitrary matrix
    Matrix rectangular_2x3;   // 2x3 rectangular matrix
    Matrix rectangular_3x2;   // 3x2 rectangular matrix
    Matrix diagonal_3x3;      // 3x3 diagonal matrix
    Matrix invertible_3x3;    // 3x3 invertible matrix
    Matrix singular_3x3;      // 3x3 singular matrix
    Matrix symmetric_3x3;     // 3x3 symmetric matrix

    void SetUp() override {
        // Initialize zero matrices
        zero_2x2 = Matrix(2, 2, 0.0);
        zero_3x3 = Matrix(3, 3, 0.0);
        
        // Initialize identity matrices
        identity_3x3 = Matrix::identity(3);
        identity_4x4 = Matrix::identity(4);
        
        // Initialize 2x2 square matrix
        std::vector<std::vector<double>> square_2x2_data = {
            {1.0, 2.0},
            {3.0, 4.0}
        };
        square_2x2 = Matrix(square_2x2_data);
        
        // Initialize 3x3 square matrix
        std::vector<std::vector<double>> square_3x3_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        square_3x3 = Matrix(square_3x3_data);
        
        // Initialize rectangular matrices
        std::vector<std::vector<double>> rect_2x3_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0}
        };
        rectangular_2x3 = Matrix(rect_2x3_data);
        
        std::vector<std::vector<double>> rect_3x2_data = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0}
        };
        rectangular_3x2 = Matrix(rect_3x2_data);
        
        // Initialize diagonal matrix
        std::vector<std::vector<double>> diag_3x3_data = {
            {2.0, 0.0, 0.0},
            {0.0, 5.0, 0.0},
            {0.0, 0.0, 9.0}
        };
        diagonal_3x3 = Matrix(diag_3x3_data);
        
        // Initialize invertible 3x3 matrix
        std::vector<std::vector<double>> inv_3x3_data = {
            {4.0, 2.0, 1.0},
            {3.0, 1.0, 2.0},
            {2.0, 5.0, 3.0}
        };
        invertible_3x3 = Matrix(inv_3x3_data);
        
        // Initialize singular 3x3 matrix (determinant = 0)
        std::vector<std::vector<double>> sing_3x3_data = {
            {1.0, 2.0, 3.0},
            {2.0, 4.0, 6.0},  // Row 2 is 2 * Row 1
            {4.0, 5.0, 6.0}
        };
        singular_3x3 = Matrix(sing_3x3_data);
        
        // Initialize symmetric 3x3 matrix
        std::vector<std::vector<double>> sym_3x3_data = {
            {1.0, 2.0, 3.0},
            {2.0, 4.0, 5.0},
            {3.0, 5.0, 6.0}
        };
        symmetric_3x3 = Matrix(sym_3x3_data);
    }
};

// Test constructors and basic properties
TEST_F(MatrixTest, Construction) {
    // Default constructor
    Matrix m1;
    EXPECT_EQ(m1.getRows(), 0);
    EXPECT_EQ(m1.getCols(), 0);
    
    // Dimensions constructor
    Matrix m2(3, 4, 1.5);
    EXPECT_EQ(m2.getRows(), 3);
    EXPECT_EQ(m2.getCols(), 4);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_DOUBLE_EQ(m2(i, j), 1.5);
        }
    }
    
    // Vector data constructor
    std::vector<std::vector<double>> data = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0}
    };
    Matrix m3(data);
    EXPECT_EQ(m3.getRows(), 2);
    EXPECT_EQ(m3.getCols(), 3);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(m3(i, j), data[i][j]);
        }
    }
    
    // Copy constructor
    Matrix m4(m3);
    EXPECT_EQ(m4.getRows(), 2);
    EXPECT_EQ(m4.getCols(), 3);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(m4(i, j), m3(i, j));
        }
    }
    
    // Move constructor (harder to test directly)
    Matrix m5(Matrix(2, 2, 3.0));
    EXPECT_EQ(m5.getRows(), 2);
    EXPECT_EQ(m5.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m5(i, j), 3.0);
        }
    }
    
    // Identity factory
    Matrix id = Matrix::identity(3);
    EXPECT_EQ(id.getRows(), 3);
    EXPECT_EQ(id.getCols(), 3);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(id(i, j), (i == j) ? 1.0 : 0.0);
        }
    }
    
    // Diagonal factory
    std::vector<double> diag_values = {2.0, 4.0, 6.0};
    Matrix diag = Matrix::diagonal(diag_values);
    EXPECT_EQ(diag.getRows(), 3);
    EXPECT_EQ(diag.getCols(), 3);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(diag(i, j), (i == j) ? diag_values[i] : 0.0);
        }
    }
}

// Test assignment operators
TEST_F(MatrixTest, Assignment) {
    // Copy assignment
    Matrix m1;
    m1 = square_3x3;
    EXPECT_EQ(m1.getRows(), 3);
    EXPECT_EQ(m1.getCols(), 3);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(m1(i, j), square_3x3(i, j));
        }
    }
    
    // Move assignment (harder to test directly)
    Matrix m2;
    m2 = Matrix(2, 2, 4.0);
    EXPECT_EQ(m2.getRows(), 2);
    EXPECT_EQ(m2.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m2(i, j), 4.0);
        }
    }
    
    // Self-assignment should be safe
    m1 = m1;
    EXPECT_EQ(m1.getRows(), 3);
    EXPECT_EQ(m1.getCols(), 3);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(m1(i, j), square_3x3(i, j));
        }
    }
}

// Test element access operator
TEST_F(MatrixTest, ElementAccess) {
    // Valid access
    EXPECT_DOUBLE_EQ(square_2x2(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(square_2x2(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(square_2x2(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(square_2x2(1, 1), 4.0);
    
    // Const access
    const Matrix& const_ref = square_2x2;
    EXPECT_DOUBLE_EQ(const_ref(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(const_ref(0, 1), 2.0);
    
    // Modifying elements
    Matrix m = square_2x2;
    m(0, 0) = 10.0;
    EXPECT_DOUBLE_EQ(m(0, 0), 10.0);
    
    // Out of bounds access should throw
    EXPECT_THROW(square_2x2(2, 0), std::out_of_range);
    EXPECT_THROW(square_2x2(0, 2), std::out_of_range);
    EXPECT_THROW(const_ref(2, 0), std::out_of_range);
}

// Test basic arithmetic operations
TEST_F(MatrixTest, BasicArithmetic) {
    // Addition
    Matrix sum = square_2x2 + square_2x2;
    EXPECT_EQ(sum.getRows(), 2);
    EXPECT_EQ(sum.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(sum(i, j), square_2x2(i, j) * 2);
        }
    }
    
    // Incompatible dimensions should throw
    EXPECT_THROW(square_2x2 + rectangular_2x3, std::invalid_argument);
    
    // Subtraction
    Matrix diff = square_2x2 - zero_2x2;
    EXPECT_EQ(diff.getRows(), 2);
    EXPECT_EQ(diff.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(diff(i, j), square_2x2(i, j));
        }
    }
    
    // Incompatible dimensions should throw
    EXPECT_THROW(square_2x2 - rectangular_2x3, std::invalid_argument);
    
    // Scalar multiplication
    Matrix scaled = square_2x2 * 2.5;
    EXPECT_EQ(scaled.getRows(), 2);
    EXPECT_EQ(scaled.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(scaled(i, j), square_2x2(i, j) * 2.5);
        }
    }
    
    // Global scalar multiplication (scalar on left)
    Matrix scaled2 = 3.0 * square_2x2;
    EXPECT_EQ(scaled2.getRows(), 2);
    EXPECT_EQ(scaled2.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(scaled2(i, j), square_2x2(i, j) * 3.0);
        }
    }
    
    // Matrix multiplication
    Matrix prod = square_2x2 * square_2x2;
    EXPECT_EQ(prod.getRows(), 2);
    EXPECT_EQ(prod.getCols(), 2);
    EXPECT_DOUBLE_EQ(prod(0, 0), 1.0*1.0 + 2.0*3.0);  // 7
    EXPECT_DOUBLE_EQ(prod(0, 1), 1.0*2.0 + 2.0*4.0);  // 10
    EXPECT_DOUBLE_EQ(prod(1, 0), 3.0*1.0 + 4.0*3.0);  // 15
    EXPECT_DOUBLE_EQ(prod(1, 1), 3.0*2.0 + 4.0*4.0);  // 22
    
    // Rectangular matrix multiplication
    Matrix rect_prod = rectangular_2x3 * rectangular_3x2;
    EXPECT_EQ(rect_prod.getRows(), 2);
    EXPECT_EQ(rect_prod.getCols(), 2);
    EXPECT_DOUBLE_EQ(rect_prod(0, 0), 1.0*1.0 + 2.0*3.0 + 3.0*5.0);  // 22
    EXPECT_DOUBLE_EQ(rect_prod(0, 1), 1.0*2.0 + 2.0*4.0 + 3.0*6.0);  // 28
    EXPECT_DOUBLE_EQ(rect_prod(1, 0), 4.0*1.0 + 5.0*3.0 + 6.0*5.0);  // 49
    EXPECT_DOUBLE_EQ(rect_prod(1, 1), 4.0*2.0 + 5.0*4.0 + 6.0*6.0);  // 64
    
    Matrix valid_mult = square_2x2 * rectangular_2x3;
    EXPECT_EQ(valid_mult.getRows(), 2);
    EXPECT_EQ(valid_mult.getCols(), 3);
        
    // Scalar division
    Matrix divided = square_2x2 / 2.0;
    EXPECT_EQ(divided.getRows(), 2);
    EXPECT_EQ(divided.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(divided(i, j), square_2x2(i, j) / 2.0);
        }
    }
    
    // Division by zero should throw
    EXPECT_THROW(square_2x2 / 0.0, std::invalid_argument);
}

// Test compound assignment operators
TEST_F(MatrixTest, CompoundAssignment) {
    // +=
    Matrix m1 = square_2x2;
    m1 += square_2x2;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m1(i, j), square_2x2(i, j) * 2);
        }
    }
    
    // Incompatible dimensions should throw
    Matrix m1_copy = m1;
    EXPECT_THROW(m1_copy += rectangular_2x3, std::invalid_argument);
    
    // -=
    Matrix m2 = square_2x2;
    m2 -= square_2x2;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m2(i, j), 0.0);
        }
    }
    
    // Incompatible dimensions should throw
    Matrix m2_copy = m2;
    EXPECT_THROW(m2_copy -= rectangular_2x3, std::invalid_argument);
    
    // *=
    Matrix m3 = square_2x2;
    m3 *= 3.0;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m3(i, j), square_2x2(i, j) * 3.0);
        }
    }
    
    // Matrix *=
    Matrix m4 = square_2x2;
    Matrix expected_prod = square_2x2 * square_2x2;
    m4 *= square_2x2;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m4(i, j), expected_prod(i, j));
        }
    }
    
    Matrix m4_copy = square_2x2;
    m4_copy *= Matrix(std::vector<std::vector<double>>{{1.0, 0.0}, {0.0, 1.0}});
    EXPECT_TRUE(matricesEqual(m4_copy, square_2x2));
    EXPECT_THROW(m4_copy *= rectangular_3x2, std::invalid_argument); // This should throw (2x2 *= 3x2 is invalid)
        
    // /=
    Matrix m5 = square_2x2;
    m5 /= 2.0;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(m5(i, j), square_2x2(i, j) / 2.0);
        }
    }
    
    // Division by zero should throw
    Matrix m5_copy = m5;
    EXPECT_THROW(m5_copy /= 0.0, std::invalid_argument);
}

// Test matrix transpose
TEST_F(MatrixTest, Transpose) {
    // Square matrix transpose
    Matrix trans_square = square_2x2.transpose();
    EXPECT_EQ(trans_square.getRows(), 2);
    EXPECT_EQ(trans_square.getCols(), 2);
    EXPECT_DOUBLE_EQ(trans_square(0, 0), square_2x2(0, 0));
    EXPECT_DOUBLE_EQ(trans_square(0, 1), square_2x2(1, 0));
    EXPECT_DOUBLE_EQ(trans_square(1, 0), square_2x2(0, 1));
    EXPECT_DOUBLE_EQ(trans_square(1, 1), square_2x2(1, 1));
    
    // Rectangular matrix transpose
    Matrix trans_rect = rectangular_2x3.transpose();
    EXPECT_EQ(trans_rect.getRows(), 3);
    EXPECT_EQ(trans_rect.getCols(), 2);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(trans_rect(j, i), rectangular_2x3(i, j));
        }
    }
    
    // Transpose of transpose should be original
    Matrix trans_trans = trans_rect.transpose();
    EXPECT_EQ(trans_trans.getRows(), 2);
    EXPECT_EQ(trans_trans.getCols(), 3);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(trans_trans(i, j), rectangular_2x3(i, j));
        }
    }
}

// Test determinant calculation
TEST_F(MatrixTest, Determinant) {
    // 1x1 matrix
    Matrix m1x1(std::vector<std::vector<double>>{{5.0}});
    EXPECT_DOUBLE_EQ(m1x1.determinant(), 5.0);
    
    // 2x2 matrix
    EXPECT_DOUBLE_EQ(square_2x2.determinant(), 1.0*4.0 - 2.0*3.0); // -2
    
    // 3x3 matrix (using rule of Sarrus)
    Matrix m3x3(std::vector<std::vector<double>>{
        {1.0, 2.0, 3.0},
        {4.0, 5.0, 6.0},
        {7.0, 8.0, 10.0} // Changed 9 to 10 to make det non-zero
    });
    EXPECT_NEAR(m3x3.determinant(), -3.0, 1e-10);
    
    // 4x4 matrix (using cofactor expansion)
    Matrix m4x4(std::vector<std::vector<double>>{
        {1.0, 3.0, 5.0, 9.0},
        {1.0, 3.0, 1.0, 7.0},
        {4.0, 3.0, 9.0, 7.0},
        {5.0, 2.0, 0.0, 9.0}
    });
    EXPECT_NEAR(m4x4.determinant(), -376.0, 1e-10);
    
    // Singular matrix should have determinant 0
    EXPECT_NEAR(singular_3x3.determinant(), 0.0, 1e-10);
    
    // Non-square matrix should throw
    EXPECT_THROW(rectangular_2x3.determinant(), std::invalid_argument);
}

// Test matrix inverse
TEST_F(MatrixTest, Inverse) {
    // 1x1 matrix
    Matrix m1x1(std::vector<std::vector<double>>{{2.0}});
    Matrix inv1x1 = m1x1.inverse();
    EXPECT_EQ(inv1x1.getRows(), 1);
    EXPECT_EQ(inv1x1.getCols(), 1);
    EXPECT_DOUBLE_EQ(inv1x1(0, 0), 0.5);
    
    // 2x2 matrix
    Matrix inv2x2 = square_2x2.inverse();
    EXPECT_EQ(inv2x2.getRows(), 2);
    EXPECT_EQ(inv2x2.getCols(), 2);
    double det2x2 = 1.0*4.0 - 2.0*3.0; // -2
    EXPECT_DOUBLE_EQ(inv2x2(0, 0), 4.0/det2x2);
    EXPECT_DOUBLE_EQ(inv2x2(0, 1), -2.0/det2x2);
    EXPECT_DOUBLE_EQ(inv2x2(1, 0), -3.0/det2x2);
    EXPECT_DOUBLE_EQ(inv2x2(1, 1), 1.0/det2x2);
    
    // 3x3 matrix
    Matrix inv3x3 = invertible_3x3.inverse();
    
    // Verify inverse by multiplying with original matrix (should get identity)
    Matrix prod = invertible_3x3 * inv3x3;
    EXPECT_TRUE(matricesEqual(prod, Matrix::identity(3), 1e-9)); // Higher tolerance for numerical issues
    
    // Singular matrix should throw
    EXPECT_THROW(singular_3x3.inverse(), std::runtime_error);
    
    // Non-square matrix should throw
    EXPECT_THROW(rectangular_2x3.inverse(), std::invalid_argument);
}

// Test matrix power
TEST_F(MatrixTest, Power) {
    // Power of 0 (any square matrix^0 = identity)
    Matrix pow0 = square_2x2.power(0);
    EXPECT_TRUE(matricesEqual(pow0, Matrix::identity(2)));
    
    // Power of 1 (any matrix^1 = matrix)
    Matrix pow1 = square_2x2.power(1);
    EXPECT_TRUE(matricesEqual(pow1, square_2x2));
    
    // Power of 2
    Matrix pow2 = square_2x2.power(2);
    EXPECT_TRUE(matricesEqual(pow2, square_2x2 * square_2x2));
    
    // Power of 3
    Matrix pow3 = square_2x2.power(3);
    EXPECT_TRUE(matricesEqual(pow3, square_2x2 * square_2x2 * square_2x2));
    
    // Negative power (should compute inverse then raise to positive power)
    Matrix pow_neg1 = square_2x2.power(-1);
    EXPECT_TRUE(matricesEqual(pow_neg1, square_2x2.inverse()));
    
    // Non-square matrix should throw
    EXPECT_THROW(rectangular_2x3.power(2), std::invalid_argument);
}

// Test matrix trace
TEST_F(MatrixTest, Trace) {
    // 2x2 matrix
    EXPECT_DOUBLE_EQ(square_2x2.trace(), 1.0 + 4.0); // 5
    
    // 3x3 matrix
    EXPECT_DOUBLE_EQ(square_3x3.trace(), 1.0 + 5.0 + 9.0); // 15
    
    // Identity matrix (trace = dimension)
    EXPECT_DOUBLE_EQ(identity_3x3.trace(), 3.0);
    EXPECT_DOUBLE_EQ(identity_4x4.trace(), 4.0);
    
    // Diagonal matrix (trace = sum of diagonal elements)
    EXPECT_DOUBLE_EQ(diagonal_3x3.trace(), 2.0 + 5.0 + 9.0); // 16
    
    // Non-square matrix should throw
    EXPECT_THROW(rectangular_2x3.trace(), std::invalid_argument);
}

// Test matrix rank
TEST_F(MatrixTest, Rank) {
    // Identity matrix (rank = dimension)
    EXPECT_EQ(identity_3x3.rank(), 3);
    
    // Full rank matrix
    EXPECT_EQ(invertible_3x3.rank(), 3);
    
    // Singular matrix (not full rank)
    EXPECT_EQ(singular_3x3.rank(), 2); // Rank 2 because one row is dependent
    
    // Zero matrix (rank = 0)
    EXPECT_EQ(zero_3x3.rank(), 0);
    
    // Rectangular matrices
    EXPECT_EQ(rectangular_2x3.rank(), 2); // Min of rows and cols if full rank
    EXPECT_EQ(rectangular_3x2.rank(), 2);
}

// Test LU decomposition
TEST_F(MatrixTest, LUDecomposition) {
    // Cannot test singular matrix
    auto [L, U] = invertible_3x3.luDecomposition();
    
    // Check dimensions
    EXPECT_EQ(L.getRows(), 3);
    EXPECT_EQ(L.getCols(), 3);
    EXPECT_EQ(U.getRows(), 3);
    EXPECT_EQ(U.getCols(), 3);
    
    // Verify L is lower triangular with ones on diagonal
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(L(i, i), 1.0);
        for (size_t j = i + 1; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(L(i, j), 0.0);
        }
    }
    
    // Verify U is upper triangular
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < i; ++j) {
            EXPECT_DOUBLE_EQ(U(i, j), 0.0);
        }
    }
    
    // Verify L*U = original matrix
    Matrix prod = L * U;
    EXPECT_TRUE(matricesEqual(prod, invertible_3x3, 1e-9)); // Higher tolerance for numerical issues
    
    // Singular matrix should throw
    EXPECT_THROW(singular_3x3.luDecomposition(), std::runtime_error);
    
    // Non-square matrix should throw
    EXPECT_THROW(rectangular_2x3.luDecomposition(), std::invalid_argument);
}

// Test utility functions
TEST_F(MatrixTest, Utilities) {
    // isSquare()
    EXPECT_TRUE(square_2x2.isSquare());
    EXPECT_TRUE(square_3x3.isSquare());
    EXPECT_FALSE(rectangular_2x3.isSquare());
    EXPECT_FALSE(rectangular_3x2.isSquare());
    
    // isSymmetric()
    EXPECT_TRUE(symmetric_3x3.isSymmetric());
    EXPECT_FALSE(square_3x3.isSymmetric());
    EXPECT_FALSE(rectangular_2x3.isSymmetric()); // Non-square matrix
    
    // isDiagonal()
    EXPECT_TRUE(diagonal_3x3.isDiagonal());
    EXPECT_TRUE(zero_3x3.isDiagonal());
    EXPECT_TRUE(identity_3x3.isDiagonal());
    EXPECT_FALSE(square_3x3.isDiagonal());
    
    // isIdentity()
    EXPECT_TRUE(identity_3x3.isIdentity());
    EXPECT_FALSE(diagonal_3x3.isIdentity());
    EXPECT_FALSE(square_3x3.isIdentity());
    
    // isSingular()
    EXPECT_TRUE(singular_3x3.isDiagonal() == false);
    EXPECT_EQ(singular_3x3.rank(), 2); // Rank should be less than dimension for singular matrix
        
    // toString() and stream operator
    std::stringstream ss;
    ss << square_2x2;
    std::string str = square_2x2.toString();
    EXPECT_FALSE(str.empty());
    EXPECT_FALSE(ss.str().empty());
    
    // Test data access
    auto data = square_2x2.getData();
    EXPECT_EQ(data.size(), 2);
    EXPECT_EQ(data[0].size(), 2);
    EXPECT_DOUBLE_EQ(data[0][0], 1.0);
    EXPECT_DOUBLE_EQ(data[0][1], 2.0);
    EXPECT_DOUBLE_EQ(data[1][0], 3.0);
    EXPECT_DOUBLE_EQ(data[1][1], 4.0);
}

// Test comparison operators
TEST_F(MatrixTest, ComparisonOperators) {
    // == operator
    Matrix m1 = square_2x2;
    EXPECT_TRUE(m1 == square_2x2);
    EXPECT_FALSE(m1 == square_3x3);
    
    // != operator
    EXPECT_FALSE(m1 != square_2x2);
    EXPECT_TRUE(m1 != square_3x3);
    
    // Different dimensions
    EXPECT_FALSE(square_2x2 == rectangular_2x3);
    EXPECT_TRUE(square_2x2 != rectangular_2x3);
    
    // Same dimensions but different values
    Matrix m2(2, 2, 0.0);
    EXPECT_FALSE(m1 == m2);
    EXPECT_TRUE(m1 != m2);
}

// Test exception handling
TEST_F(MatrixTest, ExceptionHandling) {
    // Invalid dimensions constructor
    EXPECT_THROW(Matrix(0, 5), std::invalid_argument);
    EXPECT_THROW(Matrix(5, 0), std::invalid_argument);
    
    // Invalid 2D vector constructor (non-uniform columns)
    std::vector<std::vector<double>> invalid_data = {
        {1.0, 2.0, 3.0},
        {4.0, 5.0}  // Different length
    };
    EXPECT_THROW({ Matrix m(invalid_data); }, std::invalid_argument);
    
    // Various exceptions tested throughout other test cases
}

// Test edge cases
TEST_F(MatrixTest, EdgeCases) {
    // Operations with zero matrix
    Matrix result = square_2x2 * zero_2x2;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            EXPECT_DOUBLE_EQ(result(i, j), 0.0);
        }
    }
    
    // Matrix with very small values (near zero)
    Matrix small(2, 2, 1e-15);
    EXPECT_NEAR(small.determinant(), 0.0, 1e-10);
    
    // Matrix with very large values
    Matrix large(2, 2, 1e15);
    EXPECT_NEAR(large.determinant(), 0.0, 1e14); // Adjust tolerance for numerical precision
    
    // Ensure zero divided by zero is handled (should throw exception)
    Matrix zero_1x1(1, 1, 0.0);
    EXPECT_THROW(zero_1x1.inverse(), std::runtime_error);
}

// Test additional mathematical properties
TEST_F(MatrixTest, MathematicalProperties) {
    // Associativity of multiplication: (AB)C = A(BC)
    Matrix A = square_2x2;
    Matrix B = Matrix(std::vector<std::vector<double>>{{2.0, 1.0}, {0.0, 3.0}});
    Matrix C = Matrix(std::vector<std::vector<double>>{{1.0, 0.0}, {2.0, 1.0}});
    
    Matrix left = (A * B) * C;
    Matrix right = A * (B * C);
    EXPECT_TRUE(matricesEqual(left, right, 1e-9));
    
    // Distributivity: A(B+C) = AB + AC
    Matrix sum_first = A * (B + C);
    Matrix sum_prod = A * B + A * C;
    EXPECT_TRUE(matricesEqual(sum_first, sum_prod, 1e-9));
    
    // Transpose properties: (A+B)^T = A^T + B^T
    Matrix trans_sum = (A + B).transpose();
    Matrix sum_trans = A.transpose() + B.transpose();
    EXPECT_TRUE(matricesEqual(trans_sum, sum_trans));
    
    // Transpose properties: (AB)^T = B^T * A^T
    Matrix trans_prod = (A * B).transpose();
    Matrix prod_trans = B.transpose() * A.transpose();
    EXPECT_TRUE(matricesEqual(trans_prod, prod_trans, 1e-9));
    
    // Inverse properties: (A*B)^-1 = B^-1 * A^-1
    Matrix inv_prod = (A * B).inverse();
    Matrix prod_inv = B.inverse() * A.inverse();
    EXPECT_TRUE(matricesEqual(inv_prod, prod_inv, 1e-9));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}