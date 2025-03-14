#include "math/Matrix.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <sstream>

// Default constructor
Matrix::Matrix() : rows(0), cols(0) {}

// Constructor with dimensions
Matrix::Matrix(size_t rows, size_t cols, double init_val)
    : rows(rows), cols(cols)
{
    if (rows == 0 || cols == 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    
    data.resize(rows, std::vector<double>(cols, init_val));
}

// Constructor from 2D vector
Matrix::Matrix(const std::vector<std::vector<double>>& values)
    : rows(values.size()), cols(values.empty() ? 0 : values[0].size())
{
    if (rows == 0 || cols == 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    
    // Verify that all rows have the same number of columns
    for (const auto& row : values) {
        if (row.size() != cols) {
            throw std::invalid_argument("All rows must have the same number of columns");
        }
    }
    
    data = values;
}

// Copy constructor
Matrix::Matrix(const Matrix& other)
    : rows(other.rows), cols(other.cols), data(other.data)
{
}

// Move constructor
Matrix::Matrix(Matrix&& other) noexcept
    : rows(other.rows), cols(other.cols), data(std::move(other.data))
{
    other.rows = 0;
    other.cols = 0;
}

// Identity matrix factory
Matrix Matrix::identity(size_t size)
{
    Matrix result(size, size, 0.0);
    for (size_t i = 0; i < size; ++i) {
        result(i, i) = 1.0;
    }
    return result;
}

// Diagonal matrix factory
Matrix Matrix::diagonal(const std::vector<double>& values)
{
    size_t size = values.size();
    Matrix result(size, size, 0.0);
    for (size_t i = 0; i < size; ++i) {
        result(i, i) = values[i];
    }
    return result;
}

// Destructor
Matrix::~Matrix() = default;

// Copy assignment
Matrix& Matrix::operator=(const Matrix& other)
{
    if (this != &other) {
        rows = other.rows;
        cols = other.cols;
        data = other.data;
    }
    return *this;
}

// Move assignment
Matrix& Matrix::operator=(Matrix&& other) noexcept {
    if (this != &other) {
        data = std::move(other.data);
        rows = other.rows;
        cols = other.cols;

        // Reset other to a valid empty state
        other.rows = 0;
        other.cols = 0;
        other.data.clear();
    }
    return *this;
}

// Element access
double& Matrix::operator()(size_t row, size_t col)
{
    if (row >= rows || col >= cols) {
        throw std::out_of_range("Matrix indices out of range");
    }
    return data[row][col];
}

double Matrix::operator()(size_t row, size_t col) const
{
    if (row >= rows || col >= cols) {
        throw std::out_of_range("Matrix indices out of range");
    }
    return data[row][col];
}

// Matrix addition
Matrix Matrix::operator+(const Matrix& other) const
{
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }
    
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }
    
    return result;
}

// Matrix subtraction
Matrix Matrix::operator-(const Matrix& other) const
{
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }
    
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) - other(i, j);
        }
    }
    
    return result;
}

// Scalar multiplication
Matrix Matrix::operator*(double scalar) const
{
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) * scalar;
        }
    }
    
    return result;
}

// Matrix multiplication
Matrix Matrix::operator*(const Matrix& other) const
{
    if (cols != other.rows) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
    }
    
    Matrix result(rows, other.cols, 0.0);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < other.cols; ++j) {
            for (size_t k = 0; k < cols; ++k) {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    
    return result;
}

// Scalar division
Matrix Matrix::operator/(double scalar) const
{
    if (std::abs(scalar) < 1e-10) {
        throw std::invalid_argument("Division by zero or near-zero");
    }
    
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(i, j) = (*this)(i, j) / scalar;
        }
    }
    
    return result;
}

// Compound addition
Matrix& Matrix::operator+=(const Matrix& other)
{
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            data[i][j] += other(i, j);
        }
    }
    
    return *this;
}

// Compound subtraction
Matrix& Matrix::operator-=(const Matrix& other)
{
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            data[i][j] -= other(i, j);
        }
    }
    
    return *this;
}

// Compound scalar multiplication
Matrix& Matrix::operator*=(double scalar)
{
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            data[i][j] *= scalar;
        }
    }
    
    return *this;
}

// Compound matrix multiplication
Matrix& Matrix::operator*=(const Matrix& other) {
    if (cols != other.rows) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
    }
    *this = std::move(*this * other);
    return *this;
}

// Compound scalar division
Matrix& Matrix::operator/=(double scalar)
{
    if (std::abs(scalar) < 1e-10) {
        throw std::invalid_argument("Division by zero or near-zero");
    }
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            data[i][j] /= scalar;
        }
    }
    
    return *this;
}

// Matrix transpose
Matrix Matrix::transpose() const
{
    Matrix result(cols, rows);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result(j, i) = (*this)(i, j);
        }
    }
    
    return result;
}

// Matrix determinant (only for square matrices)
double Matrix::determinant() const {
    if (!isSquare()) {
        throw std::invalid_argument("Determinant is defined only for square matrices");
    }

    auto [L, U] = this->luDecomposition();
    
    double det = 1.0;
    for (size_t i = 0; i < rows; ++i) {
        det *= U(i, i);  // Determinant is product of diagonal elements of U
    }
    
    return det;
}

// Matrix inverse (only for square matrices)
Matrix Matrix::inverse() const
{
    if (rows != cols) {
        throw std::invalid_argument("Inverse is defined only for square matrices");
    }
    
    double det = determinant();
    if (std::abs(det) < 1e-10) {
        throw std::runtime_error("Matrix is singular, cannot compute inverse");
    }
    
    // For 1x1 matrix
    if (rows == 1) {
        Matrix result(1, 1);
        result(0, 0) = 1.0 / data[0][0];
        return result;
    }
    
    // For 2x2 matrix, use the simple formula
    if (rows == 2) {
        Matrix result(2, 2);
        result(0, 0) = data[1][1] / det;
        result(0, 1) = -data[0][1] / det;
        result(1, 0) = -data[1][0] / det;
        result(1, 1) = data[0][0] / det;
        return result;
    }
    
    // For larger matrices, compute the adjugate matrix and divide by determinant
    Matrix result(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            // Create submatrix by removing row i and column j
            std::vector<std::vector<double>> submatrix_data(rows - 1, std::vector<double>(cols - 1));
            size_t row_idx = 0;
            for (size_t r = 0; r < rows; ++r) {
                if (r != i) {
                    size_t col_idx = 0;
                    for (size_t c = 0; c < cols; ++c) {
                        if (c != j) {
                            submatrix_data[row_idx][col_idx] = data[r][c];
                            ++col_idx;
                        }
                    }
                    ++row_idx;
                }
            }
            
            Matrix submatrix(submatrix_data);
            double minor = submatrix.determinant();
            double cofactor = minor * ((i + j) % 2 == 0 ? 1 : -1);
            
            // Adjugate is the transpose of the cofactor matrix
            result(j, i) = cofactor / det;
        }
    }
    
    return result;
}

// Matrix power (only for square matrices)
Matrix Matrix::power(int exponent) const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix power is defined only for square matrices");
    }
    
    if (exponent < 0) {
        return this->inverse().power(-exponent);
    }
    
    Matrix result = Matrix::identity(rows);
    Matrix base = *this;
    
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        base *= base;
        exponent /= 2;
    }
    
    return result;
}

// Matrix trace (sum of diagonal elements)
double Matrix::trace() const
{
    if (!isSquare()) {
        throw std::invalid_argument("Trace is defined only for square matrices");
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < rows; ++i) {
        sum += data[i][i];
    }
    
    return sum;
}

// Matrix rank (using Gaussian elimination)
int Matrix::rank() const
{
    // Create a copy of the matrix for row operations
    std::vector<std::vector<double>> temp_data = data;
    
    // Apply Gaussian elimination
    size_t r = 0; // Current row
    for (size_t c = 0; c < cols && r < rows; ++c) {
        // Find the pivot row (maximum absolute value in the current column)
        size_t pivot_row = r;
        double max_val = std::abs(temp_data[r][c]);
        
        for (size_t i = r + 1; i < rows; ++i) {
            if (std::abs(temp_data[i][c]) > max_val) {
                max_val = std::abs(temp_data[i][c]);
                pivot_row = i;
            }
        }
        
        // If the pivot is near zero, this column doesn't have a pivot
        if (max_val < 1e-10) {
            continue;
        }
        
        // Swap rows if needed
        if (pivot_row != r) {
            temp_data[r].swap(temp_data[pivot_row]);
        }
        
        // Eliminate other rows
        for (size_t i = 0; i < rows; ++i) {
            if (i != r) {
                double factor = temp_data[i][c] / temp_data[r][c];
                for (size_t j = c; j < cols; ++j) {
                    temp_data[i][j] -= factor * temp_data[r][j];
                }
            }
        }
        
        ++r; // Move to the next row
    }
    
    // The rank is the number of non-zero rows
    int rank = 0;
    for (size_t i = 0; i < rows; ++i) {
        bool is_zero_row = true;
        for (size_t j = 0; j < cols; ++j) {
            if (std::abs(temp_data[i][j]) >= 1e-10) {
                is_zero_row = false;
                break;
            }
        }
        if (!is_zero_row) {
            ++rank;
        }
    }
    
    return rank;
}

// LU Decomposition
std::pair<Matrix, Matrix> Matrix::luDecomposition() const
{
    if (!isSquare()) {
        throw std::invalid_argument("LU decomposition is defined only for square matrices");
    }
    
    // Initialize L and U matrices
    Matrix L(rows, cols, 0.0);
    Matrix U(rows, cols, 0.0);
    
    for (size_t i = 0; i < rows; ++i) {
        // Upper triangular matrix
        for (size_t j = i; j < cols; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L(i, k) * U(k, j);
            }
            U(i, j) = data[i][j] - sum;
        }
        
        // Lower triangular matrix
        L(i, i) = 1.0; // Diagonal elements of L are 1
        for (size_t j = i + 1; j < rows; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < i; ++k) {
                sum += L(j, k) * U(k, i);
            }
            if (std::abs(U(i, i)) < 1e-10) {
                throw std::runtime_error("LU decomposition failed: Singular matrix");
            }
            L(j, i) = (data[j][i] - sum) / U(i, i);
        }
    }
    
    return {L, U};
}

// Check if matrix is singular
bool Matrix::isSingular(double tolerance) const {
    return std::abs(this->determinant()) < tolerance;
}

// Equality comparison
bool Matrix::operator==(const Matrix& other) const
{
    if (rows != other.rows || cols != other.cols) {
        return false;
    }
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (data[i][j] != other.data[i][j]) {
                return false;
            }
        }
    }
    
    return true;
}

// Inequality comparison
bool Matrix::operator!=(const Matrix& other) const
{
    return !(*this == other);
}

// Print matrix
void Matrix::print() const {
    for (const auto& row : data) {
        std::cout << "| ";
        for (double val : row) {
            std::cout << std::setw(8) << std::setprecision(4) << val << " ";
        }
        std::cout << "|\n";
    }
}

// Convert matrix to string
std::string Matrix::toString() const
{
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < rows; ++i) {
        if (i > 0) {
            ss << " ";
        }
        ss << "[";
        for (size_t j = 0; j < cols; ++j) {
            ss << std::setprecision(4) << data[i][j];
            if (j < cols - 1) {
                ss << ", ";
            }
        }
        ss << "]";
        if (i < rows - 1) {
            ss << "," << std::endl;
        }
    }
    ss << "]";
    return ss.str();
}

// Check if matrix is symmetric
bool Matrix::isSymmetric() const {
    if (!isSquare()) return false;
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = i + 1; j < cols; ++j) {
            if (std::abs(data[i][j] - data[j][i]) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if matrix is diagonal
bool Matrix::isDiagonal() const
{
    if (!isSquare()) {
        return false;
    }
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (i != j && std::abs(data[i][j]) > 1e-10) {
                return false;
            }
        }
    }
    
    return true;
}

// Check if matrix is identity
bool Matrix::isIdentity(double tolerance) const
{
    if (!isSquare()) {
        return false;
    }
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(data[i][j] - expected) > tolerance) {
                return false;
            }
        }
    }
    
    return true;
}

// Global operator for scalar * matrix
Matrix operator*(double scalar, const Matrix& mat)
{
    return mat * scalar;
}

// Stream operator for Matrix
std::ostream& operator<<(std::ostream& os, const Matrix& mat)
{
    os << mat.toString();
    return os;
}