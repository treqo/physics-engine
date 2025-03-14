#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>
#include <iomanip>  // for formatting output
#include <stdexcept> // for error handling
#include <cmath> // for determinant, etc.
#include <utility> // for std::pair

class Matrix
{
private:
    std::vector<std::vector<double>> data;
    size_t rows, cols;

public:
    /* Constructors */
    Matrix(); // Default constructor - creates an empty 0x0 matrix
    Matrix(size_t rows, size_t cols, double init_val = 0.0);
    Matrix(const std::vector<std::vector<double>>& values);
    Matrix(const Matrix& other); // Copy constructor
    Matrix(Matrix&& other) noexcept; // Move constructor

    /* Static factory methods */
    static Matrix identity(size_t size); // Creates an identity matrix
    static Matrix diagonal(const std::vector<double>& values); // Creates a diagonal matrix

    /* Destructor */
    ~Matrix();

    /* Assignment operators */
    Matrix& operator=(const Matrix& other); // Copy assignment
    Matrix& operator=(Matrix&& other) noexcept; // Move assignment

    /* Access elements */
    double& operator()(size_t row, size_t col);
    double operator()(size_t row, size_t col) const;

    /* Matrix Operations */
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator/(double scalar) const; // Division by scalar
    
    /* Compound Assignment Operators */
    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix& operator*=(double scalar);
    Matrix& operator*=(const Matrix& other);
    Matrix& operator/=(double scalar);

    /* Matrix mathematical operations */
    Matrix transpose() const;
    double determinant() const;  // Only for square matrices
    Matrix inverse() const;      // Only for square matrices
    Matrix power(int exponent) const; // Only for square matrices
    double trace() const; // Sum of diagonal elements
    int rank() const; // Rank of the matrix

    /* Decompositions */
    std::pair<Matrix, Matrix> luDecomposition() const; // Returns L and U matrices

    /* Comparison operators */
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;

    /* Utilities */
    void print() const;
    std::string toString() const;
    bool isSquare() const { return rows == cols; }
    bool isSymmetric() const;
    bool isDiagonal() const;
    bool isIdentity(double tolerance = 1e-10) const;
    bool isSingular(double tolerance = 1e-10) const;
    
    /* Getters */
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    std::vector<std::vector<double>> getData() const { return data; }

    /* Friend functions */
    friend Matrix operator*(double scalar, const Matrix& mat);
    friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);
};

/* Global operators */
Matrix operator*(double scalar, const Matrix& mat);
std::ostream& operator<<(std::ostream& os, const Matrix& mat);

#endif // MATRIX_HPP