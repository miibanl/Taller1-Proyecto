#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}
 
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}
 
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}
 
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}
 
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}
 
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];
 
    return *this;
}
 
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}
 
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}
 
 
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}
 
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

double Matrix::norm() const {
    double sumaCuadrados = 0.0;

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            sumaCuadrados += matrix[i][j] * matrix[i][j];
        }
    }

    return std::sqrt(sumaCuadrados);
}

int Matrix::getRows() const {
    return this->fil;
}

Matrix Matrix::operator/(double scalar) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[i][j] = matrix[i][j] / scalar;
        }
    }

    return result;
}

Matrix Matrix::operator*(double scalar) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[i][j] = matrix[i][j] * scalar;
        }
    }

    return result;
}

Matrix operator*(double scalar, const Matrix& matrix) {
    Matrix result(matrix.fil, matrix.col);

    for (int i = 0; i < matrix.fil; ++i) {
        for (int j = 0; j < matrix.col; ++j) {
            result.matrix[i][j] = scalar * matrix.matrix[i][j];
        }
    }

    return result;
}

Matrix Matrix::cross(const Matrix& a, const Matrix& b) {

    double x = a(1, 2) * b(1, 3) - a(1, 3) * b(1, 2);
    double y = a(1, 3) * b(1, 1) - a(1, 1) * b(1, 3);
    double z = a(1, 1) * b(1, 2) - a(1, 2) * b(1, 1);

    Matrix result(1, 3);
    result(1, 1) = x;
    result(1, 2) = y;
    result(1, 3) = z;

    return result;
}


double Matrix::dot(const Matrix& a, const Matrix& b) {
    double result = 0.0;

    for (int i = 1; i <= 3; ++i) {
        result += a(1, i) * b(1, i);
    }

    return result;
}



