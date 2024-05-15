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

int Matrix::getCols() const {
    return this->col;
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

Matrix Matrix::transpose() const {
    Matrix result(col, fil);


    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[j][i] = matrix[i][j];
        }
    }

    return result;
}


Matrix Matrix::operator+(double scalar) {
    Matrix result(fil, col);

    for (int i = 0; i < fil; ++i) {
        for (int j = 0; j < col; ++j) {
            result.matrix[i][j] = matrix[i][j] + scalar;
        }
    }

    return result;
}

Matrix Matrix::inverse() const {
    // Verificar si la matriz es cuadrada
    if (fil != col) {
        std::cerr << "La matriz no es cuadrada, no se puede calcular la inversa." << std::endl;
        exit(EXIT_FAILURE);
    }

    int size = fil;
    Matrix identity = createIdentityMatrix(size);
    Matrix augmentedMatrix(size, size * 2);

    // Construir la matriz aumentada [A | I]
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            augmentedMatrix(i + 1, j + 1) = matrix[i][j];
            augmentedMatrix(i + 1, j + 1 + size) = identity(i + 1, j + 1);
        }
    }

    // Aplicar eliminación de Gauss-Jordan
    for (int i = 1; i <= size; ++i) {
        // Dividir la fila i por el pivote
        double pivot = augmentedMatrix(i, i);
        for (int j = 1; j <= size * 2; ++j) {
            augmentedMatrix(i, j) /= pivot;
        }
        // Hacer ceros en las demás filas debajo del pivote
        for (int k = 1; k <= size; ++k) {
            if (k != i) {
                double factor = augmentedMatrix(k, i);
                for (int j = 1; j <= size * 2; ++j) {
                    augmentedMatrix(k, j) -= factor * augmentedMatrix(i, j);
                }
            }
        }
    }

    // Extraer la matriz inversa de la parte derecha de la matriz aumentada
    Matrix inverse(size, size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            inverse(i + 1, j + 1) = augmentedMatrix(i + 1, j + 1 + size);
        }
    }

    return inverse;
}




Matrix Matrix::createIdentityMatrix(int size) {
    Matrix identity(size, size);

    for (int i = 1; i <= size; ++i) {
        for (int j = 1; j <= size; ++j) {
            if (i == j) {
                identity(i, j) = 1.0;
            } else {
                identity(i, j) = 0.0;
            }
        }
    }

    return identity;
}

Matrix Matrix::subMatrix(int row) {
    if (row < 1 || row > fil) {
        std::cerr << "Índice de fila fuera de rango." << std::endl;
        exit(EXIT_FAILURE);
    }

    Matrix sub(1, col); // Crear una nueva matriz con una sola fila

    for (int j = 0; j < col; ++j) {
        sub(1, j + 1) = matrix[row - 1][j];
    }

    return sub;
}


Matrix Matrix::range(int start, int step, int end) {
    if ((end - start) % step != 0) {
        std::cerr << "El rango no es divisible uniformemente por el paso." << std::endl;
        exit(EXIT_FAILURE);
    }

    int numElements = ((end - start) / step) + 1; // Calcular el número de elementos en el rango
    Matrix result(1, numElements);

    int value = start;
    for (int i = 0; i < numElements; ++i) {
        result(1, i + 1) = value;
        value += step;
    }

    return result;
}

Matrix Matrix::subMatrix(int row, int startCol, int endCol) {
    if (row < 1 || row > fil || startCol < 1 || startCol > col || endCol < 1 || endCol > col || startCol > endCol) {
        std::cerr << "Índices de fila o columna fuera de rango." << std::endl;
        exit(EXIT_FAILURE);
    }

    Matrix sub(1, endCol - startCol + 1); // Crear una nueva matriz con una sola fila y el rango de columnas especificado

    for (int j = startCol - 1, subCol = 0; j < endCol; ++j, ++subCol) {
        sub(1, subCol + 1) = matrix[row - 1][j];
    }

    return sub;
}


Matrix Matrix::concatenateHorizontal(const Matrix& matrix1, const Matrix& matrix2) {
    if (matrix1.getRows() != matrix2.getRows()) {
        std::cerr << "Las matrices tienen diferente número de filas." << std::endl;
        exit(EXIT_FAILURE);
    }

    int newCols = matrix1.getCols() + matrix2.getCols();
    Matrix result(matrix1.getRows(), newCols);

    for (int i = 1; i <= matrix1.getRows(); ++i) {
        for (int j = 1; j <= matrix1.getCols(); ++j) {
            result(i, j) = matrix1(i, j);
        }
        for (int j = 1; j <= matrix2.getCols(); ++j) {
            result(i, matrix1.getCols() + j) = matrix2(i, j);
        }
    }

    return result;
}





