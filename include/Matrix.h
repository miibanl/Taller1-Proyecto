#ifndef _MATRIX_
#define _MATRIX_

#include <cmath>

/*!
 * @file Matrix.h
 * @class Matrix
 * @brief Class representing a matrix
 */
class Matrix
{
    public:
        /*!
        * @brief Constructor for creating a matrix with specified dimensions
        * @param fil Number of rows
        * @param col Number of columns
        */
        Matrix(int fil, int col);

        /*!
        * @brief Constructor for creating a matrix with specified dimensions and initial values
        * @param fil Number of rows
        * @param col Number of columns
        * @param v Array of initial values
        * @param n Number of elements in the array v
        */
        Matrix(int fil, int col, double v[], int n);

        /*!
        * @brief Copy constructor
        * @param m Another Matrix object to copy from
        */
        Matrix(const Matrix& m);

        /*!
        * @brief Destructor for cleaning up memory
        */
        ~Matrix();

        /*!
        * @brief Assignment operator
        * @param matrix2 Another Matrix object
        * @return Reference to the updated Matrix object
        */
        Matrix& operator=(const Matrix& matrix2);

        /*!
        * @brief Addition operator
        * @param matrix2 Another Matrix object
        * @return Resultant Matrix object after addition
        */
        Matrix  operator+(const Matrix& matrix2);

        /*!
        * @brief Subtraction operator
        * @param matrix2 Another Matrix object
        * @return Resultant Matrix object after subtraction
        */
        Matrix  operator-(const Matrix& matrix2);

        /*!
        * @brief Multiplication operator
        * @param matrix2 Another Matrix object
        * @return Resultant Matrix object after multiplication
        */
        Matrix  operator*(const Matrix& matrix2);

        /*!
        * @brief Division operator by scalar
        * @param scalar Scalar value
        * @return Resultant Matrix object after division by scalar
        */
        Matrix operator/(double scalar);

        /*!
        * @brief Multiplication operator by scalar(Matrix*scalar)
        * @param scalar Scalar value
        * @return Resultant Matrix object after multiplication by scalar
        */
        Matrix operator*(double scalar);

        /*!
        * @brief Multiplication operator by scalar(scalar*Matrix)
        * @param scalar Scalar value
        * @return Resultant Matrix object after multiplication by scalar
        */
        friend Matrix operator*(double scalar, const Matrix& matrix);

        /*!
        * @brief Addition operator by scalar
        * @param scalar Scalar value
        * @return Resultant Matrix object after addition by scalar
        */
        Matrix operator+(double scalar);

        /*!
        * @brief Unary negation operator
        * @return Negated Matrix object
        */
        Matrix operator-() const;

        /*!
        * @brief Function to access elements of the matrix
        * @param i Row index (1-based)
        * @param j Column index (1-based)
        * @return Reference to the element at position (i, j)
        */
        double& operator()(int i, int j) const;

        /*!
        * @brief Function to compute the cross product of two 3x1 matrices
        * @param a First Matrix object
        * @param b Second Matrix object
        * @return Resultant Matrix object containing the cross product
        */
        static Matrix cross(const Matrix& a, const Matrix& b);

        /*!
        * @brief Function to compute the dot product of two 1x3 matrices
        * @param a First Matrix object
        * @param b Second Matrix object
        * @return Resultant dot product
        */
        static double dot(const Matrix &a, const Matrix &b);

        /*!
        * @brief Function to compute the transpose of the matrix
        * @return Transposed Matrix object
        */
        Matrix transpose() const;

        /*!
        * @brief Function to compute the inverse of the matrix
        * @return Inverted Matrix object
        */
        Matrix inverse() const;

        /*!
        * @brief Function to create an identity matrix of a specified size
        * @param size Size of the identity matrix (number of rows and columns)
        * @return Identity Matrix object
        */
        static Matrix createIdentityMatrix(int size);

        /*!
        * @brief Function to generate a range of values
        * @param start Starting value of the range
        * @param step Step size between consecutive values
        * @param end Ending value of the range
        * @return Matrix object containing the specified range of values
        */
        static Matrix range(int start, int step, int end);

        /**
        * @brief Extract a submatrix consisting of a single row from the current matrix.
        * @param row The index of the row to extract (1-based indexing).
        * @return A new matrix containing the specified row of the original matrix.
        * @note If the given row index is out of range, the function will print an error message to the standard error stream and exit the program.
        */
        Matrix subMatrix(int row);

        /*!
        * @brief Function to concatenate two matrices horizontally
        * @param matrix1 First Matrix object
        * @param matrix2 Second Matrix object
        * @return Resultant Matrix object after concatenation
        */
        static Matrix concatenateHorizontal(const Matrix& matrix1, const Matrix& matrix2);

        /*!
        * @brief Function to extract a submatrix consisting of a range of rows and columns
        * @param startRow Starting index of the row range to extract (1-based)
        * @param endRow Ending index of the row range to extract (1-based)
        * @param startCol Starting index of the column range to extract (1-based)
        * @param endCol Ending index of the column range to extract (1-based)
        * @return Submatrix containing the specified range of rows and columns
        */
        Matrix subMatrix(int startRow, int endRow, int startCol, int endCol);

        /*!
        * @brief Calculates the Euclidean norm of the matrix.
        * @return The Euclidean norm of the matrix.
        */
        double norm() const;

        /*!
        * @brief Function to print the matrix
        */
        void print();

        /**
        * @brief Get the number of rows in the matrix.
        *
        * @return The number of rows in the matrix.
        */
        int getRows() const;

        /**
        * @brief Get the number of columns in the matrix.
        * @return The number of columns in the matrix.
        */
        int getCols() const;



private:
    /*!
     * @brief Function to initialize the matrix with zeros
     */
    void initMatrix();

    int fil;        //!< Number of rows
    int col;        //!< Number of columns
    double **matrix;    //!< Pointer to the matrix data

};

#endif
