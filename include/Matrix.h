#ifndef _MATRIX_
#define _MATRIX_

#include <cmath>

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, double v[], int n);
        Matrix(const Matrix& m);
        ~Matrix();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        Matrix operator/(double scalar);
        Matrix operator*(double scalar);
        friend Matrix operator*(double scalar, const Matrix& matrix);
        Matrix operator+(double scalar);
        Matrix operator-() const;




        double& operator()(int i, int j) const;

        static Matrix cross(const Matrix& a, const Matrix& b);
        static double dot(const Matrix &a, const Matrix &b);
        Matrix transpose() const;
        Matrix inverse() const;
        static Matrix createIdentityMatrix(int size);


        static Matrix range(int start, int step, int end);
        Matrix subMatrix(int row);
    static Matrix concatenateHorizontal(const Matrix& matrix1, const Matrix& matrix2);
    Matrix subMatrix(int startRow, int endRow, int startCol, int endCol);
    void copyColumnFrom(const Matrix& source, int sourceCol, int destCol);











        void print();

        double norm() const;
        int getRows() const;
        int getCols() const;



private:
        void initMatrix();

        int fil;
        int col;
        double **matrix;

};

#endif
