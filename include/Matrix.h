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


        double& operator()(int i, int j) const;

        static Matrix cross(const Matrix& a, const Matrix& b);
        static double dot(const Matrix &a, const Matrix &b);



    void print();

        double norm() const;
        int getRows() const;
 
    private:
        void initMatrix();

        int fil;
        int col;
        double **matrix;

};

#endif
