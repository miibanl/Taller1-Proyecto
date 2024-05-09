//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_CHEB3D_H
#define PROYECTO_CHEB3D_H


#include "Matrix.h"
#include <stdexcept>


Matrix Cheb3D(double t, int n, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);

#endif //PROYECTO_CHEB3D_H
