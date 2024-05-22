//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_CHEB3D_H
#define PROYECTO_CHEB3D_H


#include "Matrix.h"
#include <stdexcept>

/*!
 * @file Cheb3D.h
 * @brief Chebyshev approximation of 3-dimensional vectors
 * @param t Value of the interval where the aproximation is evaluated
 * @param N Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval
 * @param Cx Coefficients of Chebyshev polyomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polyomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polyomial (z-coordinate)
 * @return Aproximated 3-Dimensional matrix at position t
 */
Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);

#endif //PROYECTO_CHEB3D_H
