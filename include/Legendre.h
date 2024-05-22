//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_LEGENDRE_H
#define PROYECTO_LEGENDRE_H


#include "Matrix.h"

/*!
 * @file Legendre.h
 * @brief Computes the associated Legendre functions and their derivatives
 * @param n The maximum degree
 * @param m The maximum order
 * @param fi The argument of the Legendre functions [rad]
 * @param pnm The matrix that will hold the values of the associated Legendre functions
 * @param dpnm The matrix that will hold the values of the derivatives of the associated Legendre functions
 */
void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm);



#endif //PROYECTO_LEGENDRE_H
