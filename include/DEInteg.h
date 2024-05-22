//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_DEINTEG_H
#define PROYECTO_DEINTEG_H

#include <limits>
#include "Matrix.h"
#include "../include/sign.h"

/*!
 * @file DEInteg.h
 * @brief Numerical integration methods for ordinaray differential equations. This module provides implemenation of the variable order variable stepsize multistep method of Shampine & Gordon.
 * @param func function that computes the derivatives of the system state vector
 * @param t initial time
 * @param tout time at the solution is desired
 * @param relerr relative error tolerance
 * @param abserr absolute error tolerance
 * @param n_eqn number of ecuations
 * @param y state vector at the initial time t. Updated to the state vector at tout
 */
void DEInteg(Matrix (*func)(double,Matrix&), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y);


#endif //PROYECTO_DEINTEG_H
