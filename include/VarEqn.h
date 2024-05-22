//
// Created by miibanl on 16/05/2024.
//

#ifndef PROYECTO_VAREQN_H
#define PROYECTO_VAREQN_H


#include "Matrix.h"
#include "Global.h"
#include "IERS.h"
#include "TimeDiff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"

/*!
 * @file VarEqn.h
 * @brief Computes the variational equations, i.e. the derivative of the state vector and the state transition matrix
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
 * @return Derivative of yPhi
 */
Matrix VarEqn(double x,Matrix& yPhi);

#endif //PROYECTO_VAREQN_H
