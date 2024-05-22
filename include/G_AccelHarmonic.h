//
// Created by miibanl on 16/05/2024.
//

#ifndef PROYECTO_G_ACCELHARMONIC_H
#define PROYECTO_G_ACCELHARMONIC_H


#include "Matrix.h"
#include "../include/AccelHarmonic.h"



/*!
 * @file G_AccelHarmonic.h
 * @brief Computes the gradient of the Earth's harmonic gravity field
 * @param r Satellite position vector in the true-of-date system
 * @param U Transformation matrix to body-fixed system
 * @param n_max Gravity model degree
 * @param m_max Gravity model order
 * @return Gradient (G=da/dr) in the true-of-date system
 */
Matrix G_AccelHarmonic(Matrix& r,Matrix& U,int n_max,int m_max);


#endif //PROYECTO_G_ACCELHARMONIC_H
