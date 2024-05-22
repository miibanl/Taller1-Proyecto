//
// Created by miguel on 06/05/2024.
//

#ifndef PROYECTO_ACCELHARMONIC_H
#define PROYECTO_ACCELHARMONIC_H

#include "Matrix.h"
#include "../include/Legendre.h"
#include "../include/Global.h"


/*!
 * @file AccelHarmonic.h
 * @brief Computes the acceleration due to the harmonic gravity field of the central body
 * @param r Satellite position vector in the inertial system
 * @param E Transformation matrix to body-fixed system
 * @param n_max Maximum degree
 * @param m_max Maximum order (m_max<=n_max; m_max=0 for zonals, only)
 * @return Acceleration (a=d^2r/dt^2)
 */
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);


#endif //PROYECTO_ACCELHARMONIC_H
