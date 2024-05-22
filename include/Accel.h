//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_ACCEL_H
#define PROYECTO_ACCEL_H


#include "Matrix.h"
#include "../include/Global.h"
#include "../include/IERS.h"
#include "../include/TimeDiff.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/Mjday_TDB.h"
#include "../include/AccelPointMass.h"
#include "../include/AccelHarmonic.h"
#include "../include/JPL_Eph_DE430.h"




/*!
 * @file Accel.h
 * @brief Computes the acceleration of an Earth orbiting satellite due to
 *  - the Earth's harmonic gravity field,
 *  - the gravitational perturbations of the Sun and Moon
 *  - the solar radiation pressure and
 *  - the atmospheric drag
 * @param x Terrestrial Time (Modified Julian Date)
 * @param Y  Satellite state vector in the ICRF/EME2000 system
 * @return  Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
 */
Matrix Accel(double x, Matrix& Y);


#endif //PROYECTO_ACCEL_H
