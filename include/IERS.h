//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H


#include "Matrix.h"
#include <stdexcept>
#include "SAT_Const.h"

/*!
 * @file IERS.h
 * @brief Provides management of IERS time and polar motion data.
 * @param eop Matrix containing the Earth orientation parameters.
 * @param Mjd_UTC Modified Julian Date for UTC.
 * @param interp Character indicating the interpolation method ('l' for linear, 'n' for no interpolation).
 * @return Matrix containing the interpolated or non-interpolated Earth orientation parameters.
 */
Matrix IERS(const Matrix& eop, double Mjd_UTC, char interp= 'n');


#endif //PROYECTO_IERS_H
