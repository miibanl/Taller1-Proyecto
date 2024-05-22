//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_GEODETIC_H
#define PROYECTO_GEODETIC_H


#include "Matrix.h"
#include <cmath>
#include <stdexcept>
#include "SAT_Const.h"
#include <limits>


/*!
 * @file Geodetic.h
 * @brief geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
 * @param r Position vector[m]
 * @return Matrix(lon,lat,alt)
 */
Matrix Geodetic(const Matrix& r);


#endif //PROYECTO_GEODETIC_H
