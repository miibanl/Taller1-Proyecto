//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_GHAMATRIX_H
#define PROYECTO_GHAMATRIX_H

#include "Matrix.h"
#include "R_z.h"
#include "gast.h"


/*!
 * @file GHAMatrix.h
 * @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return Greenwich Hour Angle matrix
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif //PROYECTO_GHAMATRIX_H
