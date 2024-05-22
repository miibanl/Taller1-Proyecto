//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_PRECMATRIX_H
#define PROYECTO_PRECMATRIX_H

#include "Matrix.h"
#include "SAT_Const.h"
#include "R_z.h"
#include "R_y.h"


/*!
 * @file PrecMatrix.h
 * @brief Precession transformation of equatorial coordinates
 * @param Mjd_1 Epoch given (Modified Julian Date TT)
 * @param Mjd_2 Epoch to precess to (Modified Julian Date TT)
 * @return Precession transformation matrix
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);


#endif //PROYECTO_PRECMATRIX_H
