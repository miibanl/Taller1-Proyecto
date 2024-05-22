//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_NUTMATRIX_H
#define PROYECTO_NUTMATRIX_H


#include "Matrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"


/*!
 * @file NutMatrix.h
 * @brief Transformation from mean to true equator and equinox
 * @param Mjd_TT  Modified Julian Date (Terrestrial Time)
 * @return Nutation matrix
 */
Matrix NutMatrix (double Mjd_TT);

#endif //PROYECTO_NUTMATRIX_H
