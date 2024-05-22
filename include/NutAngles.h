//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_NUTANGLES_H
#define PROYECTO_NUTANGLES_H

#include "Matrix.h"
#include "SAT_Const.h"


/*!
 * @file NutAngles.h
 * @brief Nutation in longitude and obliquity
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Nutation Angles
 */
Matrix NutAngles(double Mjd_TT);

#endif //PROYECTO_NUTANGLES_H
