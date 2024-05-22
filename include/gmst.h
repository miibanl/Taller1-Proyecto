//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_GMST_H
#define PROYECTO_GMST_H

#include <cmath>
#include "SAT_Const.h"
#include "Frac.h"


/*!
 * @file gmst.h
 * @brief Greenwich Mean Sidereal Time
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return GMST in [rad]
 */
double gmst(double Mjd_UT1);


#endif //PROYECTO_GMST_H
