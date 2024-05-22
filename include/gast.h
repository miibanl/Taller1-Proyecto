//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_GAST_H
#define PROYECTO_GAST_H

#include "SAT_Const.h"
#include "gmst.h"
#include "EqnEquinox.h"

/*!
 * @file gast.h
 * @brief Greenwich Apparent Sidereal Time
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return GAST in [rad]
 */
double gast(double Mjd_UT1);

#endif //PROYECTO_GAST_H
