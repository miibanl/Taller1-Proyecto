//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_ECCANOM_H
#define PROYECTO_ECCANOM_H

#include "SAT_Const.h"
#include "cmath"
#include <stdexcept>
#include <limits>


/*!
 * @file EccAnom.h
 * @brief Computes the eccentric anomaly for elliptic orbits
 * @param M Mean anomaly in [rad]
 * @param e Eccentricity of the orbit [0,1]
 * @return Eccentric anomaly in [rad]
 */
double EccAnom(double M, double e);


#endif //PROYECTO_ECCANOM_H
