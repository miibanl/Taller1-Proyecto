//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_POSITION_H
#define PROYECTO_POSITION_H

#include "cmath"
#include "Matrix.h"
#include "SAT_Const.h"


/*!
 * @file Position.h
 * @brief Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 * @param lon (Longitude [rad])
 * @param lat (Latitude [rad]
 * @param h (Altitude [m]
 * @return Position vector (r [m])
 */
Matrix Position(double lon, double lat, double h);


#endif //PROYECTO_POSITION_H
