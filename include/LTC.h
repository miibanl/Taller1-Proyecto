//
// Created by miguel on 14/05/2024.
//

#ifndef PROYECTO_LTC_H
#define PROYECTO_LTC_H

#include "Matrix.h"




/*!
 * @file LTC.h
 * @brief Transformation from Greenwich meridian system to local tangent coordinates
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @return Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 */
Matrix LTC(double lon, double lat);

#endif //PROYECTO_LTC_H
