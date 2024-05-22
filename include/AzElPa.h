//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_AZELPA_H
#define PROYECTO_AZELPA_H

#include "Matrix.h"
#include "SAT_Const.h"

/*!
 * @file AzElPa.h
 * @brief Computes azimuth, elevation and partials from local tangent coordinates
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @param Az   Azimuth [rad]
 * @param El   Elevation [rad]
 * @param dAds Partials of azimuth w.r.t. s
 * @param dEds Partials of elevation w.r.t. s
 */
void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds);


#endif //PROYECTO_AZELPA_H
