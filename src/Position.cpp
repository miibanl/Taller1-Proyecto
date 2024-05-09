//
// Created by miguel on 26/04/2024.
//

#include "../include/Position.h"

/*
%--------------------------------------------------------------------------
%
% Position.m
%
% Purpose:
%   Position vector (r [m]) from geodetic coordinates (Longitude [rad],
%   latitude [rad], altitude [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix Position(double lon, double lat, double h) {
    const double R_equ = Constants::R_Earth;
    const double f = Constants::f_Earth;

    double e2 = f * (2.0 - f);  // Square of eccentricity
    double CosLat = cos(lat);   // (Co)sine of geodetic latitude
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    Matrix r(1, 3);
    r(1, 1) = (N + h) * CosLat * cos(lon);
    r(1, 2) = (N + h) * CosLat * sin(lon);
    r(1, 3) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}