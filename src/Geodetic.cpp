//
// Created by miguel on 26/04/2024.
//

#include "../include/Geodetic.h"


/*
%--------------------------------------------------------------------------
%
% Geodetic.m
%
% Purpose:
%   geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
%   from given position vector (r [m])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix Geodetic(const Matrix& r){
    // Constantes globales
    const double R_equ = Constants::R_Earth;
    const double f = Constants::f_Earth;

    const double epsRequ = std::numeric_limits<double>::epsilon() * R_equ; // Criterio de convergencia
    const double e2 = f * (2.0 - f); // Cuadrado de la excentricidad

    double X = r(1, 1); // Coordenadas cartesianas
    double Y = r(1, 2);
    double Z = r(1, 3);
    double rho2 = X * X + Y * Y; // Cuadrado de la distancia desde el eje z

    // Comprobar validez de los datos de entrada
    if (r.norm() == 0.0) {
        throw std::invalid_argument("Invalid input in Geodetic constructor");
    }

    // Iteración
    double dZ = e2 * Z;
    double ZdZ, Nh,N;
    while (true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        double SinPhi = ZdZ / Nh; // Seno de la latitud geodésica
        N = R_equ / std::sqrt(1.0 - e2 * SinPhi * SinPhi);
        double dZ_new = N * e2 * SinPhi;
        if (std::abs(dZ - dZ_new) < epsRequ) {
            break;
        }
        dZ = dZ_new;
    }

    // Longitud, latitud, altitud
    double lon = atan2(Y, X);
    double lat = atan2(ZdZ, std::sqrt(rho2));
    double h = Nh - N;

    Matrix result(1, 3);
    result(1, 1) = lon;
    result(1, 2) = lat;
    result(1, 3) = h;

    return result;
}