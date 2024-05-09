//
// Created by miguel on 26/04/2024.
//

#include "../include/AzElPa.h"

/*
%--------------------------------------------------------------------------
%
% Purpose:
%  Computes azimuth, elevation and partials from local tangent coordinates
%
% Input:
%   s      Topocentric local tangent coordinates (East-North-Zenith frame)
%
% Outputs:
%   A      Azimuth [rad]
%   E      Elevation [rad]
%   dAds   Partials of azimuth w.r.t. s
%   dEds   Partials of elevation w.r.t. s
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {

    double rho = sqrt(s(1, 1) * s(1, 1) + s(1, 2) * s(1, 2));

    // Angles
    Az = atan2(s(1, 1), s(1, 2));

    if (Az < 0.0) {
        Az = Az + Constants::pi2;
    }

    El = atan(s(1, 3) / rho);

    // Partials
    double s_prodEscalar = s(1, 1) * s(1, 1) + s(1, 2) * s(1, 2) + s(1, 3) * s(1, 3);
    dAds(1, 1) = s(1, 2) / (rho * rho);
    dAds(1, 2) = -s(1, 1) / (rho * rho);
    dAds(1, 3) = 0.0;

    dEds(1, 1) = (-s(1, 1) * s(1, 3) / rho) / s_prodEscalar;
    dEds(1, 2) = (-s(1, 2) * s(1, 3) / rho) /s_prodEscalar;
    dEds(1, 3) = rho / s_prodEscalar;
}
