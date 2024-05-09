//
// Created by miguel on 26/04/2024.
//

#include "../include/EccAnom.h"
/*
%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the eccentric anomaly for elliptic orbits
%
% Inputs:
%   M         Mean anomaly in [rad]
%   e         Eccentricity of the orbit [0,1]
%
% Output:
%             Eccentric anomaly in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
double EccAnom(double M, double e){

    int maxit = 15;
    int i = 1;

    // Starting value
    M = fmod(M, 2.0*Constants::pi);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = Constants::pi;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteration
    while (std::abs(f) > 1e2 * std::numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit) {
            throw std::runtime_error("Convergence problems in EccAnom");
        }
    }

    return E;
}


