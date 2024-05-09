//
// Created by miguel on 25/04/2024.
//

#include "../include/AccelPointMass.h"

/*
%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
*/

Matrix AccelPointMass(Matrix& r, Matrix& s, double GM) {
    // Vector de posición relativa del satélite con respecto a la masa puntual
    Matrix d = r - s;

    // Magnitud al cubo del vector de posición relativa y del vector de posición de la masa puntual
    double norm_d_cubed = pow(d.norm(), 3);
    double norm_s_cubed = pow(s.norm(), 3);

    // Aceleración
    Matrix a = (d / norm_d_cubed + s / norm_s_cubed) * -GM;

    return a;
}