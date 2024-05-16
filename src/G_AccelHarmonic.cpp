//
// Created by miibanl on 16/05/2024.
//

#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"

/*
%--------------------------------------------------------------------------
%
% G_AccelHarmonic.m
%
% Purpose:
%   Computes the gradient of the Earth's harmonic gravity field
%
% Inputs:
%   r           Satellite position vector in the true-of-date system
%   U           Transformation matrix to body-fixed system
%   n           Gravity model degree
%   m 			Gravity model order
%
% Output:
%   G    		Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix G_AccelHarmonic(const Matrix& r,const Matrix& U,int n_max,int m_max) {

    double d = 1.0;   //% Position increment [m]

    Matrix G(3, 3);

    //% Gradient
    for (int i = 1; i <= 3; i++) {
        //% Set offset in i-th component of the position vector
        Matrix dr(3,1);
        dr(i,1) = d;
        //% Acceleration difference
        Matrix aux1=r + dr / 2;
        Matrix aux2=r - dr / 2;
        Matrix da = AccelHarmonic(aux1, U, n_max, m_max) - AccelHarmonic(aux2, U, n_max, m_max);
        //% Derivative with respect to i-th axis
        G(1,i) = (da / d)(1,i);
        G(2,i) = (da / d)(2,i);
        G(3,i) = (da / d)(3,i);
    }

    return G;
}


