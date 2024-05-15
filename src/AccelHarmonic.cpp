//
// Created by miguel on 06/05/2024.
//

#include <iostream>
#include "../include/AccelHarmonic.h"

/*
%--------------------------------------------------------------------------
%
% AccelHarmonic.m
%
% Purpose:
%   Computes the acceleration due to the harmonic gravity field of the
%   central body
%
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max) {

    Global::GGM03S();

    double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    double gm = 398600.4415e9; // [m^3/s^2]; GGM03S
    Matrix r_bf = E * r;                 // Body-fixed position


    // Auxiliary quantities
    double d = r_bf.norm();                     // distance
    double latgc = asin(r_bf(3, 1) / d);
    double lon = atan2(r_bf(2, 1), r_bf(1, 1));

    Matrix pnm(n_max + 1, m_max + 1);
    Matrix dpnm(n_max + 1, m_max + 1);

    Legendre(n_max, m_max, latgc, pnm, dpnm);



    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0, q2 = q3, q1 = q2;
    for (int n = 0; n <= n_max; n++) {
        double b1 = (-gm / (d * d)) * pow((r_ref / d), n) * (n + 1);
        double b2 = (gm / d) * pow((r_ref / d), n);
        double b3 = (gm / d) * pow((r_ref / d), n);

        for (int m = 0; m <= m_max; ++m) {
                q1 = q1 + pnm(n + 1, m + 1) * ((*Global::Cnm)(n + 1, m + 1) * cos(m * lon) + (*Global::Snm)(n + 1, m + 1) * sin(m * lon));
                q2 = q2 + dpnm(n + 1, m + 1) * ((*Global::Cnm)(n + 1, m + 1) * cos(m * lon) + (*Global::Snm)(n + 1, m + 1) * sin(m * lon));
                q3 = q3 + m * pnm(n + 1, m + 1) * ((*Global::Snm)(n + 1, m + 1) * cos(m * lon) - (*Global::Cnm)(n + 1, m + 1) * sin(m * lon));
        }
        dUdr = dUdr + q1 * b1;
        dUdlatgc = dUdlatgc + q2 * b2;
        dUdlon = dUdlon + q3 * b3;
        q3 = 0;
        q2 = q3;
        q1 = q2;
    }



    //% Body-fixed acceleration
    double r2xy = pow(r_bf(1, 1), 2) + pow(r_bf(2, 1), 2);


    double ax = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) - (1.0 / r2xy * dUdlon) * r_bf(2, 1);
    double ay = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(2, 1) + (1.0 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1.0 / d * dUdr * r_bf(3, 1) + sqrt(r2xy) / (d * d) * dUdlatgc;



    Matrix a_bf(1,3);
    a_bf(1,1)=ax;
    a_bf(1,2)=ay;
    a_bf(1,3)=az;


    Matrix auxa_bf(a_bf.getCols(),a_bf.getCols());


    for (int i = 1; i <= a_bf.getCols(); ++i) {
        auxa_bf(i,1)=a_bf(1,i);
    }


    // Inertial acceleration
    return  E.transpose()*auxa_bf;

}
