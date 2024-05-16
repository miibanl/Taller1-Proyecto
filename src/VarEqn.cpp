//
// Created by miibanl on 16/05/2024.
//

#include <iostream>
#include "../include/VarEqn.h"
#include "../include/AccelHarmonic.h"
#include "../include/G_AccelHarmonic.h"


/*
%------------------------------------------------------------------------------
%
% VarEqn.m
%
% Purpose:
%   Computes the variational equations, i.e. the derivative of the state vector
%   and the state transition matrix
%
% Input:
%   x           Time since epoch in [s]
%   yPhi        (6+36)-dim vector comprising the state vector (y) and the
%               state transition matrix (Phi) in column wise storage order
%
% Output:
%   yPhip       Derivative of yPhi
%
% Last modified:   2015/08/12   M. Mahooti
%
%------------------------------------------------------------------------------
*/
Matrix VarEqn(double x,Matrix& yPhi) {

    Global::eop19620101(21413);


    //USAR SOLO AL HACER LOS TEST
    //Global::AuxParam::Mjd_UTC = 49746.1163541665;
    //Global::AuxParam::Mjd_TT = 49746.1170623147;
    //Global::AuxParam::n = 20;
    //Global::AuxParam::m = 20;


    Matrix iers = IERS(*Global::eopdata, Global::AuxParam::Mjd_UTC, 'l');
    Matrix timediff = timeDiff(iers(1, 3), iers(1, 9));
    double Mjd_UT1 = Global::AuxParam::Mjd_TT + (iers(1, 3) - timediff(1, 4)) / 86400.0;


    //% Transformation matrix
    Matrix P = PrecMatrix(Constants::MJD_J2000, Global::AuxParam::Mjd_TT + x / 86400.0);
    Matrix N = NutMatrix(Global::AuxParam::Mjd_TT + x / 86400.0);
    Matrix T = N * P;
    Matrix E = PoleMatrix(iers(1, 1), iers(1, 2)) * GHAMatrix(Mjd_UT1) * T;

    //% State vector components
    Matrix r = yPhi.subMatrix(1, 3, 1, 1);
    Matrix v = yPhi.subMatrix(4, 6, 1, 1);
    Matrix Phi(6, 6);



    //% State transition matrix
    for (int j = 1; j <= 6; j++) {
        Phi(1, j) = yPhi(6 * j + 1, 1);
        Phi(2, j) = yPhi(6 * j + 1 + 1, 1);
        Phi(3, j) = yPhi(6 * j + 1 + 2, 1);
        Phi(4, j) = yPhi(6 * j + 1 + 3, 1);
        Phi(5, j) = yPhi(6 * j + 1 + 4, 1);
        Phi(6, j) = yPhi(6 * j + 1 + 5, 1);
    }


    //% Acceleration and gradient
    Matrix a = AccelHarmonic(r, E, Global::AuxParam::n, Global::AuxParam::m);
    Matrix G = G_AccelHarmonic(r, E, Global::AuxParam::n, Global::AuxParam::m);

    //% Time derivative of state transition matrix
    Matrix yPhip(42, 1);
    Matrix dfdy(6, 6);

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i, j) = 0.0;                 //% dv/dr(i,j)
            dfdy(i + 3, j) = G(i, j);            //% da/dr(i,j)
            if (i == j) {
                dfdy(i, j + 3) = 1.0;
            } else {
                dfdy(i, j + 3) = 0.0;             //% dv/dv(i,j)
            }
            dfdy(i + 3, j + 3) = 0.0;             //% da/dv(i,j)
        }
    }

    Matrix Phip = dfdy * Phi;

    //% Derivative of combined state vector and state transition matrix
    for (int i = 1; i <= 3; i++) {
        yPhip(i, 1) = v(i, 1);                 //% dr/dt(i)
        yPhip(i + 3, 1) = a(i, 1);                 //% dv/dt(i)
    }

    for (int i = 1; i <= 6; i++) {
        for(int j=1;j<=6;j++){
            yPhip(6 * j + i,1) = Phip(i, j);     //% dPhi/dt(i,j)
        }
    }


    return yPhip;
}
