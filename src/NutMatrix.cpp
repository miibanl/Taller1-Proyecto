//
// Created by miguel on 11/05/2024.
//

#include "../include/NutMatrix.h"


/*
%--------------------------------------------------------------------------
%
% NutMatrix.m
%
% Purpose:
%   Transformation from mean to true equator and equinox
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   NutMat    Nutation matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix NutMatrix (double Mjd_TT) {

    //% Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    //% Nutation in longitude andobliquity
    Matrix nutAngles = NutAngles(Mjd_TT);

    //% Transformation from mean to true equator andequinox
    return (R_x(-eps - nutAngles(1,2)) * R_z(-nutAngles(1,1)) * R_x(+eps));

}