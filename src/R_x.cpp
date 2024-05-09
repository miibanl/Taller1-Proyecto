//
// Created by miibanl on 11/04/2024.
//

#include "../include/R_x.h"

/*
%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------
function [rotmat] = R_x(angle)
*/
Matrix R_x(double alpha){
    double C=cos(alpha);
    double S=sin(alpha);

    Matrix rotmat(3,3);

    rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;

    return rotmat;
}


