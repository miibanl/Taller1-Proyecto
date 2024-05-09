//
// Created by miibanl on 18/04/2024.
//


#include "../include/R_y.h"

/*
%--------------------------------------------------------------------------
%  input:
%    angle       - angle of rotation [rad]
%
%  output:
%    rotmat      - vector result
%--------------------------------------------------------------------------
*/

Matrix R_y(double alpha){
    double C=cos(alpha);
    double S=sin(alpha);

    Matrix rotmat(3,3);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;

    return rotmat;
}

