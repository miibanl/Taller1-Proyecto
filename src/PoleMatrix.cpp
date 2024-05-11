//
// Created by miguel on 11/05/2024.
//

#include "../include/PoleMatrix.h"


/*
%--------------------------------------------------------------------------
%
% PoleMatrix.m
%
% Purpose:
%   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
%   for a given date
%
% Input:
%   Pole coordinte(xp,yp)
%
% Output:
%   PoleMat   Pole matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix PoleMatrix(double xp, double yp){
    return (R_y(-xp) * R_x(-yp));
}
