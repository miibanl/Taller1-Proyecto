//
// Created by miguel on 26/04/2024.
//

#include "../include/Frac.h"

/*
%--------------------------------------------------------------------------
%
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
double Frac(double x) {
    return x - floor(x);
}
