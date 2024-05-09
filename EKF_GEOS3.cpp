//
// Created by miguel on 06/05/2024.
//

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include "./include/Matrix.h"
#include "./include/R_x.h"
#include "./include/R_y.h"
#include "./include/R_z.h"
#include "./include/SAT_Const.h"
#include "./include/Global.h"
#include "./include/sign.h"
#include "./include/TimeDiff.h"
#include "./include/Unit.h"
#include "./include/AccelPointMass.h"
#include "./include/AzElPa.h"
#include "./include/Cheb3D.h"
#include "./include/EccAnom.h"
#include "./include/Frac.h"
#include "./include/Geodetic.h"
#include "./include/Legendre.h"
#include "./include/MeanObliquity.h"
#include "./include/Mjday.h"
#include "./include/Mjday_TDB.h"
#include "./include/NutAngles.h"
#include "./include/Position.h"
#include "./include/IERS.h"

/*
%--------------------------------------------------------------------------
%
%  Initial Orbit Determination using Gauss and Extended Kalman Filter methods
%
% References:
%   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
%   Applications", Springer Verlag, Heidelberg, 2000.
%
%   D. Vallado, "Fundamentals of Astrodynamics and Applications",
%   4th Edition, 2013.
%
%   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
%
% Last modified:   2020/03/16   Meysam Mahooti
%--------------------------------------------------------------------------
*/

/*
int main(){

    Global::DE430Coeff();
    Global::GGM03S();
    Global::eop19620101(21413);


    int nobs=46;
    Matrix obs(nobs,4);

    for(int i=1;i<=nobs;i++){



        double Y = str2num(tline(1:4));
        double M = str2num(tline(6:7));
        double D = str2num(tline(9:10));
        double h = str2num(tline(13:14));
        double m = str2num(tline(16:17));
        double s = str2num(tline(19:24));
        double az = str2num(tline(26:33));
        double el = str2num(tline(36:42));
        double Dist = str2num(tline(45:54));
        obs(i,1) = Mjday(Y,M,D,h,m,s);
        obs(i,2) = Constants::Rad*az;
        obs(i,3) = Constants::Rad*el;
        obs(i,4) = 1e3*Dist;
    }





    return 0;
}*/