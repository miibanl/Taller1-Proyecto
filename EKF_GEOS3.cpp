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

    Global::GEOS3();

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


    double sigma_range = 92.5;          //% [m]
    double sigma_az = 0.0224*Constants::Rad; //% [rad]
    double sigma_el = 0.0139*Constants::Rad; //% [rad]

    //% Kaena Point station
    double lat = Constants::Rad*21.5748;     //% [rad]
    double lon = Constants::Rad*(-158.2706); //% [rad]
    double alt = 300.20;                //% [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

//   %[r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
//              %                 Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
//   % [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
//               %                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);


    //Y0_apr = [r2;v2];

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    Global::AuxParam::Mjd_UTC = Mjd_UTC;
    Global::AuxParam::n      = 20;
    Global::AuxParam::m      = 20;
    Global::AuxParam::sun     = 1;
    Global::AuxParam::moon    = 1;
    Global::AuxParam::planets = 1;

    int n_eqn  = 6;


    double Y = DEInteg(&Accel, 0, -(obs(9,1)-Mjd0)*86400.0, 1e-13, 1e-6, 6, Y0_apr);






    return 0;
}
 */