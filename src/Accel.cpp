//
// Created by miguel on 11/05/2024.
//

#include "../include/Accel.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/Mjday_TDB.h"
#include "../include/AccelPointMass.h"
#include "../include/AccelHarmonic.h"
#include "../include/JPL_Eph_DE430.h"

/*
%--------------------------------------------------------------------------
%
% Accel.m
%
% Purpose:
%   Computes the acceleration of an Earth orbiting satellite due to
%    - the Earth's harmonic gravity field,
%    - the gravitational perturbations of the Sun and Moon
%    - the solar radiation pressure and
%    - the atmospheric drag
%
% Inputs:
%   Mjd_TT      Terrestrial Time (Modified Julian Date)
%   Y           Satellite state vector in the ICRF/EME2000 system
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix Accel(double x,Matrix& Y){




    Matrix iers = IERS(*Global::eopdata,(Global::AuxParam::Mjd_UTC + x/86400.0),'l');
    Matrix timediff = timeDiff(iers(1,3),iers(1,9));
    double Mjd_UT1 = Global::AuxParam::Mjd_UTC + x/86400.0 + iers(1,3)/86400.0;
    double Mjd_TT = Global::AuxParam::Mjd_UTC + x/86400.0 + timediff(1,4)/86400.0;

    Matrix P = PrecMatrix(Constants::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(iers(1,1),iers(1,2)) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury(3, 1);
    Matrix r_Venus(3, 1);
    Matrix r_Earth(3, 1);
    Matrix r_Mars(3, 1);
    Matrix r_Jupiter(3, 1);
    Matrix r_Saturn(3, 1);
    Matrix r_Uranus(3, 1);
    Matrix r_Neptune(3, 1);
    Matrix r_Pluto(3, 1);
    Matrix r_Moon(3, 1);
    Matrix r_Sun(3, 1);
    JPL_Eph_DE430(MJD_TDB,r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun);

    Matrix auxY(3,1);
    for(int i=1;i<=3;i++){
        auxY(i,1)=Y(i,1);
    }

    //% Acceleration due to harmonic gravity field
    Matrix a = AccelHarmonic(auxY, E, Global::AuxParam::n, Global::AuxParam::m);

    //% Luni-solar perturbations
    if (Global::AuxParam::sun) {
        a = a + AccelPointMass(auxY,r_Sun,Constants::GM_Sun);
    }

    if (Global::AuxParam::moon) {
        a = a + AccelPointMass(auxY,r_Moon,Constants::GM_Moon);
    }

    //% Planetary perturbations
    if (Global::AuxParam::planets) {
        a = a + AccelPointMass(auxY,r_Mercury,Constants::GM_Mercury);
        a = a + AccelPointMass(auxY,r_Venus,Constants::GM_Venus);
        a = a + AccelPointMass(auxY,r_Mars,Constants::GM_Mars);
        a = a + AccelPointMass(auxY,r_Jupiter,Constants::GM_Jupiter);
        a = a + AccelPointMass(auxY,r_Saturn,Constants::GM_Saturn);
        a = a + AccelPointMass(auxY,r_Uranus,Constants::GM_Uranus);
        a = a + AccelPointMass(auxY,r_Neptune,Constants::GM_Neptune);
        a = a + AccelPointMass(auxY,r_Pluto,Constants::GM_Pluto);
    }


    Matrix dY(6,1);

    for(int i=4;i<=6;i++){
        dY(i-3,1)=Y(i,1);
    }
    for(int i=7;i<=9;i++){
        dY(i-3,1)=a(i-6,1);
    }

    return dY;

}
