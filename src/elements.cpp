//
// Created by miguel on 14/05/2024.
//

#include "../include/elements.h"

/*
% Input:
%    y        State vector (x,y,z,vx,vy,vz)
%
% Outputs:
%    p        semilatus rectum [m]
%    a        Semimajor axis
%    e        Eccentricity
%    i        Inclination [rad]
%    Omega    Longitude of the ascending node [rad]
%    omega    Argument of pericenter [rad]
%    M        Mean anomaly [rad]
%
% Notes:
%   The function cannot be used with state vectors describing a circular
%   or non-inclined orbit.
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix elements(const Matrix& y){

    Matrix r(1,3);
    Matrix result(1,7);


    for(int i=1;i<=3;i++) {
        r(1,i) = y(1,i);                //% Position
    }

    Matrix v(1, 3);

    for(int i=4;i<=6;i++) {
        v(1, i-3) = y(1, i);          //% Velocity
    }

    Matrix h = Matrix::cross(r,v);        //% Areal velocity


    double magh = h.norm();
    double p = (magh*magh)/Constants::GM_Earth;
    double H = h.norm();

    double Omega = atan2(h(1,1), -h(1,2) );       //% Long. ascend. node
    Omega = fmod(Omega,Constants::pi2);
    double i = atan2(sqrt((h(1,1)*h(1,1))+(h(1,2)*h(1,2))), h(1,3) ); //% Inclination
    double u = atan2 ( r(1,3)*H, -r(1,1)*h(1,2)+r(1,2)*h(1,1) );    //% Arg. of latitude
    double R  = r.norm();     //% Distance


    double a = 1.0/(2.0/R-Matrix::dot(v,v)/Constants::GM_Earth);    // % Semi-major axis
    double eCosE = 1.0-R/a;                                    // % e*cos(E)
    double eSinE = Matrix::dot(r,v)/sqrt(Constants::GM_Earth*a);           //% e*sin(E)

    double e2 = eCosE*eCosE + eSinE*eSinE;
    double e  = sqrt(e2);                                     //% Eccentricity
    double E  = atan2(eSinE,eCosE);                           //% Eccentric anomaly

    double M  = fmod(E-eSinE,Constants::pi2);             //% Mean anomaly

    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          //% True anomaly

    double omega = fmod(u-nu,Constants::pi2);                             //% Arg. of perihelion


    result(1,1)=p;
    result(1,2)=a;
    result(1,3)=e;
    result(1,4)=i;
    result(1,5)=Omega;
    result(1,6)=omega;
    result(1,7)=M;

    return result;
}
