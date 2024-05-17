//
// Created by miguel on 11/05/2024.
//

#include <iostream>
#include "../include/JPL_Eph_DE430.h"
#include "../include/Global.h"
#include "../include/Cheb3D.h"

/*
%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------
*/
void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth,
                     Matrix& r_Mars, Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus,
                     Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun){

    double JD = Mjd_TDB + 2400000.5;

    //MIRAR ESTO

    Global::DE430Coeff();
    int i;
    for (i = 1; i <= Global::PC->getCols(); i++) {
        if ((*Global::PC)(i, 1) <= JD && JD<=(*Global::PC)(i, 2) ) {
            break;
        }
    }
    Matrix PCtemp = (*Global::PC).subMatrix(i);


    double t1 = PCtemp(1,1)-2400000.5; //% MJD at start of interval

    double dt = Mjd_TDB - t1;


    Matrix temp = Matrix::range(231,13,270);


    Matrix Cx_Earth = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Earth = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Earth = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    temp = temp+39;
    Matrix Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);


    Cx_Earth = Matrix::concatenateHorizontal(Cx_Earth,Cx);
    Cy_Earth = Matrix::concatenateHorizontal(Cy_Earth,Cy);
    Cz_Earth = Matrix::concatenateHorizontal(Cz_Earth,Cz);


    double Mjd0;
    int j;
    if (0<=dt && dt<=16) {
        j = 0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32) {
        j = 1;
        Mjd0 = t1 + 16.0 * j;
    }



    r_Earth = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16.0, Cx_Earth.subMatrix(1,1,13*j+1,13*j+13),Cy_Earth.subMatrix(1,1,13*j+1,13*j+13),Cz_Earth.subMatrix(1,1,13*j+1,13*j+13)).transpose();


    temp = Matrix::range(441,13,480);
    Matrix Cx_Moon = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Moon = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Moon = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    for(int i=1;i<=7;i++) {
        temp = temp + 39;
        Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
        Cz = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
        Cx_Moon = Matrix::concatenateHorizontal(Cx_Moon,Cx);
        Cy_Moon = Matrix::concatenateHorizontal(Cy_Moon,Cy);
        Cz_Moon = Matrix::concatenateHorizontal(Cz_Moon,Cz);
    }
    if (0<=dt && dt<=4) {
        j = 0;
        Mjd0 = t1;
    }else if(4<dt && dt<=8) {
        j = 1;
        Mjd0 = t1 + 4 * j;
    }
    else if(8<dt && dt<=12) {
        j = 2;
        Mjd0 = t1 + 4 * j;
    }
    else if(12<dt && dt<=16) {
        j = 3;
        Mjd0 = t1 + 4 * j;
    }
    else if(16<dt && dt<=20) {
        j = 4;
        Mjd0 = t1 + 4 * j;
    }
    else if(20<dt && dt<=24) {
        j = 5;
        Mjd0 = t1 + 4 * j;
    }
    else if(24<dt && dt<=28) {
        j = 6;
        Mjd0 = t1 + 4 * j;
    }
    else if(28<dt && dt<=32) {
        j = 7;
        Mjd0 = t1 + 4 * j;
    }

    r_Moon = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon.subMatrix(1,1,13*j+1,13*j+13),Cy_Moon.subMatrix(1,1,13*j+1,13*j+13),Cz_Moon.subMatrix(1,1,13*j+1,13*j+13)).transpose();


    temp = Matrix::range(753,11,786);
    Matrix Cx_Sun = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Sun = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Sun = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    temp = temp+33;
    Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Cz = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    Cx_Sun = Matrix::concatenateHorizontal(Cx_Sun,Cx);
    Cy_Sun = Matrix::concatenateHorizontal(Cy_Sun,Cy);
    Cz_Sun = Matrix::concatenateHorizontal(Cz_Sun,Cz);
    if (0<=dt && dt<=16) {
        j = 0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    r_Sun = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, Cx_Sun.subMatrix(1,1,11*j+1,11*j+11),Cy_Sun.subMatrix(1,1,11*j+1,11*j+11),Cz_Sun.subMatrix(1,1,11*j+1,11*j+11)).transpose();


    temp = Matrix::range(3,14,45);
    Matrix Cx_Mercury = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Mercury = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Mercury = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    for (int i=1;i<=3;i++) {
        temp = temp + 42;
        Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
        Cz = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
        Cx_Mercury = Matrix::concatenateHorizontal(Cx_Mercury,Cx);
        Cy_Mercury = Matrix::concatenateHorizontal(Cy_Mercury,Cy);
        Cz_Mercury = Matrix::concatenateHorizontal(Cz_Mercury,Cz);
    }
    if (0<=dt && dt<=8) {
        j = 0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    }
    else if (16<dt && dt<=24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    }
    else if(24<dt && dt<=32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    r_Mercury = 1e3*Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, Cx_Mercury.subMatrix(1,1,14*j+1,14*j+14),Cy_Mercury.subMatrix(1,1,14*j+1,14*j+14),Cz_Mercury.subMatrix(1,1,14*j+1,14*j+14)).transpose();


    temp = Matrix::range(171,10,201);
    Matrix Cx_Venus = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Venus = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Venus = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    temp = temp+30;
    Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Cz = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    Cx_Venus = Matrix::concatenateHorizontal(Cx_Venus,Cx);
    Cy_Venus = Matrix::concatenateHorizontal(Cy_Venus,Cy);
    Cz_Venus = Matrix::concatenateHorizontal(Cz_Venus,Cz);
    if (0<=dt && dt<=16) {
        j = 0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32) {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }

    r_Venus = 1e3*Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, Cx_Venus.subMatrix(1,1,10*j+1,10*j+10),Cy_Venus.subMatrix(1,1,10*j+1,10*j+10),Cz_Venus.subMatrix(1,1,10*j+1,10*j+10)).transpose();


    temp = Matrix::range(309,11,342);
    Matrix Cx_Mars = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Mars = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Mars = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Mars = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, Cx_Mars.subMatrix(1,1,11*j+1,11*j+11),Cy_Mars.subMatrix(1,1,11*j+1,11*j+11),Cz_Mars.subMatrix(1,1,11*j+1,11*j+11)).transpose();


    temp = Matrix::range(342,8,366);
    Matrix Cx_Jupiter = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Jupiter = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Jupiter = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Jupiter = 1e3*Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, Cx_Jupiter.subMatrix(1,1,8*j+1,8*j+8),Cy_Jupiter.subMatrix(1,1,8*j+1,8*j+8),Cz_Jupiter.subMatrix(1,1,8*j+1,8*j+8)).transpose();


    temp = Matrix::range(366,7,387);
    Matrix Cx_Saturn = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Saturn = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Saturn = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Saturn = 1e3*Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, Cx_Saturn.subMatrix(1,1,7*j+1,7*j+7),Cy_Saturn.subMatrix(1,1,7*j+1,7*j+7),Cz_Saturn.subMatrix(1,1,7*j+1,7*j+7)).transpose();

    temp = Matrix::range(387,6,405);
    Matrix Cx_Uranus = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Uranus = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Uranus = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Uranus = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Uranus.subMatrix(1,1,6*j+1,6*j+6),Cy_Uranus.subMatrix(1,1,6*j+1,6*j+6),Cz_Uranus.subMatrix(1,1,6*j+1,6*j+6)).transpose();


    temp = Matrix::range(405,6,423);
    Matrix Cx_Neptune = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Neptune = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Neptune = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Neptune = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Neptune.subMatrix(1,1,6*j+1,6*j+6),Cy_Neptune.subMatrix(1,1,6*j+1,6*j+6),Cz_Neptune.subMatrix(1,1,6*j+1,6*j+6)).transpose();

    temp = Matrix::range(423,6,441);
    Matrix Cx_Pluto = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Pluto = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Pluto = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    j=0;
    Mjd0 = t1;
    r_Pluto = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Pluto.subMatrix(1,1,6*j+1,6*j+6),Cy_Pluto.subMatrix(1,1,6*j+1,6*j+6),Cz_Pluto.subMatrix(1,1,6*j+1,6*j+6)).transpose();


    temp = Matrix::range(819,10,839);
    Matrix Cx_Nutations = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Nutations = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    for (int i=1;i<=3;i++) {
        temp = temp + 20;
        Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
        Cx_Nutations = Matrix::concatenateHorizontal(Cx_Nutations,Cx);
        Cy_Nutations = Matrix::concatenateHorizontal(Cy_Nutations,Cy);
    }
    if (0<=dt && dt<=8) {
        j = 0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    }
    else if (16<dt && dt<=24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    }
    else if(24<dt && dt<=32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }
    Matrix zeros(10,1);
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations.subMatrix(1,1,10*j+1,10*j+10),Cy_Nutations.subMatrix(1,1,10*j+1,10*j+10),zeros).transpose();
    Nutations(3,1)=0.0;


    temp = Matrix::range(899,10,929);
    Matrix Cx_Librations = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
    Matrix Cy_Librations = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
    Matrix Cz_Librations = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
    for (int i=1;i<=3;i++) {
        temp = temp + 30;
        Cx = PCtemp.subMatrix(1,1,temp(1,1),temp(1,2)-1);
        Cy = PCtemp.subMatrix(1,1,temp(1,2),temp(1,3)-1);
        Cz = PCtemp.subMatrix(1,1,temp(1,3),temp(1,4)-1);
        Cx_Librations = Matrix::concatenateHorizontal(Cx_Librations,Cx);
        Cy_Librations = Matrix::concatenateHorizontal(Cy_Librations,Cy);
        Cz_Librations = Matrix::concatenateHorizontal(Cz_Librations,Cz);
    }
    if (0<=dt && dt<=8) {
        j = 0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16) {
        j = 1;
        Mjd0 = t1 + 8 * j;
    }
    else if (16<dt && dt<=24) {
        j = 2;
        Mjd0 = t1 + 8 * j;
    }
    else if(24<dt && dt<=32) {
        j = 3;
        Mjd0 = t1 + 8 * j;
    }

    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations.subMatrix(1,1,10*j+1,10*j+10),Cy_Librations.subMatrix(1,1,10*j+1,10*j+10),Cz_Librations.subMatrix(1,1,10*j+1,10*j+10)).transpose();


    double EMRAT = 81.30056907419062; //% DE430
    double EMRAT1 = 1.0/(1.0+EMRAT);
    r_Earth = r_Earth-EMRAT1*r_Moon;
    r_Mercury = -r_Earth+r_Mercury;
    r_Venus = -r_Earth+r_Venus;
    r_Mars = -r_Earth+r_Mars;
    r_Jupiter = -r_Earth+r_Jupiter;
    r_Saturn = -r_Earth+r_Saturn;
    r_Uranus = -r_Earth+r_Uranus;
    r_Neptune = -r_Earth+r_Neptune;
    r_Pluto = -r_Earth+r_Pluto;
    r_Sun = -r_Earth+r_Sun;

}