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


    Matrix Cx_EarthAux = Matrix::concatenateHorizontal(Cx_Earth,Cx);
    Matrix Cy_EarthAux = Matrix::concatenateHorizontal(Cy_Earth,Cy);
    Matrix Cz_EarthAux = Matrix::concatenateHorizontal(Cz_Earth,Cz);



    double Mjd0;
    int j;
    if (0<=dt && dt<=16) {
        j = 0;
        Mjd0 = t1;
    }else if(16<dt && dt<=32) {
        j = 1;
        Mjd0 = t1 + 16.0 * j;
    }


    Cx_EarthAux.subMatrix(1,1,13*j+1,13*j+13).print();

    r_Earth = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16.0, Cx_EarthAux.subMatrix(1,1,13*j+1,13*j+13),Cy_EarthAux.subMatrix(1,1,13*j+1,13*j+13),Cz_EarthAux.subMatrix(1,1,13*j+1,13*j+13)).transpose();

    r_Earth.print();

}