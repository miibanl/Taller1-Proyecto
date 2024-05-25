//
// Created by miguel on 06/05/2024.
//


/*!
 * @file EKF_GEOS3.cpp
 * @brief  Test file for Initial Orbit Determination using Gauss and Extended Kalman Filter methods.
 *
 * This test file is designed to validate the main functionalities of the orbit determination code.
 * It uses a set of predefined observations and checks the correctness of the output.
 *
 * Last modified: 2020/03/16 by Meysam Mahooti
 *
 * @author miguel
 */
#include <iostream>
#include <fstream>
#include "./include/Matrix.h"
#include "./include/R_z.h"
#include "./include/SAT_Const.h"
#include "./include/Global.h"
#include "./include/TimeDiff.h"
#include "./include/AzElPa.h"
#include "./include/Cheb3D.h"
#include "./include/Mjday.h"
#include "./include/Position.h"
#include "./include/IERS.h"
#include "./include/Accel.h"
#include "./include/DEInteg.h"
#include "./include/LTC.h"
#include "./include/VarEqn.h"
#include "./include/TimeUpdate.h"
#include "./include/MeasUpdate.h"

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

int main() {

    Global::DE430Coeff();
    Global::GGM03S();
    Global::eop19620101(21413);


    int nobs = 46;
    Matrix obs(nobs, 4);




    std::ifstream infile("./data/GEOS3.txt");
    std::string line;

    for (int i = 1; i <= nobs && std::getline(infile, line); ++i) {
        if (line.empty()) { break; }

        int Y = std::stoi(line.substr(0, 4));
        int M = std::stoi(line.substr(5, 2));
        int D = std::stoi(line.substr(8, 2));
        int h = std::stoi(line.substr(12, 2));
        int m = std::stoi(line.substr(15, 2));
        int s = std::stoi(line.substr(18, 6));
        double az = std::stod(line.substr(25, 8));
        double el = std::stod(line.substr(35, 7));
        double Dist = std::stod(line.substr(44, 10));

        obs(i, 1) = Mjday(Y, M, D, h, m, s);
        obs(i, 2) = (Constants::Rad) * az;
        obs(i, 3) = (Constants::Rad) * el;
        obs(i, 4) = 1e3 * Dist;
    }

    infile.close();


    double sigma_range = 92.5;          //% [m]
    double sigma_az = 0.0224 * Constants::Rad; //% [rad]
    double sigma_el = 0.0139 * Constants::Rad; //% [rad]

    //% Kaena Point station
    double lat = Constants::Rad * 21.5748;     //% [rad]
    double lon = Constants::Rad * (-158.2706); //% [rad]
    double alt = 300.20;                //% [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);

//   %[r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
//              %                 Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
//   % [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
//               %                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);


    Matrix Y0_apr(6, 1);

    Y0_apr(1, 1) = 6221397.62857869;
    Y0_apr(2, 1) = 2867713.77965738;
    Y0_apr(3, 1) = 3006155.98509949;
    Y0_apr(4, 1) = 4645.04725161807;
    Y0_apr(5, 1) = -2752.21591588205;
    Y0_apr(6, 1) = -7507.99940987033;


    double Mjd0 = Mjday(1995, 1, 29, 02, 38, 0);

    double Mjd_UTC = obs(9, 1);

    Global::AuxParam::Mjd_UTC = Mjd_UTC;
    Global::AuxParam::n = 20;
    Global::AuxParam::m = 20;
    Global::AuxParam::sun = 1;
    Global::AuxParam::moon = 1;
    Global::AuxParam::planets = 1;

    int n_eqn = 6;

    DEInteg(Accel, 0, -(obs(9, 1) - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);

    Matrix Y(6, 1);
    for (int i = 1; i <= 6; i++) {
        Y(i, 1) = Y0_apr(i, 1);
    }


    Matrix P(6, 6);

    for (int i = 1; i <= 3; i++) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; i++) {
        P(i, i) = 1e3;
    }

    Matrix LT = LTC(lon, lat);

    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    //% Measurement loop
    double t = 0.0;

    Matrix iers(1,9);
    Matrix timediff(1,5);

    Matrix Y_old(6, 1);
    for(int i = 1;i<=nobs;i++) {
        //% Previous step
        double t_old = t;


        for (int i = 1; i <= 6; i++) {
            Y_old(i, 1) = Y(i, 1);
        }

        //% Time increment andpropagation
        Mjd_UTC = obs(i, 1);                       //% Modified Julian Date
        t = (Mjd_UTC - Mjd0) * 86400.0;         //% Time since epoch[s]

        iers = IERS(*Global::eopdata, Mjd_UTC, 'l');
        timediff = timeDiff(iers(1,3),iers(1,9));
        double Mjd_TT = Mjd_UTC + timediff(1,4) / 86400.0;
        double Mjd_UT1 = Mjd_TT + (iers(1,3) - timediff(1,4)) / 86400.0;
        Global::AuxParam::Mjd_UTC = Mjd_UTC;
        Global::AuxParam::Mjd_TT = Mjd_TT;

        for (int ii = 1;ii<=6;ii++) {
            yPhi(ii,1) = Y_old(ii,1);
            for (int j = 1;j<=6;j++) {
                if (ii == j) {
                    yPhi(6 * j + ii,1) = 1.0;
                }else {
                    yPhi(6 * j + ii,1) = 0.0;
                }
            }
        }


        DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);



        //% Extract state transition matrices
        for(int j = 1;j<=6;j++) {
            Phi(1,j) = yPhi(6 * j + 1,1);
            Phi(2,j) = yPhi(6 * j + 2,1);
            Phi(3,j) = yPhi(6 * j + 3,1);
            Phi(4,j) = yPhi(6 * j + 4,1);
            Phi(5,j) = yPhi(6 * j + 5,1);
            Phi(6,j) = yPhi(6 * j + 6,1);
        }

        DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);


        for (int i = 1; i <= 6; i++) {
            Y(i, 1) = Y_old(i, 1);
        }



        //% Topocentric coordinates
        double theta = gmst(Mjd_UT1);                    //% Earth rotation
        Matrix U = R_z(theta);

        Matrix r(3,1);
        for(int i=1;i<=3;i++){
            r(i,1)=Y(i,1);
        }

        Matrix s(3,1);
        s = LT*(U*r-Rs.transpose());                    //% Topocentric position[m]


        //% Time update
        P = TimeUpdate(P, Phi);

        //% Azimuth andpartials
        double Az, El;
        Matrix dAds(1, 3), dEds(1, 3);

        AzElPa(s.transpose(),Az,El,dAds,dEds);     //% Azimuth, Elevation

        Matrix aux(1,3);
        Matrix dAdY(1,6);

        dAdY=Matrix::concatenateHorizontal((dAds * LT * U),aux);


        //% Measurement update
        Matrix K(6,1);
        Matrix obsaux(1,1);
        obsaux(1,1)=obs(i,2);

        Matrix Azaux(1,1);
        Azaux(1,1)=Az;

        Matrix sigma_azaux(1,1);
        sigma_azaux(1,1)=sigma_az;



        MeasUpdate(Y, obsaux, Azaux, sigma_azaux, dAdY, P, 6, K);



        //% Elevation andpartials
        for(int i=1;i<=3;i++){
            r(i,1)=Y(i,1);
        }


        s = LT * (U * r - Rs.transpose());                          //% Topocentric position[m]
        AzElPa(s.transpose(),Az,El,dAds,dEds);     //% Azimuth, Elevation
        Matrix dEdY = Matrix::concatenateHorizontal(dEds * LT * U,aux);


        obsaux(1,1)=obs(i,3);
        Matrix Elevaux(1,1);
        Elevaux(1,1)=El;
        Matrix sigma_elaux(1,1);
        sigma_elaux(1,1)=sigma_el;
        //% Measurement update

        MeasUpdate(Y, obsaux, Elevaux, sigma_elaux, dEdY, P, 6, K);



        //% Range andpartials
        for(int i=1;i<=3;i++){
            r(i,1)=Y(i,1);
        }
        s = LT * (U * r - Rs.transpose());                          //% Topocentric position[m]
        double Dist = s.norm();
        Matrix dDds = (s / Dist).transpose();            //% Range
        Matrix dDdY = Matrix::concatenateHorizontal(dDds * LT * U,aux);


        obsaux(1,1)=obs(i,4);
        Matrix Distaux(1,1);
        Distaux(1,1)=Dist;
        Matrix sigma_raaux(1,1);
        sigma_raaux(1,1)=sigma_range;
        //% Measurement update
        MeasUpdate(Y, obsaux, Distaux, sigma_raaux, dDdY, P, 6,K);



    }




    Matrix iers2 = IERS(*Global::eopdata,obs(46,1),'l');
    Matrix timediff2 = timeDiff(iers2(1,3),iers2(1,9));
    double Mjd_TT = Mjd_UTC + timediff2(1,4)/86400.0;
    Global::AuxParam::Mjd_UTC = Mjd_UTC;
    Global::AuxParam::Mjd_TT = Mjd_TT;

    DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

    double Y_true[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};

    std::cout << "Error of Position Estimation" << std::endl;
    std::cout << "dX  " << Y(1,1) - Y_true[0] << "   [m]" << std::endl;
    std::cout << "dY  " << Y(2,1) - Y_true[1] << "   [m]" << std::endl;
    std::cout << "dZ  " << Y(3,1) - Y_true[2] << "   [m]" << std::endl;
    std::cout << "\nError of Velocity Estimation" << std::endl;
    std::cout << "dVx  " << Y(4,1) - Y_true[3] << "   [m/s]" << std::endl;
    std::cout << "dVy  " << Y(5,1) - Y_true[4] << "   [m/s]" << std::endl;
    std::cout << "dVz  " << Y(6,1) - Y_true[5] << "   [m/s]" << std::endl;




    return 0;
}
