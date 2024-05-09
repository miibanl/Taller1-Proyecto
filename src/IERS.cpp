//
// Created by miguel on 26/04/2024.
//

#include "../include/IERS.h"
#include "iostream"
/*
%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%
% Last modified:   2018/02/01   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix IERS(const Matrix& eop, double Mjd_UTC, char interp){

    Matrix results(1, 9);


    if (interp == 'l') {
        // Interpolación lineal
        int mjd = static_cast<int>(floor(Mjd_UTC));
        int i;
        for (i = 1; i <= eop.getRows(); i++) {
            if (mjd == (eop)(i, 4)) {
                break;
            }
        }

        if (i == eop.getRows()+1) {
            throw std::runtime_error("MJD not found in eop.");
        }

        double mfme = 1440.0 * (Mjd_UTC - static_cast<int>(floor(Mjd_UTC)));
        double fixf = mfme / 1440.0;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        double x_pole  = (eop)(i, 5) + ((eop)(i+1, 5) - (eop)(i, 5)) * fixf;
        double y_pole  = (eop)(i, 6) + ((eop)(i+1, 6) - (eop)(i, 6)) * fixf;
        double UT1_UTC = (eop)(i, 7) + ((eop)(i+1, 7) - (eop)(i, 7)) * fixf;
        double LOD     = (eop)(i, 8) + ((eop)(i+1, 8) - (eop)(i, 8)) * fixf;
        double dpsi    = (eop)(i, 9) + ((eop)(i+1, 9) - (eop)(i, 9)) * fixf;
        double deps    = (eop)(i, 10) + ((eop)(i+1, 10) - (eop)(i, 10)) * fixf;
        double dx_pole = (eop)(i, 11) + ((eop)(i+1, 11) - (eop)(i, 11)) * fixf;
        double dy_pole = (eop)(i, 12) + ((eop)(i+1, 12) - (eop)(i, 12)) * fixf;
        double TAI_UTC = (eop)(i, 13);

        x_pole  /= Constants::Arcs;
        y_pole  /= Constants::Arcs;
        dpsi    /= Constants::Arcs;
        deps    /= Constants::Arcs;
        dx_pole /= Constants::Arcs;
        dy_pole /= Constants::Arcs;

        results(1, 1) = x_pole;
        results(1, 2) = y_pole;
        results(1, 3) = UT1_UTC;
        results(1, 4) = LOD;
        results(1, 5) = dpsi;
        results(1, 6) = deps;
        results(1, 7) = dx_pole;
        results(1, 8) = dy_pole;
        results(1, 9) = TAI_UTC;

        return results;
    } else if (interp == 'n') {
        // No hay interpolación
        int mjd = static_cast<int>(floor(Mjd_UTC));
        int i;
        for (i = 1; i <= eop.getRows(); i++) {
            if (mjd == (eop)(i, 4)) {
                break;
            }
        }
        if (i == eop.getRows()+1) {
            throw std::runtime_error("MJD not found in eop.");
        }
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        double x_pole  = (eop)(i, 5) / Constants::Arcs;   // Pole coordinate [rad]
        double y_pole  = (eop)(i, 6) / Constants::Arcs;   // Pole coordinate [rad]
        double UT1_UTC = (eop)(i, 7);                     // UT1-UTC time difference [s]
        double LOD     = (eop)(i, 8);                     // Length of day [s]
        double dpsi    = (eop)(i, 9) / Constants::Arcs;
        double deps    = (eop)(i, 10) / Constants::Arcs;
        double dx_pole = (eop)(i, 11) / Constants::Arcs;  // Pole coordinate [rad]
        double dy_pole = (eop)(i, 12) / Constants::Arcs;  // Pole coordinate [rad]
        double TAI_UTC = (eop)(i, 13);                    // TAI-UTC time difference [s]

        results(1, 1) = x_pole;
        results(1, 2) = y_pole;
        results(1, 3) = UT1_UTC;
        results(1, 4) = LOD;
        results(1, 5) = dpsi;
        results(1, 6) = deps;
        results(1, 7) = dx_pole;
        results(1, 8) = dy_pole;
        results(1, 9) = TAI_UTC;

        return results;
    } else {
        // Interpolación no especificada
        throw std::invalid_argument("Invalid interpolation type specified.");
    }
}

