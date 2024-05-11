//
// Created by miguel on 06/05/2024.
//

#include "../include/AccelHarmonic.h"

Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max) {

    Global::GGM03S();

    double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    double gm = 398600.4415e9; // [m^3/s^2]; GGM03S
    Matrix r_bf = E * r;                 // Body-fixed position

    // Auxiliary quantities
    double d = r_bf.norm();                     // distance
    double latgc = asin(r_bf(1, 3) / d);
    double lon = atan2(r_bf(1, 2), r_bf(1, 1));

    Matrix pnm(n_max + 1, m_max + 1);
    Matrix dpnm(n_max + 1, m_max + 1);

    Legendre(n_max, m_max, latgc, pnm, dpnm);

    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0, q2 = q3, q1 = q2;
    for (int n = 0; n <= n_max; n++) {
        double b1 = (-gm / (d * d)) * pow((r_ref / d), n) * (n + 1);
        double b2 = (gm / d) * pow((r_ref / d), n);
        double b3 = (gm / d) * pow((r_ref / d), n);

        for (int m = 0; m <= m_max; ++m) {
                q1 = q1 + pnm(n + 1, m + 1) * ((*Global::Cnm)(n + 1, m + 1) * cos(m * lon) + (*Global::Snm)(n + 1, m + 1) * sin(m * lon));
                q2 = q2 + dpnm(n + 1, m + 1) * ((*Global::Cnm)(n + 1, m + 1) * cos(m * lon) + (*Global::Snm)(n + 1, m + 1) * sin(m * lon));
                q3 = q3 + m * pnm(n + 1, m + 1) * ((*Global::Snm)(n + 1, m + 1) * cos(m * lon) - (*Global::Cnm)(n + 1, m + 1) * sin(m * lon));
        }
        dUdr = dUdr + q1 * b1;
        dUdlatgc = dUdlatgc + q2 * b2;
        dUdlon = dUdlon + q3 * b3;
        q3 = 0;
        q2 = q3;
        q1 = q2;
    }

    //% Body-fixed acceleration
    double r2xy = pow(r_bf(1, 1), 2) + pow(r_bf(1, 2), 2);

    double ax = (1.0 / d * dUdr - r_bf(1, 3) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) - (1.0 / r2xy * dUdlon) * r_bf(1, 2);
    double ay = (1.0 / d * dUdr - r_bf(1, 3) / (d * d * sqrt(r2xy)) * dUdlatgc) * r_bf(1, 2) + (1.0 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1.0 / d * dUdr * r_bf(1, 3) + sqrt(r2xy) / (d * d) * dUdlatgc;



    Matrix a_bf(1,3);
    a_bf(1,1)=ax;
    a_bf(1,2)=ay;
    a_bf(1,3)=az;

                     // Inertial acceleration
    return  E*a_bf;
}
