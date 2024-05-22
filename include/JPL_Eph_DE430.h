//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_JPL_EPH_DE430_H
#define PROYECTO_JPL_EPH_DE430_H

#include "Matrix.h"


/*!
 * @file JPL_Eph_DE430.h
 * @brief Computes the sun, moon, and nine major planets' equatorial position using JPL Ephemerides
 * @details Light-time is already taken into account
 * @param Mjd_TDB Modified julian date of TDB
 * @param r_Mercury geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Venus geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Earth solar system barycenter (SSB)
 * @param r_Mars geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Jupiter geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Saturn geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Uranus geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Neptune geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Pluto geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Moon geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @param r_Sun geocentric equatorial position ([m]) referred to the International Celestial Reference Frame (ICRF)
 * @note The positions are given in meters and referred to the International Celestial Reference Frame (ICRF).
 * @note The input Modified Julian Date of TDB (Mjd_TDB) is required for the calculation.
 * @note This function uses Chebyshev polynomial interpolation to compute the positions of celestial bodies.
 */
void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth,
                     Matrix& r_Mars, Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus,
                     Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun);


#endif //PROYECTO_JPL_EPH_DE430_H
