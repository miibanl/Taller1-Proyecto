//
// Created by miibanl on 18/04/2024.
//

#ifndef PROYECTO_SAT_CONST_H
#define PROYECTO_SAT_CONST_H


/*!
 * @file SAT_Const.h
 * @file Constants
 * @brief Class containing various physical and mathematical constants
 * @details This class provides static member variables representing various physical and mathematical constants used in astronomical calculations and other scientific applications.
 */
class Constants {
public:
    // Mathematical constants
    static constexpr double pi = 3.141592653589793;            //!< Value of pi (π)
    static constexpr double pi2 = 2 * pi;                     //!< Value of 2π
    static constexpr double Rad = pi / 180;                   //!< Radians per degree
    static constexpr double Deg = 180 / pi;                   //!< Degrees per radian
    static constexpr double Arcs = 3600 * 180 / pi;           //!< Arcseconds per radian

    // General
    static constexpr double MJD_J2000 = 51544.5;              //!< Modified Julian Date of J2000
    static constexpr double T_B1950 = -0.500002108;           //!< Epoch B1950
    static constexpr double c_light = 299792458.000000000;    //!< Speed of light in vacuum [m/s]
    static constexpr double AU = 149597870700.000000;         //!< Astronomical unit [m]

    // Physical parameters of the Earth, Sun, and Moon

    // Equatorial radius and flattening
    static constexpr double R_Earth = 6378.1363e3;            //!< Equatorial radius of the Earth [m]
    static constexpr double f_Earth = 1 / 298.257223563;      //!< Flattening of the Earth (1/f)
    static constexpr double R_Sun = 696000e3;                 //!< Radius of the Sun [m]
    static constexpr double R_Moon = 1738e3;                  //!< Radius of the Moon [m]

    // Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
    static constexpr double omega_Earth = 15.04106717866910 / 3600 * Rad; //!< Earth's rotation rate [rad/s]

    // Gravitational coefficients
    static constexpr double GM_Earth = 398600.435436e9;       //!< Gravitational constant of the Earth [m^3/s^2]
    static constexpr double GM_Sun = 132712440041.939400e9;   //!< Gravitational constant of the Sun [m^3/s^2]
    static constexpr double GM_Moon = GM_Earth / 81.30056907419062; //!< Gravitational constant of the Moon [m^3/s^2]
    static constexpr double GM_Mercury = 22031.780000e9;      //!< Gravitational constant of Mercury [m^3/s^2]
    static constexpr double GM_Venus = 324858.592000e9;       //!< Gravitational constant of Venus [m^3/s^2]
    static constexpr double GM_Mars = 42828.375214e9;         //!< Gravitational constant of Mars [m^3/s^2]
    static constexpr double GM_Jupiter = 126712764.800000e9;  //!< Gravitational constant of Jupiter [m^3/s^2]
    static constexpr double GM_Saturn = 37940585.200000e9;    //!< Gravitational constant of Saturn [m^3/s^2]
    static constexpr double GM_Uranus = 5794548.600000e9;     //!< Gravitational constant of Uranus [m^3/s^2]
    static constexpr double GM_Neptune = 6836527.100580e9;    //!< Gravitational constant of Neptune [m^3/s^2]
    static constexpr double GM_Pluto = 977.0000000000009e9;   //!< Gravitational constant of Pluto [m^3/s^2]

    // Solar radiation pressure at 1 AU
    static constexpr double P_Sol = 1367 / c_light;           //!< Solar radiation pressure at 1 AU [N/m^2]
};


#endif //PROYECTO_SAT_CONST_H
