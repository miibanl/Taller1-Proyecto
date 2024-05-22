//
// Created by miibanl on 24/04/2024.
//

#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "Matrix.h"
#include <cstdio>
#include "cstdlib"
/*!
 * @file Global.h
 * @class Global
 * @brief Provides global variables and functions for managing Earth orientation parameters, gravity model coefficients, planetary coefficients, and GEOS3 data.
 */
class Global {
public:
    // Earth orientation parameters data
    static Matrix* eopdata;
    /*!
     * @brief Loads Earth orientation parameters from a file.
     * @param fila Number of rows to read from the file.
     */
    static void eop19620101(int fila);

    // Gravity model coefficients
    static Matrix* Cnm;
    static Matrix* Snm;
    /*!
     * @brief Loads gravity model coefficients from a file.
     */
    static void GGM03S();

    // Auxiliary parameters for various calculations
    struct AuxParam {
        static double Mjd_UTC;  ///< Modified Julian Date for Coordinated Universal Time
        static double Mjd_TT;   ///< Modified Julian Date for Terrestrial Time
        static int n;           ///< Degree of the spherical harmonics
        static int m;           ///< Order of the spherical harmonics
        static int sun;         ///< Flag to include the Sun
        static int moon;        ///< Flag to include the Moon
        static int planets;     ///< Flag to include other planets
    };

    // Planetary coefficients
    static Matrix* PC;
    /*!
     * @brief Loads planetary coefficients from a file.
     */
    static void DE430Coeff();

    // GEOS3 data
    static Matrix* geos3;
    /*!
     * @brief Loads GEOS3 data from a file.
     * @param fila Number of rows to read from the file.
     */
    static void GEOS3(int fila);
};


#endif //PROYECTO_GLOBAL_H
