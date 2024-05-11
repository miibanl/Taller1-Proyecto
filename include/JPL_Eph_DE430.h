//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_JPL_EPH_DE430_H
#define PROYECTO_JPL_EPH_DE430_H

#include "Matrix.h"

void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth,
                     Matrix& r_Mars, Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus,
                     Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun);


#endif //PROYECTO_JPL_EPH_DE430_H
