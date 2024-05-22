//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_EQNEQUINOX_H
#define PROYECTO_EQNEQUINOX_H

#include "Matrix.h"
#include "NutAngles.h"
#include "MeanObliquity.h"


/*!
 * @file EqnEquinox.h
 * @brief Computation of the equation of the equinoxes
 * @details The equation of the equinoxes dpsi*cos(eps) is the right ascension of he mean equinox referred to the true equator and equinox and is equal to the difference between apparent and mean sidereal time.
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Equation of the equinoxes
 */
double EqnEquinox(double Mjd_TT);


#endif //PROYECTO_EQNEQUINOX_H
