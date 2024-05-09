//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H


#include "Matrix.h"
#include <stdexcept>
#include "SAT_Const.h"


Matrix IERS(const Matrix& eop, double Mjd_UTC, char interp= 'n');


#endif //PROYECTO_IERS_H
