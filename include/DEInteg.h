//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_DEINTEG_H
#define PROYECTO_DEINTEG_H

#include <limits>
#include "Matrix.h"

double DEInteg(double (*function)(double, double), double t, double tout, double relerr, double abserr, int n_eqn, double y);


#endif //PROYECTO_DEINTEG_H
