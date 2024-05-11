//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_ACCEL_H
#define PROYECTO_ACCEL_H


#include "Matrix.h"
#include "../include/Global.h"
#include "../include/IERS.h"
#include "../include/TimeDiff.h"
#include "../include/PrecMatrix.h"





Matrix Accel(double x, const Matrix& Y);


#endif //PROYECTO_ACCEL_H
