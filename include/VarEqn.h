//
// Created by miibanl on 16/05/2024.
//

#ifndef PROYECTO_VAREQN_H
#define PROYECTO_VAREQN_H


#include "Matrix.h"
#include "Global.h"
#include "IERS.h"
#include "TimeDiff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"


Matrix VarEqn(double x,Matrix& yPhi);

#endif //PROYECTO_VAREQN_H
