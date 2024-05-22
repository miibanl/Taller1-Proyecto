//
// Created by miguel on 11/05/2024.
//

#ifndef PROYECTO_POLEMATRIX_H
#define PROYECTO_POLEMATRIX_H

#include "Matrix.h"
#include "R_y.h"
#include "R_x.h"



/*!
 * @file PoleMatrix.h
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 * @param xp Pole coordinte
 * @param yp Pole coordinte
 * @return Pole matrix
 */
Matrix PoleMatrix(double xp, double yp);

#endif //PROYECTO_POLEMATRIX_H
