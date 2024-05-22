//
// Created by miguel on 14/05/2024.
//

#ifndef PROYECTO_ANGL_H
#define PROYECTO_ANGL_H

#include "Matrix.h"


/*!
 * @file angl.h
 * @brief Computes the angle between two vectors
 * @param vec1 vector 1
 * @param vec2 vector 2
 * @return angle between the two vectors  -pi to pi
 */
double angl(const Matrix& vec1, const Matrix& vec2);

#endif //PROYECTO_ANGL_H
