//
// Created by miguel on 14/05/2024.
//

#ifndef PROYECTO_ELEMENTS_H
#define PROYECTO_ELEMENTS_H

#include "Matrix.h"
#include "SAT_Const.h"


/*!
 * @file elements.h
 * @brief The function cannot be used with state vectors describing a circular or non-inclined orbit
 * @param y State vector (x,y,z,vx,vy,vz)
 * @return Matriz with(semilatus rectum [m],semilatus rectum [m],Eccentricity,Inclination [rad],Longitude of the ascending node [rad],Argument of pericenter [rad],Mean anomaly [rad])
 */
Matrix elements(const Matrix& y);


#endif //PROYECTO_ELEMENTS_H
