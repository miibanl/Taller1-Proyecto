//
// Created by miguel on 14/05/2024.
//

#ifndef PROYECTO_TIMEUPDATE_H
#define PROYECTO_TIMEUPDATE_H

#include "Matrix.h"

/*!
 * @file TimeUpdate.h
 * @brief Performs the time update (prediction) step of a Kalman filter
 * @param P The state covariance matrix
 * @param Phi The state transition matrix
 * @param Qdt The process noise covariance, defaulting to 0.0
 * @return The updated state covariance matrix
 */
Matrix TimeUpdate(Matrix& P, Matrix& Phi, double Qdt=0.0);


#endif //PROYECTO_TIMEUPDATE_H
