//
// Created by miguel on 14/05/2024.
//

#ifndef PROYECTO_MEASUPDATE_H
#define PROYECTO_MEASUPDATE_H

#include "Matrix.h"


/*!
 * @file MeasUpdate.h
 * @brief Performs the measurement update step of the Kalman filter
 * @details This function updates the state estimate and covariance matrix of a Kalman filter based on measurement data. It computes the Kalman gain, updates the state estimate, and updates the covariance matrix using the measurement data, predicted measurement, measurement noise covariance, measurement sensitivity matrix, state estimate, and covariance matrix.
 * @param x State estimate (vector) before the measurement update
 * @param z Measurement vector
 * @param g Predicted measurement vector
 * @param s Measurement noise covariance matrix (inverse weight matrix)
 * @param G Measurement sensitivity matrix
 * @param P Covariance matrix of the state estimate before the measurement update
 * @param n Dimension of the state vector
 * @param K Kalman gain matrix (output)
 */
void MeasUpdate(Matrix& x, Matrix& z, Matrix& g, Matrix& s, Matrix& G, Matrix& P, int n, Matrix& K);


#endif //PROYECTO_MEASUPDATE_H
