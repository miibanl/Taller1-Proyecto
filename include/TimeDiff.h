//
// Created by miibanl on 24/04/2024.
//

#ifndef PROYECTO_TIMEDIFF_H
#define PROYECTO_TIMEDIFF_H


#include "Matrix.h"

/*!
 * @file TimeDiff.h
 * @brief Calculates time differences between various time scales
 * @param UT1_UTC Difference between Universal Time 1 (UT1) and Coordinated Universal Time (UTC) in seconds
 * @param TAI_UTC Difference between International Atomic Time (TAI) and Coordinated Universal Time (UTC) in seconds
 * @return A Matrix object containing the calculated time differences
 */
Matrix timeDiff(double UT1_UTC, double TAI_UTC);




#endif //PROYECTO_TIMEDIFF_H
