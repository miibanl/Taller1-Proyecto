//
// Created by miguel on 26/04/2024.
//

#ifndef PROYECTO_MJDAY_H
#define PROYECTO_MJDAY_H

#include <cmath>



/*!
 * @file Mjday.h
 * @brief Computes the Modified Julian Date (MJD) from the given date and time
 * @param yr year
 * @param mon month
 * @param day day
 * @param hr Universal time hour (optional, default is 0)
 * @param min Universal time minute (optional, default is 0)
 * @param sec Universal time second (optional, default is 0)
 * @return Modified julian date
 */
double Mjday(int yr, int mon, int day, int hr=0, int min=0, int sec=0);

#endif //PROYECTO_MJDAY_H
