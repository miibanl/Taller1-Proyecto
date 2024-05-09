//
// Created by miibanl on 24/04/2024.
//

#include "../include/sign.h"

using namespace  std;

double sign_(double a, double b) {
    if (b >= 0.0) {
        return abs(a);
    } else {
        return -abs(a);
    }
}