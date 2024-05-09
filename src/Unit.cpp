//
// Created by miibanl on 24/04/2024.
//

#include "../include/Unit.h"


Matrix unit(const Matrix& vec){
    double small = 0.000001;
    double magv = vec.norm();
    Matrix outvec(1,3);

    if (magv > small) {
        for(int i = 1; i <= 3; i++) {
            outvec(1, i) = (vec(1, i) / magv);
        }
    } else {
        for(int i = 1; i <= 3; i++) {
            outvec(1, i) = 0.0;
        }
    }
    return outvec;
}