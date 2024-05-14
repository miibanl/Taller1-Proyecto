//
// Created by miguel on 14/05/2024.
//

#include "../include/angl.h"

/*
%--------------------------------------------------------------------------
%
%  inputs:
%    vec1         - vector 1
%    vec2         - vector 2
%
%  output:
%    theta        - angle between the two vectors  -pi to pi
%
%--------------------------------------------------------------------------
*/
double angl(const Matrix& vec1, const Matrix& vec2){
    double small     = 0.00000001;
    double undefined = 999999.1;
    double temp;

    double magv1 = vec1.norm();
    double magv2 = vec2.norm();

    if (magv1*magv2 > pow(small,2)) {
        temp = Matrix::dot(vec1, vec2) / (magv1 * magv2);
        if (std::abs(temp) > 1.0) {
            temp = ((temp > 0) - (temp < 0)) * 1.0;
        }
        return acos(temp);
    }
    else{
        return undefined;
    }
}
