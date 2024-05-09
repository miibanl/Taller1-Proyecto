//
// Created by miguel on 26/04/2024.
//

#include "../include/Cheb3D.h"
/*
%--------------------------------------------------------------------------
%
% Chebyshev approximation of 3-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
%     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
%     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix Cheb3D(double t, int n, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz){
    // Check validity
    if (t < Ta || Tb < t) {
        throw std::runtime_error("ERROR: Time out of range in Cheb3D::Value");
    }

    // Clenshaw algorithm
    double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1(1, 3), f2(1, 3), old_f1(1, 3), aux(1,3);



    for (int i = n; i >= 2; i--) {
        old_f1 = f1;
        aux(1, 1) = Cx(1, i);
        aux(1, 2) = Cy(1, i);
        aux(1, 3) = Cz(1, i);
        f1 = 2.0 * tau * f1 - f2 + aux;
        f2 = old_f1;
    }

    aux(1, 1) = Cx(1, 1);
    aux(1, 2) = Cy(1, 1);
    aux(1, 3) = Cz(1, 1);

    return tau * f1 - f2 + aux;
}
