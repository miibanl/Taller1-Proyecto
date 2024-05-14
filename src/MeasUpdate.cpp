//
// Created by miguel on 14/05/2024.
//

#include "../include/MeasUpdate.h"


void MeasUpdate(Matrix& x, Matrix& z, Matrix& g, Matrix& s, Matrix& G, Matrix& P, int n, Matrix& K){
    int m = z.getCols();
    Matrix Inv_W(m,m);

    for(int i=1;i<=m;i++) {
        Inv_W(i, i) = s(1,i) * s(1,i);    //% Inverse weight(measurement covariance)
    }

    //% Kalman gain
    K = P*G.transpose()*(Inv_W+G*P*G.transpose()).inverse();

    //% State update
    x = x + K*(z-g);

    //% Covariance update
    P = (Matrix::createIdentityMatrix(n)-K*G)*P;
}