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

    Matrix a= P*G.transpose();

    Matrix aux(6,1);
    for (int i = 1; i<=6 ; ++i) {
        aux(i,1)=a(i,1);
    }

    Matrix b= (Inv_W+G*P*G.transpose()).inverse();


    //% Kalman gain
    if(b.getRows()==1 && b.getCols()==1) {
        for (int i = 1; i <= aux.getRows(); i++) {
            K(i, 1) = aux(i, 1) * b(1, 1);
        }
    }else{
        K = aux*b;
    }

    //% State update
    x = x + K*(z-g);

    Matrix auxK(K.getRows(),K.getRows());
    Matrix auxG(G.getCols(),G.getCols());

    for (int i = 1; i <= K.getRows(); ++i) {
        auxK(i,1)=K(i,1);
    }
    for (int i = 1; i <= G.getCols(); ++i) {
        auxG(1,i)=G(1,i);
    }


    Matrix aux2=auxK*auxG;
    //% Covariance update
    P = (Matrix::createIdentityMatrix(n)-aux2)*P;
    


}