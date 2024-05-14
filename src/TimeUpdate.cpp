//
// Created by miguel on 14/05/2024.
//

#include "../include/TimeUpdate.h"



Matrix TimeUpdate(Matrix& P, Matrix& Phi, double Qdt){

    return  Phi*P*Phi.transpose() + Qdt;
}