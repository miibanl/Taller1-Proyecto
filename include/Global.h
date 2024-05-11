//
// Created by miibanl on 24/04/2024.
//

#ifndef PROYECTO_GLOBAL_H
#define PROYECTO_GLOBAL_H

#include "Matrix.h"
#include <cstdio>
#include "cstdlib"

class Global{
public:
    static Matrix *eopdata;
    static void eop19620101(int fila);

    static Matrix *Cnm;
    static Matrix *Snm;
    static void GGM03S();

    struct AuxParam {
        static double Mjd_UTC;
        static int n;
        static int m;
        static int sun;
        static int moon;
        static int planets;
    };

    static Matrix *PC;
    static void DE430Coeff();

    static Matrix *geos3;
    static void GEOS3(int fila);




};


#endif //PROYECTO_GLOBAL_H
