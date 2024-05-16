//
// Created by miibanl on 24/04/2024.
//


#include "../include/Global.h"

Matrix *Global::eopdata;

void Global::eop19620101(int fila) {
    Global::eopdata = new Matrix(fila, 13);

    FILE *fid = fopen("../data/eop19620101.txt","r");
    if (fid == nullptr) {
        printf("Error al abrir el fichero.");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= fila; i++) {
        fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &((*eopdata)(i, 1)), &((*eopdata)(i, 2)), &((*eopdata)(i, 3)),
               &((*eopdata)(i, 4)), &((*eopdata)(i, 5)), &((*eopdata)(i, 6)),
               &((*eopdata)(i, 7)), &((*eopdata)(i, 8)),
               &((*eopdata)(i, 9)), &((*eopdata)(i, 10)),
               &((*eopdata)(i, 11)), &((*eopdata)(i, 12)),
               &((*eopdata)(i, 13)));
    }

    fclose(fid);
}

Matrix *Global::Cnm;
Matrix *Global::Snm;

void Global::GGM03S () {
    Global::Cnm = new Matrix(181, 181);
    Global::Snm = new Matrix(181, 181);

    Matrix *aux = new Matrix(6, 1);


    FILE *fid = fopen("../data/GGM03S.txt","r");
    if (fid == NULL) {
        printf("Error al abrir el fichero.");
        exit(EXIT_FAILURE);
    }

    for (int n = 0; n <= 180; n++) {
        for (int m = 0; m <= n; m++) {
            fscanf(fid, "%lf %lf %lf %lf %lf %lf", &((*aux)(1,1)), &((*aux)(2,1)), &((*aux)(3,1)),
                   &((*aux)(4,1)), &((*aux)(5,1)), &((*aux)(6,1)));

            (*Global::Cnm)(n + 1, m + 1) = (*aux)(3,1);
            (*Global::Snm)(n + 1, m + 1) = (*aux)(4,1);
        }
    }

    fclose(fid);
}


Matrix *Global::PC;

void Global::DE430Coeff(){
    Global::PC = new Matrix(2285, 1020);


    FILE *fid = fopen("../data/DE430Coeff.txt","r");
    if (fid == nullptr) {
        printf("Error al abrir el fichero.");
        exit(EXIT_FAILURE);
    }

    for(int i=1;i<=2285;i++){
        for(int j=1;j<=1020;j++){
            fscanf(fid, "%lf", &((*PC)(i,j)));
        }
    }


    fclose(fid);
}

Matrix *Global::geos3;

void Global::GEOS3(int fila) {
    Global::geos3 = new Matrix(1, 1);

    FILE *fid = fopen("../data/eop19620101.txt", "r");
    if (fid == nullptr) {
        printf("Error al abrir el fichero.");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= 46; i++) {
        fscanf(fid, "%s",
               &((*geos3)(i,1)));
    }
}



double Global::AuxParam::Mjd_UTC = 0.0;
double Global::AuxParam::Mjd_TT = 0.0;
int Global::AuxParam::n = 0;
int Global::AuxParam::m = 0;
int Global::AuxParam::sun = 0;
int Global::AuxParam::moon = 0;
int Global::AuxParam::planets = 0;

