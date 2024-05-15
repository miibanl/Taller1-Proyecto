#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include "./include/Matrix.h"
#include "./include/R_x.h"
#include "./include/R_y.h"
#include "./include/R_z.h"
#include "./include/SAT_Const.h"
#include "./include/Global.h"
#include "./include/sign.h"
#include "./include/TimeDiff.h"
#include "./include/Unit.h"
#include "./include/AccelPointMass.h"
#include "./include/AzElPa.h"
#include "./include/Cheb3D.h"
#include "./include/EccAnom.h"
#include "./include/Frac.h"
#include "./include/Geodetic.h"
#include "./include/Legendre.h"
#include "./include/MeanObliquity.h"
#include "./include/Mjday.h"
#include "./include/Mjday_TDB.h"
#include "./include/NutAngles.h"
#include "./include/Position.h"
#include "./include/IERS.h"
#include "./include/PoleMatrix.h"
#include "./include/NutMatrix.h"
#include "./include/PrecMatrix.h"
#include "./include/gmst.h"
#include "./include/EqnEquinox.h"
#include "./include/gast.h"
#include "./include/GHAMatrix.h"
#include "./include/AccelHarmonic.h"
#include "./include/LTC.h"
#include "./include/elements.h"
#include "./include/angl.h"
#include "./include/TimeUpdate.h"
#include "./include/MeasUpdate.h"
#include "./include/JPL_Eph_DE430.h"
















int tests_run = 0;

#define TOL_ 10e-10
#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;





int proMat_01()
{
    double v1[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {1.0, 0.0, 0.0, 1.0};
    Matrix m1(2, 2, v1, 4);
    Matrix m2(2, 2, v2, 4);
    Matrix sol(2, 2);


    sol = m1 * m2;

    //m1.print();
    //m2.print();
    //sol.print();

    _assert(sol(1,1) == m1(1,1) && sol(1,2) == m1(1,2) && sol(2,1) == m1(2,1) && sol(2,2) == m1(2,2));
    
    return 0;
}

int R_x_01() {
    Matrix sol(3, 3);

    sol = R_x(1.0);

    _assert(fabs(sol(1, 1) - 1) < TOL_);
    _assert(fabs(sol(1, 2) - 0) < TOL_);
    _assert(fabs(sol(1, 3) - 0) < TOL_);
    _assert(fabs(sol(2, 1) - 0) < TOL_);
    _assert(fabs(sol(2, 2) - 0.540302305868140)  < TOL_);
    _assert(fabs(sol(2, 3) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(3, 1) - 0) < TOL_);
    _assert(fabs(sol(3, 2) - -0.841470984807897) < TOL_);
    _assert(fabs(sol(3, 3) - 0.540302305868140) < TOL_);

    return 0;
}

int R_y_01() {
    Matrix sol(3, 3);

    sol = R_y(1.0);

    _assert(fabs(sol(1, 1) - 0.540302305868140) < TOL_);
    _assert(fabs(sol(1, 2) - 0) < TOL_);
    _assert(fabs(sol(1, 3) - -0.841470984807897) < TOL_);
    _assert(fabs(sol(2, 1) - 0) < TOL_);
    _assert(fabs(sol(2, 2) - 1.0)  < TOL_);
    _assert(fabs(sol(2, 3) - 0) < TOL_);
    _assert(fabs(sol(3, 1) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(3, 2) - 0) < TOL_);
    _assert(fabs(sol(3, 3) - 0.540302305868140) < TOL_);

    return 0;
}


int R_z_01() {
    Matrix sol(3, 3);

    sol = R_z(1.0);

    _assert(fabs(sol(1, 1) - 0.540302305868140) < TOL_);
    _assert(fabs(sol(1, 2) - 0.841470984807897) < TOL_);
    _assert(fabs(sol(1, 3) - 0) < TOL_);
    _assert(fabs(sol(2, 1) - -0.841470984807897) < TOL_);
    _assert(fabs(sol(2, 2) - 0.540302305868140)  < TOL_);
    _assert(fabs(sol(2, 3) - 0) < TOL_);
    _assert(fabs(sol(3, 1) - 0) < TOL_);
    _assert(fabs(sol(3, 2) - 0) < TOL_);
    _assert(fabs(sol(3, 3) - 1.0) < TOL_);

    return 0;
}


int SAT_Const_01() {
    Constants constants;

    _assert(fabs(constants.pi - 3.141592653589793) < TOL_);

    return 0;
}

int Sign_01() {

    _assert(sign_(5.0, 2.0) == 5.0);
    _assert(sign_(-5.0, 2.0) == 5.0);
    _assert(sign_(5.0, -2.0) == -5.0);
    _assert(sign_(-5.0, -2.0) == -5.0);
    _assert(sign_(0.0, 2.0) == 0.0);
    _assert(sign_(5.0, 0.0) == 5.0);

    return 0;
}

int TimeDiff_01() {

    Matrix sol(1, 5);

    sol = timeDiff(10,20);

    _assert(fabs(sol(1, 1) - -10) < TOL_);
    _assert(fabs(sol(1, 2) - -1) < TOL_);
    _assert(fabs(sol(1, 3) - 9) < TOL_);
    _assert(fabs(sol(1, 4) - 52.1840) < TOL_);
    _assert(fabs(sol(1, 5) - 1) < TOL_);

    return 0;
}

int Unit_01() {
    Matrix vec1(1, 3);
    vec1(1, 1) = 1.0;
    vec1(1, 2) = 2.0;
    vec1(1, 3) = 3.0;

    Matrix outvec1 = unit(vec1);
    _assert(fabs(outvec1(1, 1) - 0.267261241912424) < TOL_);
    _assert(fabs(outvec1(1, 2) - 0.534522483824849) < TOL_);
    _assert(fabs(outvec1(1, 3) - 0.801783725737273) < TOL_);

    return 0;
}

int AccelPointMassTest() {
    Matrix r(1, 3);
    r(1, 1) = 1.0;
    r(1, 2) = 2.0;
    r(1, 3) = 3.0;

    Matrix s(1, 3);
    s(1, 1) = 10.0;
    s(1, 2) = 10.0;
    s(1, 3) = 10.0;

    Matrix accel = AccelPointMass(r, s, 1.0);

    _assert(fabs(accel(1, 1) - 0.001406232828247) < TOL_);
    _assert(fabs(accel(1, 2) - 0.001036151303187) < TOL_);
    _assert(fabs(accel(1, 3) - 0.000666069778126) < TOL_);

    return 0;
}

int AzElPa() {
    Matrix s(3, 1);
    s(1, 1) = 1.0;
    s(1, 2) = 1.0;
    s(1, 3) = 1.0;

    double Az, El;
    Matrix dAds(1, 3), dEds(1, 3);

    AzElPa(s, Az, El, dAds, dEds);

    _assert(fabs(Az - 0.785398163397448) < TOL_);
    _assert(fabs(El - 0.615479708670387) < TOL_);

    _assert(fabs(dAds(1, 1) - 0.5) < TOL_);
    _assert(fabs(dAds(1, 2) - -0.5) < TOL_);
    _assert(fabs(dAds(1, 3) - 0.0) < TOL_);

    _assert(fabs(dEds(1, 1) - -0.235702260395516) < TOL_);
    _assert(fabs(dEds(1, 2) - -0.235702260395516) < TOL_);
    _assert(fabs(dEds(1, 3) - 0.471404520791032) < TOL_);

    return 0;
}

int Cheb3D() {

    Matrix Cx(1, 3), Cy(1, 3), Cz(1,3), sol(1,3);

    Cx(1, 1) = 1.0;
    Cx(1, 2) = 0.5;
    Cx(1, 3) = 0.2;
    Cy(1, 1) = 0.8;
    Cy(1, 2) = 0.4;
    Cy(1, 3) = 0.1;
    Cz(1, 1) = 0.6;
    Cz(1, 2) = 0.3;
    Cz(1, 3) = 0.15;

    sol=Cheb3D(0.5, 3, 0, 1, Cx, Cy, Cz);

    _assert(fabs(sol(1, 1) - 0.8) < TOL_);
    _assert(fabs(sol(1, 2) - 0.7) < TOL_);
    _assert(fabs(sol(1, 3) - 0.45) < TOL_);


    return 0;
}

int EccAnom(){

    _assert(fabs(EccAnom(1.5,0.5) - 1.962189287578571) < TOL_);
    _assert(fabs(EccAnom(56.12,3.24) - 3.820352946837100) < TOL_);


    return 0;
}

int Frac(){
    _assert(fabs(Frac(10.256) - 0.256) < TOL_);
    _assert(fabs(Frac(0.0023654) - 0.0023654) < TOL_);

    return 0;
}

int Geodetic(){

    Matrix input(1,3), sol(1,3);

    input(1,1)=6378137;
    input(1,2)=0;
    input(1,3)=0;

    sol=Geodetic(input);

    _assert(fabs(sol(1, 1) - 0.0) < TOL_);
    _assert(fabs(sol(1, 2) - 0.0) < TOL_);
    _assert(fabs(sol(1, 3) - 0.7) < TOL_);

    input(1,1)=6378137;
    input(1,2)=4564;
    input(1,3)=1452;

    sol=Geodetic(input);

    _assert(fabs(sol(1, 1) - 7.155693302005374e-04) < TOL_);
    _assert(fabs(sol(1, 2) - 2.291868841276059e-04) < TOL_);
    _assert(fabs(sol(1, 3) - 2.499318960122764) < TOL_);

    return 0;
}

int Legendre(){

    int n=2,m=2;

    Matrix pnm(n + 1, m + 1);
    Matrix dpnm(n + 1, m + 1);

    Legendre(n,m,Constants::pi/4,pnm,dpnm);


    _assert(fabs(pnm(1, 1) - 1.0) < TOL_);
    _assert(fabs(pnm(1, 2) - 0.0) < TOL_);
    _assert(fabs(pnm(1, 3) - 0.0) < TOL_);
    _assert(fabs(pnm(2, 1) - 1.224744871391589) < TOL_);
    _assert(fabs(pnm(2, 2) - 1.224744871391589) < TOL_);
    _assert(fabs(pnm(2, 3) - 0.0) < TOL_);
    _assert(fabs(pnm(3, 1) - 0.559016994374947) < TOL_);
    _assert(fabs(pnm(3, 2) - 1.936491673103709) < TOL_);
    _assert(fabs(pnm(3, 3) - 0.968245836551854) < TOL_);


    _assert(fabs(dpnm(1, 1) - 0.0) < TOL_);
    _assert(fabs(dpnm(1, 2) - 0.0) < TOL_);
    _assert(fabs(dpnm(1, 3) - 0.0) < TOL_);
    _assert(fabs(dpnm(2, 1) - 1.224744871391589) < TOL_);
    _assert(fabs(dpnm(2, 2) - -1.224744871391589) < TOL_);
    _assert(fabs(dpnm(2, 3) - 0.0) < TOL_);
    _assert(fabs(dpnm(3, 1) - 3.354101966249685) < TOL_);
    _assert(fabs(dpnm(3, 2) - 0.0) < TOL_);
    _assert(fabs(dpnm(3, 3) - -1.936491673103709) < TOL_);


    return 0;
}

int MeanObliquity(){
    _assert(fabs(MeanObliquity(58849.5) - 0.409047411073268) < TOL_);
    _assert(fabs(MeanObliquity(30000000) - 5.066422195176343) < TOL_);

    return 0;
}

int Mjday(){
    _assert((Mjday(1590,12,12,12,25,25) - -9.786248234953685e+04) < TOL_);
    _assert(fabs(Mjday(2024,4,27) - 60427) < TOL_);

    return 0;
}

int Mjday_TDB(){
    _assert(fabs(Mjday_TDB(59580) - 5.957999999999870e+04) < TOL_);
    _assert(fabs(Mjday_TDB(2451.54) - 2.451539999990105e+03) < TOL_);

    return 0;
}

int NutAngles() {
    _assert(fabs(NutAngles(10.652)(1,1) - 3.1820e-05) < TOL_);
    _assert(fabs(NutAngles(10.652)(1,2) - 3.8395e-05) < TOL_);

    return 0;
}

int Position() {
    Matrix sol(1,3);
    sol=Position(-1.3119, 0.6973, 0);

    _assert(fabs(sol(1,1) - 1.253470895636730e+06) < TOL_);
    _assert(fabs(sol(1,2) -  -4.732934494089562e+06) < TOL_);
    _assert(fabs(sol(1,3) - 4.073930489273672e+06) < TOL_);

    return 0;
}

int IERS() {
    Matrix sol(1,9);

    Global::eop19620101(21413);

    //Global::eopdata->print();

    sol=IERS(*Global::eopdata, 49746.0, 'l');

    _assert(fabs(sol(1,1) - -5.59518621231704e-07) < TOL_);
    _assert(fabs(sol(1,2) -  2.33458634442529e-06) < TOL_);
    _assert(fabs(sol(1,3) -  0.3260677) < TOL_);
    _assert(fabs(sol(1,4) -  0.0027213) < TOL_);
    _assert(fabs(sol(1,5) -  -1.16864337831454e-07) < TOL_);
    _assert(fabs(sol(1,6) -  -2.48709418409192e-08) < TOL_);
    _assert(fabs(sol(1,7) -  -8.19335121075116e-10) < TOL_);
    _assert(fabs(sol(1,8) -  -1.53201123230613e-09) < TOL_);
    _assert(fabs(sol(1,9) -  29) < TOL_);

    sol=IERS(*Global::eopdata, 49665.0);

    _assert(fabs(sol(1,1) - -7.00420021372569e-07) < TOL_);
    _assert(fabs(sol(1,2) -  1.34009773663892e-06) < TOL_);
    _assert(fabs(sol(1,3) -  0.5309998) < TOL_);
    _assert(fabs(sol(1,4) -  0.0023392) < TOL_);
    _assert(fabs(sol(1,5) -  -1.24980118853227e-07) < TOL_);
    _assert(fabs(sol(1,6) -  -2.97530156096922e-08) < TOL_);
    _assert(fabs(sol(1,7) -  2.37558703743673e-10) < TOL_);
    _assert(fabs(sol(1,8) -  -9.21145994108118e-11) < TOL_);
    _assert(fabs(sol(1,9) -  29) < TOL_);



    return 0;
}

int PoleMatrix() {
    Matrix sol(3,3);
    sol=PoleMatrix(0.2551451,4515);

    _assert(fabs(sol(1,1) - 0.967626684693315) < TOL_);
    _assert(fabs(sol(1,2) -  -0.127884781524882) < TOL_);
    _assert(fabs(sol(1,3) -  -0.217586952099058) < TOL_);
    _assert(fabs(sol(2,1) -  0) < TOL_);
    _assert(fabs(sol(2,2) -  -0.862120373238167) < TOL_);
    _assert(fabs(sol(2,3) -  0.506703524802901) < TOL_);
    _assert(fabs(sol(3,1) -  -0.252385813922701) < TOL_);
    _assert(fabs(sol(3,2) -  -0.490299851827448) < TOL_);
    _assert(fabs(sol(3,3) -  -0.834210678563011) < TOL_);

    return 0;
}

int PrecMatrix() {
    Matrix sol(3,3);
    sol=PrecMatrix(0.2551451,4515);

    _assert(fabs(sol(1,1) - 0.999995464055671) < TOL_);
    _assert(fabs(sol(1,2) -  -0.00276180170860104) < TOL_);
    _assert(fabs(sol(1,3) -  -0.00120179840492362) < TOL_);
    _assert(fabs(sol(2,1) -  0.00276180170855866) < TOL_);
    _assert(fabs(sol(2,2) -  0.999996186217012) < TOL_);
    _assert(fabs(sol(2,3) -  -1.65960347120066e-06) < TOL_);
    _assert(fabs(sol(3,1) -  0.00120179840502101) < TOL_);
    _assert(fabs(sol(3,2) -  -1.65953294472943e-06) < TOL_);
    _assert(fabs(sol(3,3) -  0.999999277838659) < TOL_);

    return 0;
}

int NutMatrix() {
    Matrix sol(3,3);
    sol=NutMatrix(4515.25);

    _assert(fabs(sol(1,1) - 0.999999996447048) < TOL_);
    _assert(fabs(sol(1,2) -  7.73307456770637e-05) < TOL_);
    _assert(fabs(sol(1,3) -  3.35538296369724e-05) < TOL_);
    _assert(fabs(sol(2,1) -  -7.73307990794727e-05) < TOL_);
    _assert(fabs(sol(2,2) -  0.999999997008709) < TOL_);
    _assert(fabs(sol(2,3) -  1.5902499256315e-06) < TOL_);
    _assert(fabs(sol(3,1) -  -3.35537065613906e-05) < TOL_);
    _assert(fabs(sol(3,2) -  -1.5928446644442e-06) < TOL_);
    _assert(fabs(sol(3,3) -  0.999999999435806) < TOL_);

    return 0;
}

int gmst() {

    _assert(fabs(gmst(6523.25) - 1.66477739228949) < TOL_);
    _assert(fabs(gmst(0.652325) -  5.08310881449846) < TOL_);
    _assert(fabs(gmst(123.456) -  5.96212844883578) < TOL_);

    return 0;
}

int EqnEquinox() {

    _assert(fabs(EqnEquinox(6523.25) - 1.32310916096489e-05) < TOL_);
    _assert(fabs(EqnEquinox(0.652325) -  2.52749545765384e-05) < TOL_);
    _assert(fabs(EqnEquinox(123.456) -  4.00438480107561e-05) < TOL_);

    return 0;
}

int gast() {

    _assert(fabs(gast(6523.25) - 1.6647906233811) < TOL_);
    _assert(fabs(gast(0.652325) -  5.08313408945303) < TOL_);
    _assert(fabs(gast(123.456) -  5.96216849268379) < TOL_);

    return 0;
}

int GHAMatrix() {
    Matrix sol(3,3);
    sol=GHAMatrix(4515.25);

    _assert(fabs(sol(1,1) - 0.108006815576271) < TOL_);
    _assert(fabs(sol(1,2) -  -0.99415015354275) < TOL_);
    _assert(fabs(sol(1,3) -  0) < TOL_);
    _assert(fabs(sol(2,1) -  0.99415015354275) < TOL_);
    _assert(fabs(sol(2,2) -  0.108006815576271) < TOL_);
    _assert(fabs(sol(2,3) -  0) < TOL_);
    _assert(fabs(sol(3,1) -  0) < TOL_);
    _assert(fabs(sol(3,2) -  0) < TOL_);
    _assert(fabs(sol(3,3) -  1) < TOL_);

    return 0;
}


int AccelHarmonic() {/*
    Matrix sol(1,3);

    Matrix r(1,3);
    r(1,1)=7000000;
    r(1,2)=0;
    r(1,3)=0;

    Matrix E(3,3);
    E(1,1)=1;
    E(1,2)=0;
    E(1,3)=0;
    E(2,1)=0;
    E(2,2)=1;
    E(2,3)=0;
    E(3,1)=0;
    E(3,2)=0;
    E(3,3)=1;


    sol=AccelHarmonic(r,E,2,2);

    _assert(fabs(sol(1,1) - -8.14576607065686) < TOL_);
    _assert(fabs(sol(1,2) -  -3.66267894892037e-05) < TOL_);
    _assert(fabs(sol(1,3) -  -5.84508413583961e-09) < TOL_);
*/

    return 0;
}

int LTC() {
    Matrix sol(3,3);
    sol=LTC(0.785398163397448,0.523598775598299);

    _assert(fabs(sol(1,1) - -0.707106781186547) < TOL_);
    _assert(fabs(sol(1,2) -  0.707106781186548) < TOL_);
    _assert(fabs(sol(1,3) -  0) < TOL_);
    _assert(fabs(sol(2,1) -  -0.353553390593274) < TOL_);
    _assert(fabs(sol(2,2) -  -0.353553390593274) < TOL_);
    _assert(fabs(sol(2,3) -  0.866025403784439) < TOL_);
    _assert(fabs(sol(3,1) -  0.612372435695795) < TOL_);
    _assert(fabs(sol(3,2) -  0.612372435695794) < TOL_);
    _assert(fabs(sol(3,3) -  0.5) < TOL_);

    return 0;
}

int elements() {
    Matrix y(1,6);
    y(1,1)=7000;
    y(1,2)=5000;
    y(1,3)=10;
    y(1,4)=0;
    y(1,5)=7000;
    y(1,6)=50;

    Matrix sol(1,7);
    sol=elements(y);

    _assert(fabs(sol(1,1) - 6.02396456836168) < TOL_);
    _assert(fabs(sol(1,2) -  4303.44106929666) < TOL_);
    _assert(fabs(sol(1,3) -  0.999299853996234) < TOL_);
    _assert(fabs(sol(1,4) -  0.00803193693149894) < TOL_);
    _assert(fabs(sol(1,5) -  0.475010756265097) < TOL_);
    _assert(fabs(sol(1,6) -  3.28733653782776) < TOL_);
    _assert(fabs(sol(1,7) -  3.08812323724874) < TOL_);


    return 0;
}

int TimeUpdate() {
    Matrix P(2,2);
    P(1,1)=1;
    P(1,2)=2;
    P(2,1)=3;
    P(2,2)=4;

    Matrix Phi(2,2);
    Phi(1,1)=2;
    Phi(1,2)=1;
    Phi(2,1)=1;
    Phi(2,2)=2;

    Matrix sol(2,2);
    sol = TimeUpdate(P,Phi,0.1);

    _assert(fabs(sol(1,1) - 18.1) < TOL_);
    _assert(fabs(sol(1,2) - 21.1) < TOL_);
    _assert(fabs(sol(2,1) - 24.1) < TOL_);
    _assert(fabs(sol(2,2) - 27.1) < TOL_);

    return 0;
}

int angl() {
    Matrix vec1(1,3);
    vec1(1,1)=1;
    vec1(1,2)=0;
    vec1(1,3)=0;

    Matrix vec2(1,3);
    vec2(1,1)=0;
    vec2(1,2)=1;
    vec2(1,3)=0;

    _assert(fabs(angl(vec1,vec2) - 1.5707963267949) < TOL_);

    return 0;
}

int JPL_Eph_DE430() {

    double Mjd_TDB=49746.1166215902;

    Matrix r_Mercury(3,1);
    Matrix r_Venus(3,1);
    Matrix r_Earth(3,1);
    Matrix r_Mars(3,1);
    Matrix r_Jupiter(3,1);
    Matrix r_Saturn(3,1);
    Matrix r_Uranus(3,1);
    Matrix r_Neptune(3,1);
    Matrix r_Pluto(3,1);
    Matrix r_Moon(3,1);
    Matrix r_Sun(3,1);

    JPL_Eph_DE430(Mjd_TDB,r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun);


    _assert(fabs(r_Mercury(1,1) - 83779075648.5297) < TOL_);
    _assert(fabs(r_Mercury(2,1) - -65292051794.8272) < TOL_);
    _assert(fabs(r_Mercury(3,1) - -23392255339.1343) < TOL_);

    _assert(fabs(r_Venus(1,1) - -15232231703.7507) < TOL_);
    _assert(fabs(r_Venus(2,1) - -110133428027.475) < TOL_);
    _assert(fabs(r_Venus(3,1) - -41021066115.0631) < TOL_);

    _assert(fabs(r_Earth(1,1) - -92468458237.547) < TOL_);
    _assert(fabs(r_Earth(2,1) - 106396734968.241) < TOL_);
    _assert(fabs(r_Earth(3,1) - 46130927542.9284) < TOL_);

    _assert(fabs(r_Mars(1,1) - -88279262157.0358) < TOL_);
    _assert(fabs(r_Mars(2,1) - 46964465202.8572) < TOL_);
    _assert(fabs(r_Mars(3,1) - 29070887645.1543) < TOL_);

    _assert(fabs(r_Jupiter(1,1) - -298389625469.917) < TOL_);
    _assert(fabs(r_Jupiter(2,1) - -754499528838.573) < TOL_);
    _assert(fabs(r_Jupiter(3,1) - -314411042989.654) < TOL_);

    _assert(fabs(r_Saturn(1,1) - 1482031269946.56) < TOL_);
    _assert(fabs(r_Saturn(2,1) - -453875617598.378) < TOL_);
    _assert(fabs(r_Saturn(3,1) - -249403400202.913) < TOL_);

    _assert(fabs(r_Uranus(1,1) - 1412364844217.22) < TOL_);
    _assert(fabs(r_Uranus(2,1) - -2511357129966.91) < TOL_);
    _assert(fabs(r_Uranus(3,1) - -1118109547729.56) < TOL_);

    _assert(fabs(r_Neptune(1,1) - 1871247743893.51) < TOL_);
    _assert(fabs(r_Neptune(2,1) - -3928978347157.18) < TOL_);
    _assert(fabs(r_Neptune(3,1) - -1655021340136.64) < TOL_);

    _assert(fabs(r_Pluto(1,1) - -2171417806084.57) < TOL_);
    _assert(fabs(r_Pluto(2,1) - -3915434641172.15) < TOL_);
    _assert(fabs(r_Pluto(3,1) - -552716789894.163) < TOL_);

    _assert(fabs(r_Moon(1,1) - 89273255.8623995) < TOL_);
    _assert(fabs(r_Moon(2,1) - -336626237.822431) < TOL_);
    _assert(fabs(r_Moon(3,1) - -114663430.62365) < TOL_);

    _assert(fabs(r_Sun(1,1) - 92295749979.7609) < TOL_);
    _assert(fabs(r_Sun(2,1) - -105377012912.08) < TOL_);
    _assert(fabs(r_Sun(3,1) - -45687155004.5967) < TOL_);


    return 0;
}

int MeasUpdate() {
    Matrix x(6,1);
    x(1,1)=7101597.83995656;
    x(2,1)=1295244.79268408;
    x(3,1)=12755.6074812152;
    x(4,1)=576.097761535079;
    x(5,1)=-3084.51250378929;
    x(6,1)=-6736.03539066925;

    Matrix z(1,1);
    x(1,1)=2653472;

    Matrix g(1,1);
    g(1,1)=2653524.97225556;

    Matrix s(1,1);
    s(1,1)=92.5;

    Matrix G(1,6);
    G(1,1)=0.484876050827302;
    G(1,2)=0.042004563746839;
    G(1,3)=-0.873573598478432;
    G(1,4)=0.0;
    G(1,5)=0.0;
    G(1,6)=0.0;

    Matrix P(6,6);
    P(1,1)=15877.8629772695;
    P(1,2)=-5671.36398805618;
    P(1,3)=8540.82446015526;
    P(1,4)=48.5595982378893;
    P(1,5)=-13.3700368251816;
    P(1,6)=22.4252217502378;

    P(2,1)=-5671.36398805618;
    P(2,2)=24467.3370609133;
    P(2,3)=-1495.7021644119;
    P(2,4)=-4.64148177308955;
    P(2,5)=59.550440390861;
    P(2,6)=25850.3882249601;

    P(3,1)=8540.82446015527;
    P(3,2)=-1495.7021644119;
    P(3,3)=6196.68034436514;
    P(3,4)=26.7016364847452;
    P(3,5)=-4.13117368629972;
    P(3,6)=16.0867101164779;

    P(4,1)=48.5595982378894;
    P(4,2)=-4.64148177308954;
    P(4,3)= 26.7016364847452;
    P(4,4)=0.178295954712907;
    P(4,5)=-0.0189362504486505;
    P(4,6)=0.0668595526152223;

    P(5,1)=-13.3700368251816;
    P(5,2)=59.550440390861;
    P(5,3)=-4.1311736862997;
    P(5,4)=-0.0189362504486506;
    P(5,5)=0.152611171443287;
    P(5,6)=-0.0750233063279036;

    P(6,1)=22.4252217502378;
    P(6,2)=-26.9419372559196;
    P(6,3)=16.0867101164779;
    P(6,4)=0.0668595526152224;
    P(6,5)=-0.0750233063279037;
    P(6,6)=0.0845546951479235;




    Matrix K(1,6);

    MeasUpdate(x,z,g,s,G,P,6,K);
    //K.print();


    //_assert(fabs(angl(vec1,vec2) - 1.5707963267949) < TOL_);



    return 0;
}



int all_tests()
{
    _verify(proMat_01);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(SAT_Const_01);
    _verify(Sign_01);
    _verify(TimeDiff_01);
    _verify(Unit_01);
    _verify(AccelPointMassTest);
    _verify(AzElPa);
    _verify(Cheb3D);
    _verify(EccAnom);
    _verify(Frac);
    _verify(Geodetic);
    _verify(Legendre);
    _verify(MeanObliquity);
    _verify(Mjday);
    _verify(Mjday_TDB);
    _verify(NutAngles);
    _verify(Position);
    _verify(IERS);
    _verify(PoleMatrix);
    _verify(PrecMatrix);
    _verify(NutMatrix);
    _verify(gmst);
    _verify(EqnEquinox);
    _verify(gast);
    _verify(GHAMatrix);
    _verify(AccelHarmonic);
    _verify(LTC);
    _verify(elements);
    _verify(angl);
    _verify(TimeUpdate);
    _verify(MeasUpdate);
    _verify(JPL_Eph_DE430);





    return 0;
}


int main()
{

    //Global::eop19620101(6);

    //Global::eopdata->print();



    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}

