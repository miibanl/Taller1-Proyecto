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

