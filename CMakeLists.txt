cmake_minimum_required(VERSION 3.23.2)
project(Proyecto)

set(CMAKE_CXX_STANDARD 17)

add_executable(Proyecto EKF_Test.cpp include/AccelPointMass.h include/AzElPa.h include/Cheb3D.h include/EccAnom.h include/Frac.h include/Geodetic.h include/Global.h include/IERS.h include/Legendre.h include/Matrix.h include/MeanObliquity.h include/Mjday.h include/Mjday_TDB.h include/NutAngles.h include/Position.h include/R_x.h include/R_y.h include/R_z.h include/SAT_Const.h include/sign.h include/TimeDiff.h include/Unit.h src/AccelPointMass.cpp src/AzElPa.cpp src/Cheb3D.cpp src/EccAnom.cpp src/Frac.cpp src/Geodetic.cpp src/Global.cpp src/IERS.cpp src/Legendre.cpp src/Matrix.cpp src/MeanObliquity.cpp src/Mjday.cpp src/Mjday_TDB.cpp src/NutAngles.cpp src/Position.cpp src/R_x.cpp src/R_y.cpp src/R_z.cpp src/sign.cpp src/TimeDiff.cpp src/Unit.cpp src/AccelHarmonic.cpp include/AccelHarmonic.h EKF_GEOS3.cpp)
