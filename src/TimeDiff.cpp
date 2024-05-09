//
// Created by miibanl on 24/04/2024.
//

#include "../include/TimeDiff.h"


Matrix timeDiff(double UT1_UTC, double TAI_UTC){
    double TT_TAI  = +32.184;          // TT-TAI time difference [s]

    double GPS_TAI = -19.0;            // GPS-TAI time difference [s]

    double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

    double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    double UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]

    double UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

    double UT1_GPS = UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]

    double TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

    double GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]


    Matrix resultado(1,5);

    resultado(1,1)=UT1_TAI;
    resultado(1,2)=UTC_GPS;
    resultado(1,3)=UT1_GPS;
    resultado(1,4)=TT_UTC;
    resultado(1,5)=GPS_UTC;


    return resultado;
}