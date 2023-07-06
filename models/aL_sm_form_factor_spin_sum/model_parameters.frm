#define PSI "3370,3371,3372,3373,3374,3375,3376,3377,3378,3379"
#define GLU "21"
#define PHO "22"
#define EP "-11"
#define EM "11"
#define H "25"
#define GHO "82"
#define GHOBAR "-82"
#define FERM "-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,-11,11,-12,12,-13,13"
#define Q "1,2,3,4,5,6"
#define QBAR "-1,-2,-3,-4,-5,-6"
#define QMASSLESS "1,2,3,4,5"
#define QBARMASSLESS "-1,-2,-3,-4,-5"
#define QMASSIVE "6,"
#define QBARMASSIVE "-6,"
#define L "11,12,13"
#define LBAR "-11,-12,-13"
#define Z "23"
#define FORMFACTORS "1"
#define HASAMP "0"

**************************************************
* START SE PDGs
**************************************************
#define PHOPRIME "1022"
#define QMASSIVEPRIME "1006,"
#define QBARMASSIVEPRIME "-1006,"
#define SDUMMY "1122"
#define GLUPRIME "1021"
#define GHOPRIME "1066"
#define GHOPRIMEBAR "-1066"
**************************************************
* END SE PDGs
**************************************************

**************************************************
* START amp PDGs
**************************************************
#define PHOAMPPRIME "2022"
CF  FFS, FFT, FFU;
CF APHOAMPFFSTU;
CF APHOAMPFFTSU;
CF APHOAMPFFUST;
CF BPHOAMPFFSTU;
CF BPHOAMPFFTSU;
CF BPHOAMPFFUST;
CF CPHOAMPFFSTU;
**************************************************
* END amp PDGs
**************************************************

S vev, pi, cw ,sw2 , gw;

Auto S mass;
Auto S yukawa;
CTable masses(-10000:10000);
CTable gyq(-10000:10000);
CTable logmasses(-10000:10000);
CTable charges(-10000:10000);
CTable zVcoupling(-10000:10000);
CTable zAcoupling(-10000:10000);

#ifndef `OPTIMLVL'
    #define OPTIMLVL "4"
#endif

#ifndef `OPTIMITERATIONS'
    #define OPTIMITERATIONS "100"
#endif

Fill gyq(1) = yukawad; * d
Fill gyq(2) = yukawau; * u
Fill gyq(3) = yukawas; * s
Fill gyq(4) = yukawac; * c
Fill gyq(5) = yukawab; * b
Fill gyq(6) = yukawat; * t
Fill gyq(11) = 0; * e-
Fill gyq(12) = 0; * mu-
Fill gyq(13) = 0; * ta-
Fill gyq(-1) = yukawad; * d
Fill gyq(-2) = yukawau; * u
Fill gyq(-3) = yukawas; * s
Fill gyq(-4) = yukawac; * c
Fill gyq(-5) = yukawab; * b
Fill gyq(-6) = yukawat; * t
Fill gyq(-11) = 0; * e+
Fill gyq(-12) = 0; * mu+
Fill gyq(-13) = 0; * ta+

#ifndef `HEAVYFERMIONS'
Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = massb;
Fill masses(6) = masst;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = massb;
Fill masses(-6) = masst;
Fill masses(11) = 0;
Fill masses(12) = 0;
Fill masses(13) = 0;
Fill masses(-11) = 0;
Fill masses(-12) = 0;
Fill masses(-13) = 0;
#else
Fill masses(1) = massd;
Fill masses(2) = massu;
Fill masses(3) = massc;
Fill masses(4) = masss;
Fill masses(5) = massb;
Fill masses(6) = masst;
Fill masses(-1) = massd;
Fill masses(-2) = massu;
Fill masses(-3) = massc;
Fill masses(-4) = masss;
Fill masses(-5) = massb;
Fill masses(-6) = masst;
Fill masses(11) = masse;
Fill masses(12) = massmu;
Fill masses(13) = masstau;
Fill masses(-11) = masse;
Fill masses(-12) = massmu;
Fill masses(-13) = masstau;
#endif

Fill masses(-82) = 0;
Fill masses(82) = 0;
Fill masses(21) = 0;
Fill masses(22) = 0;
Fill masses(25) = massh;
Fill masses(3370) = 0;
Fill masses(3371) = massdummya;
Fill masses(3372) = massdummyb;
Fill masses(3373) = massdummyc;
Fill masses(3374) = massdummyd;
Fill masses(3375) = massdummye;
Fill masses(3376) = massdummyf;
Fill masses(3377) = massdummyg;
Fill masses(3378) = massdummyh;
Fill masses(3379) = massdummyi;

**************************************************
* START EW parameters
**************************************************

Fill masses(23) = massz;

Fill zVcoupling(1) = -1/2+2/3*sw2; * d
Fill zVcoupling(2) = 1/2-4/3*sw2; * u
Fill zVcoupling(3) = -1/2+2/3*sw2; * s
Fill zVcoupling(4) = 1/2-4/3*sw2; * c
Fill zVcoupling(5) = -1/2+2/3*sw2; * b
Fill zVcoupling(6) = 1/2-4/3*sw2; * t
Fill zVcoupling(11) = -1/2+2*sw2; * e-
Fill zVcoupling(12) = -1/2+2*sw2; * mu-
Fill zVcoupling(13) = -1/2+2*sw2; * ta-

Fill zAcoupling(1) = -1/2; * d
Fill zAcoupling(2) = 1/2; * u
Fill zAcoupling(3) = -1/2; * s
Fill zAcoupling(4) = 1/2; * c
Fill zAcoupling(5) = -1/2; * b
Fill zAcoupling(6) = 1/2; * t
Fill zAcoupling(11) = -1/2; * e-
Fill zAcoupling(12) = -1/2; * mu-
Fill zAcoupling(13) = -1/2; * ta-

**************************************************
* END EW parameters
**************************************************

**************************************************
* START SE parameters
**************************************************
Fill masses(-1006) = masst;
Fill masses(1006) = masst;
Fill charges(-1006) = -2/3;
Fill charges(1006) = 2/3;
Fill gyq(-1006) = yukawat;
Fill gyq(1006) = yukawat;
Fill masses(1122) = 0;
Fill masses(1021) = 0;
Fill masses(1066) = 0;
Fill masses(-1066) = 0;
**************************************************
* END SE parameters
**************************************************

**************************************************
* START amp parameters
**************************************************
Fill masses(2022) = 0;
**************************************************
* END amp parameters
**************************************************

* note: this set needs to be complete for the UV expansion
Set allmasses: massu, massd, massc, masss, masst, massb, masse, massmu, masstau, massh, massw, massz;

Fill charges(1) = -1/3; * d
Fill charges(2) = 2/3; * u
Fill charges(3) = -1/3; * s
Fill charges(4) = 2/3; * c
Fill charges(5) = -1/3; * b
Fill charges(6) = 2/3; * t
Fill charges(11) = -1; * e
Fill charges(-1) = 1/3; * d
Fill charges(-2) = -2/3; * u
Fill charges(-3) = 1/3; * s
Fill charges(-4) = -2/3; * c
Fill charges(-5) = 1/3; * b
Fill charges(-6) = -2/3; * t
Fill charges(-11) = 1; * e


Set ffvecs:p1, p2, k1, k2, k3 ,k4,-p1, -p2, -k1, -k2, -k3, -k4;    