#include <cmath>
#include <complex>
#include <math.h>
#include <mp++/complex.hpp>

using namespace std;

#include "cavh_olo.h"

// this is a constant in here for now, later this should be an argument of the
// functions
complex<double> M_f64(173.0, 0.0);

complex<double> top_Q8_f64(16.0 / 81.0, 0);

double FORM_FACTOR_PHASE = 1.0;

complex<double> zero_f64(0.0, 0.0);
double scale = 10.0;
int zero = 0;

void APHOAMPFFSTU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64);
  complex<double> b02s = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> astu;
  complex<double> Z[18];

  Z[0] = 3.0 * c02u;
  Z[1] = Z[0] + 3.0 * b02u;
  Z[2] = 16.0 * d02tu;
  Z[3] = 3.0 * b02s;
  Z[4] = 24.0 * c02s - Z[2] - 16.0 * d02su + Z[3] + 2.0 - Z[1];
  Z[5] = 4.0 * c02s;
  Z[6] = 9.0 * c02u;
  Z[7] = -Z[5] - Z[6];
  Z[8] = 3.0 * xu;
  Z[7] = Z[7] * Z[8];
  Z[4] = 4.0 * Z[4] + Z[7];
  Z[7] = xu * xu;
  Z[4] = Z[4] * Z[7];
  Z[9] = 3.0 * c02t;
  Z[10] = 4.0 * d02st;
  Z[11] = -Z[9] - Z[10];
  Z[12] = c02t - c02s;
  Z[13] = xu * d02tu;
  Z[12] = 3.0 * Z[12] - 5.0 * Z[13];
  Z[12] = xu * Z[12];
  Z[11] = 8.0 * Z[11] + Z[12];
  Z[11] = xt * Z[11];
  Z[12] = 9.0 * c02s;
  Z[2] = -Z[12] + Z[2] - 9.0 * c02t;
  Z[2] = xu * Z[2];
  Z[2] = 4.0 * b02t + Z[2];
  Z[2] = xu * Z[2];
  Z[2] = Z[2] + Z[11];
  Z[11] = 2.0 * xt;
  Z[2] = Z[2] * Z[11];
  Z[2] = Z[4] + Z[2];
  Z[2] = xt * Z[2];
  Z[4] = Z[3] - b02u;
  Z[6] = Z[6] + Z[4];
  Z[13] = xu * xu * xu;
  Z[6] = Z[6] * Z[13];
  Z[2] = 4.0 * Z[6] + Z[2];
  Z[2] = xt * Z[2];
  Z[6] = 3.0 * c02s;
  Z[14] = -Z[10] - Z[6];
  Z[15] = -2.0 * d02st - c02s;
  Z[16] = xu * d02su;
  Z[15] = 12.0 * Z[15] - 11.0 * Z[16];
  Z[15] = xu * Z[15];
  Z[16] = -xt + Z[8];
  Z[16] = d02st * Z[16];
  Z[17] = Z[10] + c02s;
  Z[16] = 2.0 * Z[17] + Z[16];
  Z[16] = Z[16] * Z[11];
  Z[14] = Z[16] + 8.0 * Z[14] + Z[15];
  Z[14] = xt * Z[14];
  Z[15] = d02su - c02s;
  Z[15] = 8.0 * Z[15] + Z[0];
  Z[16] = 2.0 * Z[7];
  Z[15] = Z[15] * Z[16];
  Z[14] = Z[15] + Z[14];
  Z[14] = xt * Z[14];
  Z[15] = 2.0 * xs;
  Z[17] = -Z[15] - 10.0 * xt;
  Z[17] = d02su * Z[17];
  Z[5] = Z[17] - Z[5] + Z[0];
  Z[5] = xs * Z[7] * Z[5];
  Z[17] = Z[13] * c02u;
  Z[5] = Z[5] - Z[17] + Z[14];
  Z[5] = Z[5] * Z[15];
  Z[9] = Z[9] + 8.0 * d02st;
  Z[6] = -Z[6] - Z[9];
  Z[14] = d02st - d02su - d02tu;
  Z[14] = xu * Z[14];
  Z[9] = Z[14] - c02s - Z[9];
  Z[8] = Z[9] * Z[8];
  Z[9] = c02t + Z[10];
  Z[9] = xt * Z[9];
  Z[6] = 4.0 * Z[9] + 8.0 * Z[6] + Z[8];
  Z[6] = Z[6] * Z[11];
  Z[8] = 5.0 * d02su + 3.0 * d02tu;
  Z[8] = 4.0 * Z[8] - Z[12];
  Z[8] = 4.0 * Z[8] + Z[0];
  Z[8] = xu * Z[8];
  Z[8] = 8.0 * b02s + Z[8];
  Z[8] = xu * Z[8];
  Z[6] = Z[8] + Z[6];
  Z[6] = xt * Z[6];
  Z[8] = xu * c02u;
  Z[1] = 12.0 * c02s + 5.0 * b02s - Z[1];
  Z[1] = 2.0 * Z[1] - 13.0 * Z[8];
  Z[1] = Z[1] * Z[16];
  Z[1] = Z[1] + Z[6];
  Z[1] = xt * Z[1];
  Z[1] = Z[1] + Z[5];
  Z[1] = xs * Z[1];
  Z[1] = Z[2] + Z[1];
  Z[1] = xs * Z[1];
  Z[2] = 6.0 * c02t + b02t - Z[3];
  Z[2] = 8.0 * Z[2] - 5.0 * Z[8];
  Z[2] = Z[2] * Z[7];
  Z[3] = xt * d02tu;
  Z[5] = Z[3] - 4.0 * d02tu + c02t;
  Z[0] = Z[0] - 8.0 * Z[5];
  Z[5] = Z[7] * xt;
  Z[0] = Z[0] * Z[5];
  Z[0] = Z[2] + Z[0];
  Z[0] = xt * Z[0];
  Z[2] = 6.0 * c02u - Z[4];
  Z[2] = Z[2] * Z[13];
  Z[0] = 8.0 * Z[2] + Z[0];
  Z[0] = Z[0] * xt * xt;
  Z[0] = Z[0] + Z[1];
  Z[0] = xs * Z[0];
  Z[1] = -2.0 * c02t - Z[3];
  Z[1] = Z[1] * Z[5];
  Z[1] = -2.0 * Z[17] + Z[1];
  Z[1] = Z[1] * xt * xt * xt * xt;
  Z[0] = 4.0 * Z[1] + Z[0];

  astu = M_f64 * M_f64 * M_f64 * M_f64 * Z[0];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * astu /
         (96.0 * xs * xs * xt * xt * xu * xu * E1 * E2 * E3 * E42 * s * s);
}

void APHOAMPFFTSU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64);
  complex<double> b02s = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> atsu;
  complex<double> Z[18];

  Z[0] = 3.0 * c02u;
  Z[1] = Z[0] + 3.0 * b02u;
  Z[2] = 16.0 * d02su;
  Z[3] = 3.0 * b02t;
  Z[4] = 24.0 * c02t - Z[2] - 16.0 * d02tu + Z[3] + 2.0 - Z[1];
  Z[5] = 4.0 * c02t;
  Z[6] = 9.0 * c02u;
  Z[7] = -Z[5] - Z[6];
  Z[8] = 3.0 * xu;
  Z[7] = Z[7] * Z[8];
  Z[4] = 4.0 * Z[4] + Z[7];
  Z[7] = xu * xu;
  Z[4] = Z[4] * Z[7];
  Z[9] = 3.0 * c02s;
  Z[10] = 4.0 * d02st;
  Z[11] = -Z[9] - Z[10];
  Z[12] = c02s - c02t;
  Z[13] = xu * d02su;
  Z[12] = 3.0 * Z[12] - 5.0 * Z[13];
  Z[12] = xu * Z[12];
  Z[11] = 8.0 * Z[11] + Z[12];
  Z[11] = xs * Z[11];
  Z[12] = 9.0 * c02t;
  Z[2] = -Z[12] + Z[2] - 9.0 * c02s;
  Z[2] = xu * Z[2];
  Z[2] = 4.0 * b02s + Z[2];
  Z[2] = xu * Z[2];
  Z[2] = Z[2] + Z[11];
  Z[11] = 2.0 * xs;
  Z[2] = Z[2] * Z[11];
  Z[2] = Z[4] + Z[2];
  Z[2] = xs * Z[2];
  Z[4] = Z[3] - b02u;
  Z[6] = Z[6] + Z[4];
  Z[13] = xu * xu * xu;
  Z[6] = Z[6] * Z[13];
  Z[2] = 4.0 * Z[6] + Z[2];
  Z[2] = xs * Z[2];
  Z[6] = 3.0 * c02t;
  Z[14] = -Z[10] - Z[6];
  Z[15] = -2.0 * d02st - c02t;
  Z[16] = xu * d02tu;
  Z[15] = 12.0 * Z[15] - 11.0 * Z[16];
  Z[15] = xu * Z[15];
  Z[16] = -xs + Z[8];
  Z[16] = d02st * Z[16];
  Z[17] = Z[10] + c02t;
  Z[16] = 2.0 * Z[17] + Z[16];
  Z[16] = Z[16] * Z[11];
  Z[14] = Z[16] + 8.0 * Z[14] + Z[15];
  Z[14] = xs * Z[14];
  Z[15] = d02tu - c02t;
  Z[15] = 8.0 * Z[15] + Z[0];
  Z[16] = 2.0 * Z[7];
  Z[15] = Z[15] * Z[16];
  Z[14] = Z[15] + Z[14];
  Z[14] = xs * Z[14];
  Z[15] = 2.0 * xt;
  Z[17] = -Z[15] - 10.0 * xs;
  Z[17] = d02tu * Z[17];
  Z[5] = Z[17] - Z[5] + Z[0];
  Z[5] = xt * Z[7] * Z[5];
  Z[17] = Z[13] * c02u;
  Z[5] = Z[5] - Z[17] + Z[14];
  Z[5] = Z[5] * Z[15];
  Z[9] = Z[9] + 8.0 * d02st;
  Z[6] = -Z[6] - Z[9];
  Z[14] = d02st - d02tu - d02su;
  Z[14] = xu * Z[14];
  Z[9] = Z[14] - c02t - Z[9];
  Z[8] = Z[9] * Z[8];
  Z[9] = c02s + Z[10];
  Z[9] = xs * Z[9];
  Z[6] = 4.0 * Z[9] + 8.0 * Z[6] + Z[8];
  Z[6] = Z[6] * Z[11];
  Z[8] = 5.0 * d02tu + 3.0 * d02su;
  Z[8] = 4.0 * Z[8] - Z[12];
  Z[8] = 4.0 * Z[8] + Z[0];
  Z[8] = xu * Z[8];
  Z[8] = 8.0 * b02t + Z[8];
  Z[8] = xu * Z[8];
  Z[6] = Z[8] + Z[6];
  Z[6] = xs * Z[6];
  Z[8] = xu * c02u;
  Z[1] = 12.0 * c02t + 5.0 * b02t - Z[1];
  Z[1] = 2.0 * Z[1] - 13.0 * Z[8];
  Z[1] = Z[1] * Z[16];
  Z[1] = Z[1] + Z[6];
  Z[1] = xs * Z[1];
  Z[1] = Z[1] + Z[5];
  Z[1] = xt * Z[1];
  Z[1] = Z[2] + Z[1];
  Z[1] = xt * Z[1];
  Z[2] = 6.0 * c02s + b02s - Z[3];
  Z[2] = 8.0 * Z[2] - 5.0 * Z[8];
  Z[2] = Z[2] * Z[7];
  Z[3] = xs * d02su;
  Z[5] = Z[3] - 4.0 * d02su + c02s;
  Z[0] = Z[0] - 8.0 * Z[5];
  Z[5] = Z[7] * xs;
  Z[0] = Z[0] * Z[5];
  Z[0] = Z[2] + Z[0];
  Z[0] = xs * Z[0];
  Z[2] = 6.0 * c02u - Z[4];
  Z[2] = Z[2] * Z[13];
  Z[0] = 8.0 * Z[2] + Z[0];
  Z[0] = Z[0] * xs * xs;
  Z[0] = Z[0] + Z[1];
  Z[0] = xt * Z[0];
  Z[1] = -2.0 * c02s - Z[3];
  Z[1] = Z[1] * Z[5];
  Z[1] = -2.0 * Z[17] + Z[1];
  Z[1] = Z[1] * xs * xs * xs * xs;
  Z[0] = 4.0 * Z[1] + Z[0];

  atsu = M_f64 * M_f64 * M_f64 * M_f64 * Z[0];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * atsu /
         (96.0 * xs * xs * xt * xt * xu * xu * E1 * E2 * E3 * E42 * t * t);
}

void APHOAMPFFUST_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64);
  complex<double> b02s = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> aust;
  complex<double> Z[17];

  Z[0] = 3.0 * c02t;
  Z[1] = Z[0] + 3.0 * b02t;
  Z[2] = 16.0 * d02st;
  Z[3] = 3.0 * b02u;
  Z[4] = 24.0 * c02u - Z[2] - 16.0 * d02tu + Z[3] + 2.0 - Z[1];
  Z[5] = 9.0 * c02t;
  Z[6] = 4.0 * c02u;
  Z[7] = -Z[5] - Z[6];
  Z[8] = 3.0 * xt;
  Z[7] = Z[7] * Z[8];
  Z[4] = 4.0 * Z[4] + Z[7];
  Z[7] = xt * xt;
  Z[4] = Z[4] * Z[7];
  Z[9] = 3.0 * c02s;
  Z[10] = 4.0 * d02su;
  Z[11] = -Z[9] - Z[10];
  Z[12] = c02s - c02u;
  Z[13] = xt * d02st;
  Z[12] = 3.0 * Z[12] - 5.0 * Z[13];
  Z[12] = xt * Z[12];
  Z[11] = 8.0 * Z[11] + Z[12];
  Z[11] = xs * Z[11];
  Z[12] = c02u + c02s;
  Z[2] = Z[2] - 9.0 * Z[12];
  Z[2] = xt * Z[2];
  Z[2] = 4.0 * b02s + Z[2];
  Z[2] = xt * Z[2];
  Z[2] = Z[2] + Z[11];
  Z[11] = 2.0 * xs;
  Z[2] = Z[2] * Z[11];
  Z[2] = Z[4] + Z[2];
  Z[2] = xs * Z[2];
  Z[4] = Z[3] - b02t;
  Z[5] = Z[5] + Z[4];
  Z[12] = xt * xt * xt;
  Z[5] = Z[5] * Z[12];
  Z[2] = 4.0 * Z[5] + Z[2];
  Z[2] = xs * Z[2];
  Z[5] = 3.0 * c02u;
  Z[13] = -Z[10] - Z[5];
  Z[14] = -2.0 * d02su - c02u;
  Z[15] = xt * d02tu;
  Z[14] = 12.0 * Z[14] - 11.0 * Z[15];
  Z[14] = xt * Z[14];
  Z[15] = -xs + Z[8];
  Z[15] = d02su * Z[15];
  Z[16] = Z[10] + c02u;
  Z[15] = 2.0 * Z[16] + Z[15];
  Z[15] = Z[15] * Z[11];
  Z[13] = Z[15] + 8.0 * Z[13] + Z[14];
  Z[13] = xs * Z[13];
  Z[14] = c02u - d02tu;
  Z[14] = Z[0] - 8.0 * Z[14];
  Z[15] = 2.0 * Z[7];
  Z[14] = Z[14] * Z[15];
  Z[13] = Z[14] + Z[13];
  Z[13] = xs * Z[13];
  Z[14] = 2.0 * xu;
  Z[16] = -Z[14] - 10.0 * xs;
  Z[16] = d02tu * Z[16];
  Z[6] = Z[16] + Z[0] - Z[6];
  Z[6] = xu * Z[7] * Z[6];
  Z[16] = Z[12] * c02t;
  Z[6] = Z[6] - Z[16] + Z[13];
  Z[6] = Z[6] * Z[14];
  Z[9] = Z[9] + 8.0 * d02su;
  Z[5] = -Z[5] - Z[9];
  Z[13] = d02su - d02tu - d02st;
  Z[13] = xt * Z[13];
  Z[9] = Z[13] - c02u - Z[9];
  Z[8] = Z[9] * Z[8];
  Z[9] = c02s + Z[10];
  Z[9] = xs * Z[9];
  Z[5] = 4.0 * Z[9] + 8.0 * Z[5] + Z[8];
  Z[5] = Z[5] * Z[11];
  Z[8] = 5.0 * d02tu + 3.0 * d02st;
  Z[8] = -36.0 * c02u + 16.0 * Z[8] + Z[0];
  Z[8] = xt * Z[8];
  Z[8] = 8.0 * b02u + Z[8];
  Z[8] = xt * Z[8];
  Z[5] = Z[8] + Z[5];
  Z[5] = xs * Z[5];
  Z[8] = xt * c02t;
  Z[1] = 12.0 * c02u + 5.0 * b02u - Z[1];
  Z[1] = 2.0 * Z[1] - 13.0 * Z[8];
  Z[1] = Z[1] * Z[15];
  Z[1] = Z[1] + Z[5];
  Z[1] = xs * Z[1];
  Z[1] = Z[1] + Z[6];
  Z[1] = xu * Z[1];
  Z[1] = Z[2] + Z[1];
  Z[1] = xu * Z[1];
  Z[2] = 6.0 * c02s + b02s - Z[3];
  Z[2] = 8.0 * Z[2] - 5.0 * Z[8];
  Z[2] = Z[2] * Z[7];
  Z[3] = xs * d02st;
  Z[5] = Z[3] - 4.0 * d02st + c02s;
  Z[0] = Z[0] - 8.0 * Z[5];
  Z[5] = Z[7] * xs;
  Z[0] = Z[0] * Z[5];
  Z[0] = Z[2] + Z[0];
  Z[0] = xs * Z[0];
  Z[2] = 6.0 * c02t - Z[4];
  Z[2] = Z[2] * Z[12];
  Z[0] = 8.0 * Z[2] + Z[0];
  Z[0] = Z[0] * xs * xs;
  Z[0] = Z[0] + Z[1];
  Z[0] = xu * Z[0];
  Z[1] = -2.0 * c02s - Z[3];
  Z[1] = Z[1] * Z[5];
  Z[1] = -2.0 * Z[16] + Z[1];
  Z[1] = Z[1] * xs * xs * xs * xs;
  Z[0] = 4.0 * Z[1] + Z[0];

  aust = M_f64 * M_f64 * M_f64 * M_f64 * Z[0];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * aust /
         (96.0 * xs * xs * xt * xt * xu * xu * E1 * E2 * E3 * E42 * u * u);
}

void BPHOAMPFFSTU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64);
  complex<double> b02s = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> bstu;
  complex<double> Z[29];

  Z[0] = 2.0 * d02st;
  Z[1] = 3.0 * c02t;
  Z[2] = Z[0] + Z[1];
  Z[3] = 3.0 * xu;
  Z[4] = c02s + d02st;
  Z[5] = c02t - Z[4];
  Z[5] = Z[5] * Z[3];
  Z[2] = 4.0 * Z[2] + Z[5];
  Z[5] = d02st * Z[3];
  Z[6] = -xt * Z[0];
  Z[4] = Z[6] + 4.0 * Z[4] + Z[5];
  Z[4] = xt * Z[4];
  Z[2] = 2.0 * Z[2] + Z[4];
  Z[4] = 2.0 * xt;
  Z[2] = Z[2] * Z[4];
  Z[5] = xu * d02su;
  Z[6] = 37.0 * Z[5];
  Z[7] = 26.0 * d02su - 33.0 * c02s;
  Z[7] = 2.0 * Z[7] + Z[6];
  Z[8] = xu * xu;
  Z[7] = Z[7] * Z[8];
  Z[2] = Z[7] + Z[2];
  Z[2] = xt * Z[2];
  Z[7] = b02s + b02u;
  Z[9] = 3.0 * c02s;
  Z[10] = Z[9] + Z[7];
  Z[11] = 2.0 * d02su;
  Z[12] = Z[11] - Z[9];
  Z[13] = 7.0 * c02u;
  Z[12] = 4.0 * Z[12] - Z[13];
  Z[12] = 5.0 * Z[12] + 7.0 * Z[5];
  Z[12] = xu * Z[12];
  Z[10] = 8.0 * Z[10] + Z[12];
  Z[12] = 2.0 * Z[8];
  Z[10] = Z[10] * Z[12];
  Z[2] = Z[10] + Z[2];
  Z[2] = xt * Z[2];
  Z[10] = 3.0 * c02u;
  Z[14] = -Z[5] - 4.0 * c02s - Z[10];
  Z[15] = xu * xu * xu;
  Z[16] = 4.0 * Z[15];
  Z[14] = Z[14] * Z[16];
  Z[17] = 10.0 * d02su;
  Z[18] = -4.0 * Z[5] + Z[17] - 19.0 * c02s;
  Z[12] = Z[18] * Z[12];
  Z[18] = Z[8] * xt;
  Z[19] = Z[18] * d02su;
  Z[12] = Z[12] - 7.0 * Z[19];
  Z[12] = xt * Z[12];
  Z[20] = -c02s - Z[5];
  Z[21] = xs * d02su;
  Z[20] = -4.0 * Z[21] + 8.0 * Z[20];
  Z[20] = Z[8] * Z[20];
  Z[19] = -11.0 * Z[19] + Z[20];
  Z[19] = xs * Z[19];
  Z[12] = Z[19] + Z[14] + Z[12];
  Z[12] = xs * Z[12];
  Z[14] = 2.0 * c02u;
  Z[19] = Z[14] + c02s;
  Z[20] = xu * xu * xu * xu;
  Z[21] = Z[19] * Z[20];
  Z[2] = 2.0 * Z[12] - 16.0 * Z[21] + Z[2];
  Z[2] = xs * Z[2];
  Z[12] = 9.0 * c02t;
  Z[21] = -Z[9] + 4.0 * d02st + Z[12];
  Z[22] = xu + 1.0;
  Z[22] = d02st * Z[22];
  Z[22] = c02t + Z[22];
  Z[22] = Z[22] * Z[3];
  Z[23] = Z[9] + 7.0 * d02st + c02t;
  Z[24] = xt * d02st;
  Z[23] = 2.0 * Z[23] - 3.0 * Z[24];
  Z[23] = xt * Z[23];
  Z[21] = Z[23] + 4.0 * Z[21] + Z[22];
  Z[21] = xt * Z[21];
  Z[22] = Z[1] - Z[9];
  Z[23] = b02t - b02s;
  Z[24] = Z[22] - Z[23];
  Z[25] = -Z[0] + c02t;
  Z[25] = 2.0 * Z[25] - 5.0 * c02s;
  Z[25] = 6.0 * Z[25] + Z[6];
  Z[25] = xu * Z[25];
  Z[24] = 16.0 * Z[24] + Z[25];
  Z[24] = xu * Z[24];
  Z[21] = Z[24] + 4.0 * Z[21];
  Z[21] = xt * Z[21];
  Z[24] = 5.0 * b02s;
  Z[25] = Z[24] + 3.0 * b02t;
  Z[26] = 8.0 * d02su;
  Z[27] = 24.0 * c02s - Z[26] + 9.0 * b02u + Z[25];
  Z[28] = 21.0 * c02u;
  Z[6] = Z[6] - 46.0 * c02s - Z[28];
  Z[6] = xu * Z[6];
  Z[6] = 4.0 * Z[27] + Z[6];
  Z[6] = xu * Z[6];
  Z[6] = -32.0 * b02t + Z[6];
  Z[6] = xu * Z[6];
  Z[6] = Z[6] + Z[21];
  Z[6] = xt * Z[6];
  Z[21] = 9.0 * c02u;
  Z[27] = 6.0 * c02s + Z[7];
  Z[27] = 2.0 * Z[27] + Z[21];
  Z[17] = -Z[28] + Z[17] - 7.0 * c02s;
  Z[17] = xu * Z[17];
  Z[17] = 2.0 * Z[27] + Z[17];
  Z[16] = Z[17] * Z[16];
  Z[6] = Z[16] + Z[6];
  Z[6] = xt * Z[6];
  Z[16] = xu * xu * xu * xu * xu;
  Z[17] = c02u * Z[16];
  Z[2] = Z[2] - 8.0 * Z[17] + Z[6];
  Z[2] = xs * Z[2];
  Z[6] = -2.0 * Z[23] + Z[22];
  Z[17] = Z[0] - d02tu;
  Z[11] = Z[11] + Z[17];
  Z[11] = xu * Z[11];
  Z[17] = Z[17] + c02s;
  Z[17] = -c02t - 2.0 * Z[17];
  Z[11] = 2.0 * Z[17] + Z[11];
  Z[11] = Z[11] * Z[3];
  Z[6] = 16.0 * Z[6] + Z[11];
  Z[6] = xu * Z[6];
  Z[11] = d02st + Z[1];
  Z[11] = 2.0 * Z[11] - Z[9];
  Z[0] = Z[0] + c02t;
  Z[17] = d02st - d02tu;
  Z[17] = xu * Z[17];
  Z[17] = 2.0 * Z[0] + Z[17];
  Z[17] = Z[17] * Z[3];
  Z[0] = xt * Z[0];
  Z[0] = 12.0 * Z[0] + 8.0 * Z[11] + Z[17];
  Z[0] = Z[0] * Z[4];
  Z[0] = Z[6] + Z[0];
  Z[0] = xt * Z[0];
  Z[4] = -Z[9] + 3.0 * d02tu - 23.0 * d02su;
  Z[4] = 3.0 * Z[5] + 2.0 * Z[4] + 5.0 * c02u;
  Z[4] = xu * Z[4];
  Z[5] = Z[24] + 1.0 + 6.0 * b02u;
  Z[5] = 12.0 * c02s + 2.0 * Z[5] - 5.0 * b02t;
  Z[4] = 2.0 * Z[5] + Z[4];
  Z[4] = xu * Z[4];
  Z[5] = -b02s - b02t;
  Z[4] = 16.0 * Z[5] + Z[4];
  Z[4] = xu * Z[4];
  Z[0] = 2.0 * Z[4] + Z[0];
  Z[0] = xt * Z[0];
  Z[4] = 2.0 * b02t;
  Z[5] = -Z[4] - Z[7];
  Z[6] = 18.0 * c02u + 36.0 * c02s - 16.0 * d02su + b02u + Z[25];
  Z[9] = -92.0 * d02su - c02u;
  Z[9] = xu * Z[9];
  Z[6] = 4.0 * Z[6] + Z[9];
  Z[6] = xu * Z[6];
  Z[5] = 32.0 * Z[5] + Z[6];
  Z[5] = Z[5] * Z[8];
  Z[0] = Z[5] + Z[0];
  Z[0] = xt * Z[0];
  Z[5] = xu * c02u;
  Z[6] = 24.0 * Z[19] - 7.0 * Z[5];
  Z[6] = Z[6] * Z[20];
  Z[0] = 2.0 * Z[6] + Z[0];
  Z[0] = xt * Z[0];
  Z[0] = Z[0] + Z[2];
  Z[0] = xs * Z[0];
  Z[2] = 8.0 * d02tu;
  Z[6] = 24.0 * c02u;
  Z[9] = 3.0 * b02u + b02s;
  Z[9] = Z[6] - Z[26] - Z[2] + 3.0 * Z[9] - 4.0 * b02t;
  Z[9] = 4.0 * Z[9] - 27.0 * Z[5];
  Z[9] = xu * Z[9];
  Z[11] = -b02t - b02u - 2.0 * b02s;
  Z[9] = 32.0 * Z[11] + Z[9];
  Z[9] = Z[9] * Z[8];
  Z[2] = 12.0 * c02t - Z[2] + Z[4] + 7.0 * b02u + 3.0 * b02s;
  Z[4] = 4.0 * d02tu;
  Z[11] = -Z[4] - Z[12];
  Z[3] = -d02tu * Z[3];
  Z[3] = Z[3] + 2.0 * Z[11] - 39.0 * c02u;
  Z[3] = xu * Z[3];
  Z[2] = 4.0 * Z[2] + Z[3];
  Z[2] = xu * Z[2];
  Z[2] = -32.0 * b02s + Z[2];
  Z[2] = xu * Z[2];
  Z[3] = Z[4] - Z[1];
  Z[11] = xu * d02tu;
  Z[3] = 2.0 * Z[3] - Z[11];
  Z[3] = Z[3] * Z[8];
  Z[12] = Z[18] * d02tu;
  Z[3] = 5.0 * Z[3] - 14.0 * Z[12];
  Z[3] = xt * Z[3];
  Z[2] = Z[2] + Z[3];
  Z[2] = xt * Z[2];
  Z[2] = Z[9] + Z[2];
  Z[2] = xt * Z[2];
  Z[3] = -b02t - Z[7];
  Z[7] = -4.0 * d02su + Z[21];
  Z[7] = xu * Z[7];
  Z[3] = 4.0 * Z[3] + Z[7];
  Z[3] = Z[3] * Z[15];
  Z[2] = 8.0 * Z[3] + Z[2];
  Z[2] = xt * Z[2];
  Z[3] = Z[16] * Z[6];
  Z[2] = Z[3] + Z[2];
  Z[2] = xt * Z[2];
  Z[0] = Z[2] + Z[0];
  Z[0] = xs * Z[0];
  Z[2] = Z[1] + b02t;
  Z[3] = Z[10] + Z[2];
  Z[1] = -Z[1] - Z[13];
  Z[1] = xu * Z[1];
  Z[1] = 4.0 * Z[3] + Z[1];
  Z[1] = Z[1] * Z[15];
  Z[3] = d02tu - c02t;
  Z[3] = 10.0 * Z[3] - Z[11];
  Z[3] = Z[3] * Z[8];
  Z[4] = -Z[18] * Z[4];
  Z[3] = Z[3] + Z[4];
  Z[3] = xt * Z[3];
  Z[4] = c02u - d02tu;
  Z[4] = -7.0 * c02t - 10.0 * Z[4];
  Z[4] = xu * Z[4];
  Z[2] = 4.0 * Z[2] + Z[4];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] + Z[3];
  Z[2] = xt * Z[2];
  Z[1] = Z[1] + Z[2];
  Z[1] = xt * Z[1];
  Z[2] = 4.0 * c02u - Z[5];
  Z[2] = Z[2] * Z[20];
  Z[1] = 3.0 * Z[2] + Z[1];
  Z[1] = Z[1] * xt * xt * xt;
  Z[0] = 4.0 * Z[1] + Z[0];
  Z[0] = xs * Z[0];
  Z[1] = -c02t - c02u;
  Z[1] = Z[1] * Z[15];
  Z[2] = -2.0 * c02t - Z[11];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] - Z[12];
  Z[2] = xt * Z[2];
  Z[1] = 2.0 * Z[1] + Z[2];
  Z[1] = xt * Z[1];
  Z[2] = -Z[20] * Z[14];
  Z[1] = Z[2] + Z[1];
  Z[1] = Z[1] * xt * xt * xt * xt * xt;
  Z[0] = 8.0 * Z[1] + Z[0];

  bstu = M_f64 * M_f64 * Z[0];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * bstu /
         (48.0 * xs * xs * xs * xt * xt * xt * xu * xu * xu *
          (E1 * E2 * E32 * E42 * s * s * t * t * u));
}

void BPHOAMPFFTSU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64);
  complex<double> b02s = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> btsu;
  complex<double> Z[28];

  Z[0] = 2.0 * d02st;
  Z[1] = 3.0 * c02s;
  Z[2] = Z[0] + Z[1];
  Z[3] = 3.0 * xu;
  Z[4] = c02t + d02st;
  Z[5] = c02s - Z[4];
  Z[5] = Z[5] * Z[3];
  Z[2] = 4.0 * Z[2] + Z[5];
  Z[5] = d02st * Z[3];
  Z[6] = -xs * Z[0];
  Z[4] = Z[6] + 4.0 * Z[4] + Z[5];
  Z[4] = xs * Z[4];
  Z[2] = 2.0 * Z[2] + Z[4];
  Z[4] = 2.0 * xs;
  Z[2] = Z[2] * Z[4];
  Z[5] = xu * d02tu;
  Z[6] = 37.0 * Z[5];
  Z[7] = 26.0 * d02tu - 33.0 * c02t;
  Z[7] = 2.0 * Z[7] + Z[6];
  Z[8] = xu * xu;
  Z[7] = Z[7] * Z[8];
  Z[2] = Z[7] + Z[2];
  Z[2] = xs * Z[2];
  Z[7] = b02t + b02u;
  Z[9] = 3.0 * c02t;
  Z[10] = Z[9] + Z[7];
  Z[11] = 2.0 * d02tu;
  Z[12] = Z[11] - Z[9];
  Z[13] = 7.0 * c02u;
  Z[12] = 4.0 * Z[12] - Z[13];
  Z[12] = 5.0 * Z[12] + 7.0 * Z[5];
  Z[12] = xu * Z[12];
  Z[10] = 8.0 * Z[10] + Z[12];
  Z[12] = 2.0 * Z[8];
  Z[10] = Z[10] * Z[12];
  Z[2] = Z[10] + Z[2];
  Z[2] = xs * Z[2];
  Z[10] = 3.0 * c02u;
  Z[14] = -Z[5] - 4.0 * c02t - Z[10];
  Z[15] = xu * xu * xu;
  Z[16] = 4.0 * Z[15];
  Z[14] = Z[14] * Z[16];
  Z[17] = 10.0 * d02tu;
  Z[18] = -4.0 * Z[5] + Z[17] - 19.0 * c02t;
  Z[12] = Z[18] * Z[12];
  Z[18] = Z[8] * xs;
  Z[19] = Z[18] * d02tu;
  Z[12] = Z[12] - 7.0 * Z[19];
  Z[12] = xs * Z[12];
  Z[20] = -c02t - Z[5];
  Z[21] = xt * d02tu;
  Z[20] = -4.0 * Z[21] + 8.0 * Z[20];
  Z[20] = Z[8] * Z[20];
  Z[19] = -11.0 * Z[19] + Z[20];
  Z[19] = xt * Z[19];
  Z[12] = Z[19] + Z[14] + Z[12];
  Z[12] = xt * Z[12];
  Z[14] = 2.0 * c02u;
  Z[19] = Z[14] + c02t;
  Z[20] = xu * xu * xu * xu;
  Z[21] = Z[19] * Z[20];
  Z[2] = 2.0 * Z[12] - 16.0 * Z[21] + Z[2];
  Z[2] = xt * Z[2];
  Z[12] = 9.0 * c02s;
  Z[21] = -Z[9] + 4.0 * d02st + Z[12];
  Z[22] = xu + 1.0;
  Z[22] = d02st * Z[22];
  Z[22] = c02s + Z[22];
  Z[22] = Z[22] * Z[3];
  Z[23] = Z[9] + 7.0 * d02st + c02s;
  Z[24] = xs * d02st;
  Z[23] = 2.0 * Z[23] - 3.0 * Z[24];
  Z[23] = xs * Z[23];
  Z[21] = Z[23] + 4.0 * Z[21] + Z[22];
  Z[21] = xs * Z[21];
  Z[22] = Z[1] - Z[9];
  Z[23] = b02t - b02s;
  Z[24] = Z[22] + Z[23];
  Z[25] = -Z[0] + c02s;
  Z[25] = 2.0 * Z[25] - 5.0 * c02t;
  Z[25] = 6.0 * Z[25] + Z[6];
  Z[25] = xu * Z[25];
  Z[24] = 16.0 * Z[24] + Z[25];
  Z[24] = xu * Z[24];
  Z[21] = Z[24] + 4.0 * Z[21];
  Z[21] = xs * Z[21];
  Z[24] = 5.0 * b02t;
  Z[25] = 8.0 * d02tu;
  Z[26] = 3.0 * b02u + b02s;
  Z[26] = 24.0 * c02t - Z[25] + 3.0 * Z[26] + Z[24];
  Z[27] = 21.0 * c02u;
  Z[6] = Z[6] - 46.0 * c02t - Z[27];
  Z[6] = xu * Z[6];
  Z[6] = 4.0 * Z[26] + Z[6];
  Z[6] = xu * Z[6];
  Z[6] = -32.0 * b02s + Z[6];
  Z[6] = xu * Z[6];
  Z[6] = Z[6] + Z[21];
  Z[6] = xs * Z[6];
  Z[21] = 9.0 * c02u;
  Z[26] = 6.0 * c02t + Z[7];
  Z[26] = 2.0 * Z[26] + Z[21];
  Z[17] = -Z[27] + Z[17] - 7.0 * c02t;
  Z[17] = xu * Z[17];
  Z[17] = 2.0 * Z[26] + Z[17];
  Z[16] = Z[17] * Z[16];
  Z[6] = Z[16] + Z[6];
  Z[6] = xs * Z[6];
  Z[16] = c02u * xu * xu * xu * xu * xu;
  Z[2] = Z[2] - 8.0 * Z[16] + Z[6];
  Z[2] = xt * Z[2];
  Z[6] = 2.0 * Z[23] + Z[22];
  Z[17] = Z[0] - d02su;
  Z[11] = Z[11] + Z[17];
  Z[11] = xu * Z[11];
  Z[17] = Z[17] + c02t;
  Z[17] = -c02s - 2.0 * Z[17];
  Z[11] = 2.0 * Z[17] + Z[11];
  Z[11] = Z[11] * Z[3];
  Z[6] = 16.0 * Z[6] + Z[11];
  Z[6] = xu * Z[6];
  Z[11] = d02st + Z[1];
  Z[11] = 2.0 * Z[11] - Z[9];
  Z[0] = Z[0] + c02s;
  Z[17] = -d02su + d02st;
  Z[17] = xu * Z[17];
  Z[17] = 2.0 * Z[0] + Z[17];
  Z[17] = Z[17] * Z[3];
  Z[0] = xs * Z[0];
  Z[0] = 12.0 * Z[0] + 8.0 * Z[11] + Z[17];
  Z[0] = Z[0] * Z[4];
  Z[0] = Z[6] + Z[0];
  Z[0] = xs * Z[0];
  Z[4] = -Z[9] + 3.0 * d02su - 23.0 * d02tu;
  Z[4] = 3.0 * Z[5] + 2.0 * Z[4] + 5.0 * c02u;
  Z[4] = xu * Z[4];
  Z[5] = 1.0 + 6.0 * b02u;
  Z[5] = 12.0 * c02t + 10.0 * b02t + 2.0 * Z[5] - 5.0 * b02s;
  Z[4] = 2.0 * Z[5] + Z[4];
  Z[4] = xu * Z[4];
  Z[5] = -b02s - b02t;
  Z[4] = 16.0 * Z[5] + Z[4];
  Z[4] = xu * Z[4];
  Z[0] = 2.0 * Z[4] + Z[0];
  Z[0] = xs * Z[0];
  Z[4] = 2.0 * b02s;
  Z[5] = -Z[4] - Z[7];
  Z[6] = 18.0 * c02u + 36.0 * c02t - 16.0 * d02tu + Z[24] + b02u + 3.0 * b02s;
  Z[9] = -92.0 * d02tu - c02u;
  Z[9] = xu * Z[9];
  Z[6] = 4.0 * Z[6] + Z[9];
  Z[6] = xu * Z[6];
  Z[5] = 32.0 * Z[5] + Z[6];
  Z[5] = Z[5] * Z[8];
  Z[0] = Z[5] + Z[0];
  Z[0] = xs * Z[0];
  Z[5] = xu * c02u;
  Z[6] = 24.0 * Z[19] - 7.0 * Z[5];
  Z[6] = Z[6] * Z[20];
  Z[0] = 2.0 * Z[6] + Z[0];
  Z[0] = xs * Z[0];
  Z[0] = Z[0] + Z[2];
  Z[0] = xt * Z[0];
  Z[2] = -8.0 * d02su + 3.0 * b02t;
  Z[4] = 12.0 * c02s + 7.0 * b02u + Z[4] + Z[2];
  Z[6] = 4.0 * d02su;
  Z[9] = -Z[6] - Z[12];
  Z[3] = -d02su * Z[3];
  Z[3] = Z[3] + 2.0 * Z[9] - 39.0 * c02u;
  Z[3] = xu * Z[3];
  Z[3] = 4.0 * Z[4] + Z[3];
  Z[3] = xu * Z[3];
  Z[3] = -32.0 * b02t + Z[3];
  Z[3] = xu * Z[3];
  Z[4] = Z[6] - Z[1];
  Z[9] = xu * d02su;
  Z[4] = 2.0 * Z[4] - Z[9];
  Z[4] = Z[4] * Z[8];
  Z[11] = Z[18] * d02su;
  Z[4] = 5.0 * Z[4] - 14.0 * Z[11];
  Z[4] = xs * Z[4];
  Z[3] = Z[3] + Z[4];
  Z[3] = xs * Z[3];
  Z[2] = -Z[25] + 9.0 * b02u - 4.0 * b02s + Z[2] + 24.0 * c02u;
  Z[2] = 4.0 * Z[2] - 27.0 * Z[5];
  Z[2] = xu * Z[2];
  Z[4] = -2.0 * b02t - b02u - b02s;
  Z[2] = 32.0 * Z[4] + Z[2];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] + Z[3];
  Z[2] = xs * Z[2];
  Z[3] = -b02s - Z[7];
  Z[4] = -4.0 * d02tu + Z[21];
  Z[4] = xu * Z[4];
  Z[3] = 4.0 * Z[3] + Z[4];
  Z[3] = Z[3] * Z[15];
  Z[2] = 8.0 * Z[3] + Z[2];
  Z[2] = xs * Z[2];
  Z[2] = 24.0 * Z[16] + Z[2];
  Z[2] = xs * Z[2];
  Z[0] = Z[2] + Z[0];
  Z[0] = xt * Z[0];
  Z[2] = Z[1] + b02s;
  Z[3] = Z[10] + Z[2];
  Z[1] = -Z[1] - Z[13];
  Z[1] = xu * Z[1];
  Z[1] = 4.0 * Z[3] + Z[1];
  Z[1] = Z[1] * Z[15];
  Z[3] = d02su - c02s;
  Z[3] = 10.0 * Z[3] - Z[9];
  Z[3] = Z[3] * Z[8];
  Z[4] = -Z[18] * Z[6];
  Z[3] = Z[3] + Z[4];
  Z[3] = xs * Z[3];
  Z[4] = c02u - d02su;
  Z[4] = -7.0 * c02s - 10.0 * Z[4];
  Z[4] = xu * Z[4];
  Z[2] = 4.0 * Z[2] + Z[4];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] + Z[3];
  Z[2] = xs * Z[2];
  Z[1] = Z[1] + Z[2];
  Z[1] = xs * Z[1];
  Z[2] = 4.0 * c02u - Z[5];
  Z[2] = Z[2] * Z[20];
  Z[1] = 3.0 * Z[2] + Z[1];
  Z[1] = Z[1] * xs * xs * xs;
  Z[0] = 4.0 * Z[1] + Z[0];
  Z[0] = xt * Z[0];
  Z[1] = -c02s - c02u;
  Z[1] = Z[1] * Z[15];
  Z[2] = -2.0 * c02s - Z[9];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] - Z[11];
  Z[2] = xs * Z[2];
  Z[1] = 2.0 * Z[1] + Z[2];
  Z[1] = xs * Z[1];
  Z[2] = -Z[20] * Z[14];
  Z[1] = Z[2] + Z[1];
  Z[1] = Z[1] * xs * xs * xs * xs * xs;
  Z[0] = 8.0 * Z[1] + Z[0];

  btsu = M_f64 * M_f64 * Z[0];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * btsu /
         (48.0 * xs * xs * xs * xt * xt * xt * xu * xu * xu *
          (E12 * E2 * E3 * E42 * s * s * t * t * u));
}

void BPHOAMPFFUST_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_B0cc(&rslt, &s_f64, &M2_f64, &M2_f64);
  complex<double> b02s = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> bust;
  complex<double> Z[28];

  Z[0] = 2.0 * d02su;
  Z[1] = 3.0 * c02s;
  Z[2] = Z[0] + Z[1];
  Z[3] = 3.0 * xt;
  Z[4] = c02u + d02su;
  Z[5] = c02s - Z[4];
  Z[5] = Z[5] * Z[3];
  Z[2] = 4.0 * Z[2] + Z[5];
  Z[5] = d02su * Z[3];
  Z[6] = -xs * Z[0];
  Z[4] = Z[6] + 4.0 * Z[4] + Z[5];
  Z[4] = xs * Z[4];
  Z[2] = 2.0 * Z[2] + Z[4];
  Z[4] = 2.0 * xs;
  Z[2] = Z[2] * Z[4];
  Z[5] = xt * d02tu;
  Z[6] = 37.0 * Z[5];
  Z[7] = 26.0 * d02tu - 33.0 * c02u;
  Z[7] = 2.0 * Z[7] + Z[6];
  Z[8] = xt * xt;
  Z[7] = Z[7] * Z[8];
  Z[2] = Z[7] + Z[2];
  Z[2] = xs * Z[2];
  Z[7] = 8.0 * d02tu;
  Z[9] = 7.0 * c02t;
  Z[10] = 12.0 * c02u;
  Z[11] = -Z[10] + Z[7] - Z[9];
  Z[11] = 5.0 * Z[11] + 7.0 * Z[5];
  Z[11] = xt * Z[11];
  Z[12] = b02u + b02t;
  Z[13] = 3.0 * c02u;
  Z[14] = Z[13] + Z[12];
  Z[11] = 8.0 * Z[14] + Z[11];
  Z[14] = 2.0 * Z[8];
  Z[11] = Z[11] * Z[14];
  Z[2] = Z[11] + Z[2];
  Z[2] = xs * Z[2];
  Z[11] = 3.0 * c02t;
  Z[15] = -Z[5] - Z[11] - 4.0 * c02u;
  Z[16] = xt * xt * xt;
  Z[17] = 4.0 * Z[16];
  Z[15] = Z[15] * Z[17];
  Z[18] = 10.0 * d02tu;
  Z[19] = -4.0 * Z[5] + Z[18] - 19.0 * c02u;
  Z[14] = Z[19] * Z[14];
  Z[19] = Z[8] * xs;
  Z[20] = Z[19] * d02tu;
  Z[14] = Z[14] - 7.0 * Z[20];
  Z[14] = xs * Z[14];
  Z[21] = -c02u - Z[5];
  Z[22] = xu * d02tu;
  Z[21] = -4.0 * Z[22] + 8.0 * Z[21];
  Z[21] = Z[8] * Z[21];
  Z[20] = -11.0 * Z[20] + Z[21];
  Z[20] = xu * Z[20];
  Z[14] = Z[20] + Z[15] + Z[14];
  Z[14] = xu * Z[14];
  Z[15] = 2.0 * c02t;
  Z[20] = Z[15] + c02u;
  Z[21] = xt * xt * xt * xt;
  Z[22] = Z[20] * Z[21];
  Z[2] = 2.0 * Z[14] - 16.0 * Z[22] + Z[2];
  Z[2] = xu * Z[2];
  Z[14] = 9.0 * c02s;
  Z[22] = -Z[13] + 4.0 * d02su + Z[14];
  Z[23] = xt + 1.0;
  Z[23] = d02su * Z[23];
  Z[23] = c02s + Z[23];
  Z[23] = Z[23] * Z[3];
  Z[24] = Z[13] + 7.0 * d02su + c02s;
  Z[25] = xs * d02su;
  Z[24] = 2.0 * Z[24] - 3.0 * Z[25];
  Z[24] = xs * Z[24];
  Z[22] = Z[24] + 4.0 * Z[22] + Z[23];
  Z[22] = xs * Z[22];
  Z[23] = Z[13] - Z[1];
  Z[24] = b02u - b02s;
  Z[25] = -Z[23] + Z[24];
  Z[26] = -Z[0] + c02s;
  Z[26] = 2.0 * Z[26] - 5.0 * c02u;
  Z[26] = 6.0 * Z[26] + Z[6];
  Z[26] = xt * Z[26];
  Z[25] = 16.0 * Z[25] + Z[26];
  Z[25] = xt * Z[25];
  Z[22] = Z[25] + 4.0 * Z[22];
  Z[22] = xs * Z[22];
  Z[25] = 5.0 * b02u;
  Z[26] = 3.0 * b02t + b02s;
  Z[26] = 24.0 * c02u - Z[7] + 3.0 * Z[26] + Z[25];
  Z[27] = 21.0 * c02t;
  Z[6] = Z[6] - Z[27] - 46.0 * c02u;
  Z[6] = xt * Z[6];
  Z[6] = 4.0 * Z[26] + Z[6];
  Z[6] = xt * Z[6];
  Z[6] = -32.0 * b02s + Z[6];
  Z[6] = xt * Z[6];
  Z[6] = Z[6] + Z[22];
  Z[6] = xs * Z[6];
  Z[22] = 9.0 * c02t;
  Z[26] = Z[10] + 2.0 * Z[12] + Z[22];
  Z[18] = -7.0 * c02u + Z[18] - Z[27];
  Z[18] = xt * Z[18];
  Z[18] = 2.0 * Z[26] + Z[18];
  Z[17] = Z[18] * Z[17];
  Z[6] = Z[17] + Z[6];
  Z[6] = xs * Z[6];
  Z[17] = c02t * xt * xt * xt * xt * xt;
  Z[2] = Z[2] - 8.0 * Z[17] + Z[6];
  Z[2] = xu * Z[2];
  Z[6] = 2.0 * Z[24] - Z[23];
  Z[18] = Z[0] - d02st;
  Z[23] = Z[18] + c02u;
  Z[23] = -c02s - 2.0 * Z[23];
  Z[18] = 2.0 * d02tu + Z[18];
  Z[18] = xt * Z[18];
  Z[18] = 2.0 * Z[23] + Z[18];
  Z[18] = Z[18] * Z[3];
  Z[6] = 16.0 * Z[6] + Z[18];
  Z[6] = xt * Z[6];
  Z[18] = d02su + Z[1];
  Z[13] = 2.0 * Z[18] - Z[13];
  Z[0] = Z[0] + c02s;
  Z[18] = d02su - d02st;
  Z[18] = xt * Z[18];
  Z[18] = 2.0 * Z[0] + Z[18];
  Z[18] = Z[18] * Z[3];
  Z[0] = xs * Z[0];
  Z[0] = 12.0 * Z[0] + 8.0 * Z[13] + Z[18];
  Z[0] = Z[0] * Z[4];
  Z[0] = Z[6] + Z[0];
  Z[0] = xs * Z[0];
  Z[4] = 1.0 + 6.0 * b02t;
  Z[4] = Z[10] + 10.0 * b02u + 2.0 * Z[4] - 5.0 * b02s;
  Z[6] = 3.0 * d02st - 23.0 * d02tu;
  Z[5] = 3.0 * Z[5] - 6.0 * c02u + 2.0 * Z[6] + 5.0 * c02t;
  Z[5] = xt * Z[5];
  Z[4] = 2.0 * Z[4] + Z[5];
  Z[4] = xt * Z[4];
  Z[5] = -b02s - b02u;
  Z[4] = 16.0 * Z[5] + Z[4];
  Z[4] = xt * Z[4];
  Z[0] = 2.0 * Z[4] + Z[0];
  Z[0] = xs * Z[0];
  Z[4] = 2.0 * b02s;
  Z[5] = -Z[4] - Z[12];
  Z[6] = 36.0 * c02u + 18.0 * c02t - 16.0 * d02tu + Z[25] + b02t + 3.0 * b02s;
  Z[10] = -92.0 * d02tu - c02t;
  Z[10] = xt * Z[10];
  Z[6] = 4.0 * Z[6] + Z[10];
  Z[6] = xt * Z[6];
  Z[5] = 32.0 * Z[5] + Z[6];
  Z[5] = Z[5] * Z[8];
  Z[0] = Z[5] + Z[0];
  Z[0] = xs * Z[0];
  Z[5] = xt * c02t;
  Z[6] = 24.0 * Z[20] - 7.0 * Z[5];
  Z[6] = Z[6] * Z[21];
  Z[0] = 2.0 * Z[6] + Z[0];
  Z[0] = xs * Z[0];
  Z[0] = Z[0] + Z[2];
  Z[0] = xu * Z[0];
  Z[2] = -8.0 * d02st + 3.0 * b02u;
  Z[4] = 12.0 * c02s + 7.0 * b02t + Z[4] + Z[2];
  Z[6] = 4.0 * d02st;
  Z[10] = -Z[6] - Z[14];
  Z[3] = -d02st * Z[3];
  Z[3] = Z[3] + 2.0 * Z[10] - 39.0 * c02t;
  Z[3] = xt * Z[3];
  Z[3] = 4.0 * Z[4] + Z[3];
  Z[3] = xt * Z[3];
  Z[3] = -32.0 * b02u + Z[3];
  Z[3] = xt * Z[3];
  Z[4] = Z[6] - Z[1];
  Z[6] = xt * d02st;
  Z[4] = 2.0 * Z[4] - Z[6];
  Z[4] = Z[4] * Z[8];
  Z[10] = Z[19] * d02st;
  Z[4] = 5.0 * Z[4] - 14.0 * Z[10];
  Z[4] = xs * Z[4];
  Z[3] = Z[3] + Z[4];
  Z[3] = xs * Z[3];
  Z[2] = -Z[7] + 9.0 * b02t - 4.0 * b02s + Z[2] + 24.0 * c02t;
  Z[2] = 4.0 * Z[2] - 27.0 * Z[5];
  Z[2] = xt * Z[2];
  Z[4] = -2.0 * b02u - b02t - b02s;
  Z[2] = 32.0 * Z[4] + Z[2];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] + Z[3];
  Z[2] = xs * Z[2];
  Z[3] = -b02s - Z[12];
  Z[4] = -4.0 * d02tu + Z[22];
  Z[4] = xt * Z[4];
  Z[3] = 4.0 * Z[3] + Z[4];
  Z[3] = Z[3] * Z[16];
  Z[2] = 8.0 * Z[3] + Z[2];
  Z[2] = xs * Z[2];
  Z[2] = 24.0 * Z[17] + Z[2];
  Z[2] = xs * Z[2];
  Z[0] = Z[2] + Z[0];
  Z[0] = xu * Z[0];
  Z[2] = Z[1] + b02s;
  Z[3] = Z[11] + Z[2];
  Z[1] = -Z[1] - Z[9];
  Z[1] = xt * Z[1];
  Z[1] = 4.0 * Z[3] + Z[1];
  Z[1] = Z[1] * Z[16];
  Z[3] = c02t - d02st;
  Z[3] = -7.0 * c02s - 10.0 * Z[3];
  Z[3] = xt * Z[3];
  Z[2] = 4.0 * Z[2] + Z[3];
  Z[2] = Z[2] * Z[8];
  Z[3] = d02st - c02s;
  Z[3] = 10.0 * Z[3] - Z[6];
  Z[3] = Z[3] * Z[8];
  Z[3] = Z[3] - 4.0 * Z[10];
  Z[3] = xs * Z[3];
  Z[2] = Z[2] + Z[3];
  Z[2] = xs * Z[2];
  Z[1] = Z[1] + Z[2];
  Z[1] = xs * Z[1];
  Z[2] = 4.0 * c02t - Z[5];
  Z[2] = Z[2] * Z[21];
  Z[1] = 3.0 * Z[2] + Z[1];
  Z[1] = Z[1] * xs * xs * xs;
  Z[0] = 4.0 * Z[1] + Z[0];
  Z[0] = xu * Z[0];
  Z[1] = -c02s - c02t;
  Z[1] = Z[1] * Z[16];
  Z[2] = -2.0 * c02s - Z[6];
  Z[2] = Z[2] * Z[8];
  Z[2] = Z[2] - Z[10];
  Z[2] = xs * Z[2];
  Z[1] = 2.0 * Z[1] + Z[2];
  Z[1] = xs * Z[1];
  Z[2] = -Z[21] * Z[15];
  Z[1] = Z[2] + Z[1];
  Z[1] = Z[1] * xs * xs * xs * xs * xs;
  Z[0] = 8.0 * Z[1] + Z[0];

  bust = M_f64 * M_f64 * Z[0];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * bust /
         (48.0 * xs * xs * xs * xt * xt * xt * xu * xu * xu *
          (E1 * E22 * E3 * E42 * s * s * t * u * u));
}

void CPHOAMPFFSTU_f64(complex<double> E1, complex<double> E2,
                      complex<double> E3, complex<double> p1_p2,
                      complex<double> p1_p3, complex<double> p2_p3,
                      complex<double> *out) {
  complex<double> s, t, u, E4, E12, E22, E32, E42;
  complex<double> s_f64, t_f64, u_f64;
  complex<double> xs, xt, xu;

  s = 2.0 * p1_p2;
  t = 2.0 * p2_p3;
  u = 2.0 * p1_p3;
  E4 = -E1 - E2 - E3;
  E12 = E1 * E1;
  E22 = E2 * E2;
  E32 = E3 * E3;
  E42 = E4 * E4;

  s_f64 = (complex<double>)s;
  t_f64 = (complex<double>)t;
  u_f64 = (complex<double>)u;
  OLO_SCALE(&scale);

  xs = s / (M_f64 * M_f64);
  xt = t / (M_f64 * M_f64);
  xu = u / (M_f64 * M_f64);

  complex<double> M2_f64 = M_f64 * M_f64;

  complex<double> rslt[3];
  OLO_A0c(&rslt, &M2_f64);
  complex<double> a02 = (complex<double>)rslt[0] * 1.0 / (M_f64 * M_f64);

  OLO_B0cc(&rslt, &t_f64, &M2_f64, &M2_f64);
  complex<double> b02t = (complex<double>)rslt[0] * 1.0;

  OLO_B0cc(&rslt, &u_f64, &M2_f64, &M2_f64);
  complex<double> b02u = (complex<double>)rslt[0] * 1.0;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &s_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02s = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &t_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02t = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_C0cc(&rslt, &zero_f64, &zero_f64, &u_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> c02u = (complex<double>)rslt[0] * M_f64 * M_f64;

  OLO_B0cc(&rslt, &zero_f64, &M2_f64, &M2_f64);
  complex<double> b020 = (complex<double>)rslt[0] * 1.0;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02su =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &s_f64, &t_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02st =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  OLO_D0cc(&rslt, &zero_f64, &zero_f64, &zero_f64, &zero_f64, &t_f64, &u_f64,
           &M2_f64, &M2_f64, &M2_f64, &M2_f64);
  complex<double> d02tu =
      (complex<double>)rslt[0] * M_f64 * M_f64 * M_f64 * M_f64;

  complex<double> cstu;
  complex<double> Z[35];

  Z[0] = 3.0 * xu;
  Z[1] = -d02su - c02u;
  Z[1] = Z[1] * Z[0];
  Z[2] = 3.0 * c02u;
  Z[3] = Z[2] + b02t;
  Z[4] = 3.0 * d02su;
  Z[5] = a02 - 1.0;
  Z[6] = -c02t - Z[5];
  Z[1] = Z[1] - Z[4] + 3.0 * Z[6] + 7.0 * b02u - Z[3];
  Z[1] = xu * Z[1];
  Z[6] = 3.0 * b02u;
  Z[7] = Z[6] - a02;
  Z[1] = Z[1] - 2.0 - Z[7];
  Z[8] = xu * xu;
  Z[9] = 2.0 * Z[8];
  Z[1] = Z[1] * Z[9];
  Z[10] = 3.0 * b02t;
  Z[11] = 4.0 * c02s;
  Z[12] = 2.0 * c02t;
  Z[13] = -Z[11] - Z[12] + Z[10];
  Z[13] = xu * Z[13];
  Z[14] = xu + xt;
  Z[14] = d02st * Z[14];
  Z[15] = 2.0 * c02s;
  Z[16] = d02st - c02t;
  Z[16] = 3.0 * Z[16] - Z[15] + Z[14];
  Z[17] = 2.0 * xt;
  Z[16] = Z[16] * Z[17];
  Z[13] = Z[13] + Z[16];
  Z[13] = xt * Z[13];
  Z[16] = xu * c02u;
  Z[18] = b02u - Z[5] + b02t;
  Z[18] = 2.0 * Z[18] + Z[16];
  Z[18] = xu * Z[18];
  Z[19] = 8.0 * b02t;
  Z[18] = Z[19] + Z[18];
  Z[18] = xu * Z[18];
  Z[13] = Z[18] + Z[13];
  Z[18] = 3.0 * xt;
  Z[13] = Z[13] * Z[18];
  Z[1] = Z[1] + Z[13];
  Z[1] = xt * Z[1];
  Z[13] = 2.0 * b020;
  Z[19] = -Z[19] + Z[13] - Z[6];
  Z[20] = 4.0 * c02u;
  Z[21] = -Z[4] - Z[20];
  Z[21] = Z[21] * Z[0];
  Z[21] = Z[21] + 30.0 * c02u - 18.0 * d02su + 13.0 * b02u + 12.0 * c02s;
  Z[21] = xu * Z[21];
  Z[19] = 2.0 * Z[19] + Z[21];
  Z[21] = xu * xu * xu;
  Z[19] = Z[19] * Z[21];
  Z[1] = Z[19] + Z[1];
  Z[1] = xt * Z[1];
  Z[19] = Z[2] + c02s;
  Z[22] = 4.0 * Z[19] - Z[16];
  Z[23] = xu * xu * xu * xu * xu;
  Z[22] = Z[22] * Z[23];
  Z[1] = 3.0 * Z[22] + Z[1];
  Z[1] = xt * Z[1];
  Z[22] = xu * d02su;
  Z[24] = Z[22] + Z[2] - Z[11] + 7.0 * d02su;
  Z[24] = xu * Z[24];
  Z[25] = Z[5] - Z[6];
  Z[26] = Z[25] - b02t;
  Z[24] = Z[24] - Z[26];
  Z[24] = xu * Z[24];
  Z[24] = 10.0 * b02t + Z[24];
  Z[24] = xu * Z[24];
  Z[27] = -Z[15] - 4.0 * c02t + Z[10];
  Z[27] = xu * Z[27];
  Z[28] = -c02s + 3.0 * d02st - 8.0 * c02t;
  Z[28] = 2.0 * Z[28] + Z[14];
  Z[28] = xt * Z[28];
  Z[27] = Z[27] + Z[28];
  Z[27] = xt * Z[27];
  Z[24] = Z[24] + Z[27];
  Z[24] = Z[24] * Z[18];
  Z[27] = 3.0 * c02s;
  Z[28] = -Z[27] + 4.0 * d02su;
  Z[28] = Z[22] + 2.0 * Z[28] - 9.0 * c02u;
  Z[28] = Z[28] * Z[0];
  Z[29] = 6.0 * c02s;
  Z[30] = 5.0 * b02u + Z[29];
  Z[28] = Z[28] - Z[2] + 5.0 * Z[30] - 12.0 * d02su;
  Z[28] = xu * Z[28];
  Z[28] = Z[13] + Z[28];
  Z[28] = Z[28] * Z[8];
  Z[24] = Z[28] + Z[24];
  Z[24] = xt * Z[24];
  Z[28] = 5.0 * c02u;
  Z[30] = Z[28] + Z[29];
  Z[31] = 2.0 * d02su;
  Z[28] = -Z[28] - c02s + Z[31];
  Z[28] = xu * Z[28];
  Z[28] = Z[28] + b02u + Z[30];
  Z[32] = xu * xu * xu * xu;
  Z[28] = Z[28] * Z[32];
  Z[24] = 6.0 * Z[28] + Z[24];
  Z[24] = xt * Z[24];
  Z[27] = -Z[27] + Z[31];
  Z[27] = -5.0 * Z[22] + 4.0 * Z[27] + c02u;
  Z[27] = Z[27] * Z[21];
  Z[4] = -Z[21] * Z[4];
  Z[28] = -xt * Z[12];
  Z[4] = Z[4] + Z[28];
  Z[4] = Z[4] * Z[17];
  Z[4] = Z[27] + Z[4];
  Z[4] = xt * Z[4];
  Z[27] = -xs * Z[31];
  Z[27] = Z[27] - Z[11] - 3.0 * Z[22];
  Z[27] = Z[21] * Z[27];
  Z[28] = 6.0 * d02su;
  Z[33] = Z[21] * xt;
  Z[34] = -Z[33] * Z[28];
  Z[27] = Z[34] + Z[27];
  Z[27] = xs * Z[27];
  Z[30] = -Z[22] - Z[30];
  Z[30] = Z[30] * Z[32];
  Z[4] = Z[27] + Z[30] + Z[4];
  Z[4] = xs * Z[4];
  Z[27] = -Z[8] * Z[31];
  Z[27] = Z[27] - Z[12] + b02t;
  Z[27] = xu * Z[27];
  Z[30] = d02st - 7.0 * c02t;
  Z[17] = Z[30] * Z[17];
  Z[17] = Z[27] + Z[17];
  Z[17] = xt * Z[17];
  Z[27] = c02u - Z[11] + 5.0 * d02su;
  Z[27] = 3.0 * Z[27] - Z[22];
  Z[27] = xu * Z[27];
  Z[30] = 2.0 * b02u;
  Z[27] = Z[30] + Z[27];
  Z[27] = xu * Z[27];
  Z[27] = 4.0 * b02t + Z[27];
  Z[27] = xu * Z[27];
  Z[17] = Z[27] + Z[17];
  Z[17] = xt * Z[17];
  Z[27] = c02u - d02su;
  Z[34] = -c02s - Z[27];
  Z[34] = xu * Z[34];
  Z[11] = 6.0 * Z[34] + b02u + Z[11];
  Z[34] = 2.0 * Z[21];
  Z[11] = Z[11] * Z[34];
  Z[11] = Z[11] + Z[17];
  Z[11] = xt * Z[11];
  Z[17] = Z[19] * Z[23];
  Z[4] = Z[4] - 2.0 * Z[17] + Z[11];
  Z[4] = xs * Z[4];
  Z[11] = xu * xu * xu * xu * xu * xu;
  Z[17] = Z[11] * Z[2];
  Z[4] = 3.0 * Z[4] - Z[17] + Z[24];
  Z[4] = xs * Z[4];
  Z[1] = Z[1] + Z[4];
  Z[1] = xs * Z[1];
  Z[4] = 6.0 * Z[22];
  Z[19] = 2.0 * b02t;
  Z[5] = Z[4] + Z[28] - 18.0 * c02s - Z[19] - 3.0 * Z[5] - Z[30];
  Z[5] = xu * Z[5];
  Z[6] = Z[19] + Z[6];
  Z[22] = 3.0 * b020 - a02 - Z[6];
  Z[5] = 2.0 * Z[22] + Z[5];
  Z[5] = Z[5] * Z[21];
  Z[22] = 3.0 - 2.0 * a02;
  Z[3] = Z[4] - Z[29] + 3.0 * Z[22] + b02u + Z[3];
  Z[3] = xu * Z[3];
  Z[4] = -3.0 + b020 - Z[7];
  Z[3] = 2.0 * Z[4] + Z[3];
  Z[3] = Z[3] * Z[8];
  Z[4] = -xu * Z[26];
  Z[4] = Z[19] + Z[4];
  Z[4] = xu * Z[4];
  Z[19] = d02st - c02s;
  Z[14] = 2.0 * Z[19] + Z[14];
  Z[14] = xt * Z[14];
  Z[15] = b02t - Z[15];
  Z[15] = xu * Z[15];
  Z[14] = Z[15] + Z[14];
  Z[14] = xt * Z[14];
  Z[4] = Z[4] + Z[14];
  Z[4] = Z[4] * Z[18];
  Z[3] = Z[3] + Z[4];
  Z[3] = xt * Z[3];
  Z[3] = Z[5] + Z[3];
  Z[3] = xt * Z[3];
  Z[4] = -c02s + Z[27];
  Z[4] = Z[4] * Z[0];
  Z[4] = Z[4] + b020 - Z[6];
  Z[5] = 2.0 * Z[32];
  Z[4] = Z[4] * Z[5];
  Z[3] = Z[4] + Z[3];
  Z[3] = xt * Z[3];
  Z[4] = c02u * Z[11];
  Z[3] = 6.0 * Z[4] + Z[3];
  Z[3] = xt * Z[3];
  Z[1] = Z[3] + Z[1];
  Z[1] = xs * Z[1];
  Z[3] = 2.0 * d02tu;
  Z[4] = -Z[20] + Z[31] + Z[3] + b02t;
  Z[4] = Z[4] * Z[0];
  Z[6] = Z[10] + 4.0 * b020 + Z[7];
  Z[4] = 2.0 * Z[6] + Z[4];
  Z[4] = Z[4] * Z[21];
  Z[6] = c02t - d02tu;
  Z[7] = b02u - Z[6];
  Z[11] = 2.0 * c02u;
  Z[14] = -d02tu + Z[11];
  Z[14] = xu * Z[14];
  Z[7] = Z[14] + 2.0 * Z[7] - b02t;
  Z[7] = Z[7] * Z[0];
  Z[14] = -1.0 + Z[13];
  Z[7] = 2.0 * Z[14] + Z[7];
  Z[7] = Z[7] * Z[8];
  Z[8] = Z[3] - c02t;
  Z[14] = -xu * Z[8];
  Z[14] = b02u + Z[14];
  Z[9] = Z[14] * Z[9];
  Z[14] = Z[33] * d02tu;
  Z[9] = Z[9] + Z[14];
  Z[9] = Z[9] * Z[18];
  Z[7] = Z[7] + Z[9];
  Z[7] = xt * Z[7];
  Z[4] = Z[4] + Z[7];
  Z[4] = xt * Z[4];
  Z[7] = Z[10] + Z[13] - Z[25];
  Z[2] = Z[31] - Z[2];
  Z[0] = Z[2] * Z[0];
  Z[0] = 2.0 * Z[7] + Z[0];
  Z[0] = Z[0] * Z[32];
  Z[0] = Z[0] + Z[4];
  Z[0] = xt * Z[0];
  Z[0] = -Z[17] + Z[0];
  Z[0] = Z[0] * xt * xt;
  Z[0] = Z[0] + Z[1];
  Z[0] = xs * Z[0];
  Z[1] = Z[12] + b02t;
  Z[2] = Z[16] - Z[11] - Z[1];
  Z[2] = Z[2] * Z[5];
  Z[4] = xu * d02tu;
  Z[6] = 4.0 * Z[6] + Z[4];
  Z[6] = Z[6] * Z[21];
  Z[3] = Z[33] * Z[3];
  Z[3] = Z[6] + Z[3];
  Z[3] = xt * Z[3];
  Z[6] = Z[11] - Z[8];
  Z[6] = xu * Z[6];
  Z[1] = Z[6] - Z[1];
  Z[1] = Z[1] * Z[34];
  Z[1] = Z[1] + Z[3];
  Z[1] = xt * Z[1];
  Z[1] = Z[2] + Z[1];
  Z[1] = xt * Z[1];
  Z[2] = -Z[23] * Z[20];
  Z[1] = Z[2] + Z[1];
  Z[1] = Z[1] * xt * xt * xt * xt;
  Z[0] = 3.0 * Z[1] + Z[0];
  Z[0] = xs * Z[0];
  Z[1] = c02t + c02u;
  Z[1] = Z[1] * Z[5];
  Z[2] = Z[12] + Z[4];
  Z[2] = Z[2] * Z[21];
  Z[2] = Z[2] + Z[14];
  Z[2] = xt * Z[2];
  Z[1] = Z[1] + Z[2];
  Z[1] = xt * Z[1];
  Z[2] = Z[23] * Z[11];
  Z[1] = Z[2] + Z[1];
  Z[1] = Z[1] * xt * xt * xt * xt * xt * xt;

  cstu = Z[0] + 3.0 * Z[1];

  *out = 2.0 * 3.0 * top_Q8_f64 /
         (M_f64 * M_f64 * M_f64 * M_f64 * M_PI * M_PI) * cstu /
         (3.0 * xs * xs * xs * xs * xt * xt * xt * xt * xu * xu * xu * xu *
          (E1 * E2 * E3 * E4 * s * t * t * u));
}
