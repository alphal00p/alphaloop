#include "form_factors_f64.h"
#include <complex>
#include <mp++/complex128.hpp>

typedef mppp::complex128 complex128;

using namespace std;
void APHOAMPFFSTU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  APHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}

void APHOAMPFFTSU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  APHOAMPFFTSU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void APHOAMPFFUST_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  APHOAMPFFUST_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void BPHOAMPFFSTU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  BPHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void BPHOAMPFFTSU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  BPHOAMPFFTSU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void BPHOAMPFFUST_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  BPHOAMPFFUST_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void CPHOAMPFFSTU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  CPHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
