#include "form_factors_f64.h"
#include <complex>
#include <mp++/complex.hpp>

using namespace std;

void APHOAMPFFSTU_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  APHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}

void APHOAMPFFTSU_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  APHOAMPFFTSU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}

void APHOAMPFFUST_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  APHOAMPFFUST_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}

void BPHOAMPFFSTU_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  BPHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}

void BPHOAMPFFTSU_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  BPHOAMPFFTSU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}

void BPHOAMPFFUST_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  BPHOAMPFFUST_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}

void CPHOAMPFFSTU_mpfr(mppp::complex E1, mppp::complex E2, mppp::complex E3,
                       mppp::complex p1_p2, mppp::complex p1_p3,
                       mppp::complex p2_p3, mppp::complex *out) {
  complex<double> out_f64;
  CPHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (mppp::complex)out_f64;
}