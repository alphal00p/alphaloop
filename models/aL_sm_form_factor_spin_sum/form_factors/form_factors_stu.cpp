#include "form_factors.h"
#include <cmath>
#include <complex>
#include <math.h>
#include <mp++/complex.hpp>

using namespace std;

#include "cavh_olo.h"

// this is a constant in here for now, later this should be an argument of the
// functions
complex<double> M_f64(173.0, 0.0);
complex128 M_f128(173.0, 0.0);
mppp::complex M_mpfr(173.0, 0.0);

complex<double> top_Q8_f64(16.0 / 81.0, 0);
complex128 top_Q8_f128(16.0 / 81.0, 0);
mppp::complex top_Q8_mpfr(16.0 / 81.0, 0);

double FORM_FACTOR_PHASE = 1.0;

complex<double> zero_f64(0.0, 0.0);
double scale = 1000.0;
int zero = 0;

void FFS_f64(complex<double> p1_p1, complex<double> p2_p2,
             complex<double> p1_p2, complex<double> *out) {
  *out = p1_p1 + p2_p2 + 2.0 * p1_p2;
}

void FFS_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2,
              complex128 *out) {
  *out = p1_p1 + p2_p2 + 2.0 * p1_p2;
}

void FFS_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2,
              mppp::complex *out) {
  *out = p1_p1 + p2_p2 + 2.0 * p1_p2;
}

void FFT_f64(complex<double> p1_p1, complex<double> p3_p3,
             complex<double> p1_p3, complex<double> *out) {
  *out = p1_p1 + p3_p3 + 2.0 * p1_p3;
}

void FFT_f128(complex128 p1_p1, complex128 p3_p3, complex128 p1_p3,
              complex128 *out) {
  *out = p1_p1 + p3_p3 + 2.0 * p1_p3;
}

void FFT_mpfr(mppp::complex p1_p1, mppp::complex p3_p3, mppp::complex p1_p3,
              mppp::complex *out) {
  *out = p1_p1 + p3_p3 + 2.0 * p1_p3;
}

void FFU_f64(complex<double> p1_p1, complex<double> p2_p2,
             complex<double> p1_p2, complex<double> *out) {
  *out = p1_p1 + p2_p2 + 2.0 * p1_p2;
}

void FFU_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2,
              complex128 *out) {
  *out = p1_p1 + p2_p2 + 2.0 * p1_p2;
}

void FFU_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2,
              mppp::complex *out) {
  *out = p1_p1 + p2_p2 + 2.0 * p1_p2;
}

void FFSINV_f64(complex<double> p1_p1, complex<double> p2_p2,
                complex<double> p1_p2, complex<double> *out) {
  *out = 1.0 / (p1_p1 + p2_p2 + 2.0 * p1_p2);
}

void FFSINV_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2,
                 complex128 *out) {
  *out = 1.0 / (p1_p1 + p2_p2 + 2.0 * p1_p2);
}

void FFSINV_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2,
                 mppp::complex *out) {
  *out = 1.0 / (p1_p1 + p2_p2 + 2.0 * p1_p2);
}

void FFTINV_f64(complex<double> p1_p1, complex<double> p3_p3,
                complex<double> p1_p3, complex<double> *out) {
  *out = 1.0 / (p1_p1 + p3_p3 + 2.0 * p1_p3);
}

void FFTINV_f128(complex128 p1_p1, complex128 p3_p3, complex128 p1_p3,
                 complex128 *out) {
  *out = 1.0 / (p1_p1 + p3_p3 + 2.0 * p1_p3);
}

void FFTINV_mpfr(mppp::complex p1_p1, mppp::complex p3_p3, mppp::complex p1_p3,
                 mppp::complex *out) {
  *out = 1.0 / (p1_p1 + p3_p3 + 2.0 * p1_p3);
}

void FFUINV_f64(complex<double> p1_p1, complex<double> p2_p2,
                complex<double> p1_p2, complex<double> *out) {
  *out = 1.0 / (p1_p1 + p2_p2 + 2.0 * p1_p2);
}

void FFUINV_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2,
                 complex128 *out) {
  *out = 1.0 / (p1_p1 + p2_p2 + 2.0 * p1_p2);
}

void FFUINV_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2,
                 mppp::complex *out) {
  *out = 1.0 / (p1_p1 + p2_p2 + 2.0 * p1_p2);
}