#include "form_factors.h"
#include <mp++/complex.hpp>

// this is a constant in here for now, later this should be an argument of the functions
const complex<double> MT_f64 = 1.0;
const complex128 MT_f128 = 1.0;
const mppp::complex MT_mpfr = 1.0;

void FFS_f64(complex<double> p1_p1,complex<double> p2_p2, complex<double> p1_p2, complex<double>* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFS_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2, complex128* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFS_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2, mppp::complex* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFT_f64(complex<double> p1_p1,complex<double> p2_p2, complex<double> p1_p2, complex<double>* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFT_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2, complex128* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFT_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2, mppp::complex* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFU_f64(complex<double> p1_p1,complex<double> p2_p2, complex<double> p1_p2, complex<double>* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFU_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2, complex128* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}

void FFU_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2, mppp::complex* out) {
    *out = p1_p1 + p2_p2 + 2.0*p1_p2;
}




