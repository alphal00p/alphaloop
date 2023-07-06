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
    *out = p1_p1 + p2_p2 - 2.0*p1_p2;
}

void FFT_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2, complex128* out) {
    *out = p1_p1 + p2_p2 - 2.0*p1_p2;
}

void FFT_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2, mppp::complex* out) {
    *out = p1_p1 + p2_p2 - 2.0*p1_p2;
}

void FFU_f64(complex<double> p1_p1,complex<double> p2_p2, complex<double> p1_p2, complex<double>* out) {
    *out = p1_p1 + p2_p2 - 2.0*p1_p2;
}

void FFU_f128(complex128 p1_p1, complex128 p2_p2, complex128 p1_p2, complex128* out) {
    *out = p1_p1 + p2_p2 - 2.0*p1_p2;
}

void FFU_mpfr(mppp::complex p1_p1, mppp::complex p2_p2, mppp::complex p1_p2, mppp::complex* out) {
    *out = p1_p1 + p2_p2 - 2.0*p1_p2;
}

void APHOAMPFFSTU_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);

    complex<double> BO2, CO2, DO2;
    BO2 = 1.0;
    CO2 = 1.0;
    DO2 = 1.0;

    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;
    
    *out = res / (s*s);
}

void APHOAMPFFTSU_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);

    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;

    *out = res / (t*t);
}

void APHOAMPFFUST_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);

    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;

    *out = res / (u*u);
}

void BPHOAMPFFSTU_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);

    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;

    *out = res / (s*s*t);
}

void BPHOAMPFFTSU_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);
    
    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;

    *out = res / (t*t*s);
}

void BPHOAMPFFUST_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);
    
    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;

    *out = res / (u*u*s);
}

void CPHOAMPFFSTU_f64(complex<double> p1_p1, 
                      complex<double> p2_p2, 
                      complex<double> p3_p3, 
                      complex<double> p1_p2, 
                      complex<double> p1_p3, 
                      complex<double> p2_p3, 
                      complex<double>* out) {
    complex<double> s,t,u;
    FFS_f64(p1_p1, p2_p2, p1_p2,&s);
    FFT_f64(p1_p1, p3_p3, p1_p3, &t);
    FFU_f64(p2_p2, p3_p3, p2_p3, &u);
   
    complex<double> xs, xt, xu;
    xs = s / (MT_f64 * MT_f64);
    xt = t / (MT_f64 * MT_f64);
    xu = u / (MT_f64 * MT_f64);

    complex<double> res = 1.0;

    *out = res / (s*t*t*u);
}

void APHOAMPFFSTU_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);
   
    complex128 xs, xt, xu;
    xs = s / (MT_f128 * MT_f128);
    xt = t / (MT_f128 * MT_f128);
    xu = u / (MT_f128 * MT_f128);

    complex128 res = 1.0;

    *out = res / (s*s);
}

void APHOAMPFFTSU_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);
    
    complex128 res = 1.0;

    *out = res / (t*t);
}

void APHOAMPFFUST_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);
    
    complex128 res = 1.0;

    *out = res / (u*u);
}

void BPHOAMPFFSTU_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);

    complex128 res = 1.0;
    
    *out = res / (s*s*t);
}

void BPHOAMPFFTSU_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);
    
    complex128 res = 1.0;
    
    *out = res / (t*t*s);
}

void BPHOAMPFFUST_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);
     
    complex128 res = 1.0;
    
    *out = res / (u*u*s);
}

void CPHOAMPFFSTU_f128(complex128 p1_p1, 
                       complex128 p2_p2, 
                       complex128 p3_p3, 
                       complex128 p1_p2, 
                       complex128 p1_p3, 
                       complex128 p2_p3, 
                       complex128* out) {
    complex128 s,t,u;
    FFS_f128(p1_p1, p2_p2, p1_p2,&s);
    FFT_f128(p1_p1, p3_p3, p1_p3, &t);
    FFU_f128(p2_p2, p3_p3, p2_p3, &u);
    
    complex128 res = 1.0;

    *out = res / (s*t*t*u);
}

void APHOAMPFFSTU_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);

    mppp::complex res = 1.0;
    
    *out = res / (s*s);
}

void APHOAMPFFTSU_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);
    
    mppp::complex res = 1.0;
    *out = res / (t*t);
}

void APHOAMPFFUST_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);
    
    mppp::complex res = 1.0;
    *out = res / (u*u);
}

void BPHOAMPFFSTU_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);
    
    mppp::complex res = 1.0;

    *out = res / (s*s*t);
}

void BPHOAMPFFTSU_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);
    
    mppp::complex res = 1.0;

    *out = res / (t*t*s);
}

void BPHOAMPFFUST_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);
    
    mppp::complex res = 1.0;

    *out = res / (u*u*s);
}

void CPHOAMPFFSTU_mpfr(mppp::complex p1_p1, 
                       mppp::complex p2_p2, 
                       mppp::complex p3_p3, 
                       mppp::complex p1_p2, 
                       mppp::complex p1_p3, 
                       mppp::complex p2_p3, 
                       mppp::complex* out) {
    mppp::complex s,t,u;
    FFS_mpfr(p1_p1, p2_p2, p1_p2,&s);
    FFT_mpfr(p1_p1, p3_p3, p1_p3, &t);
    FFU_mpfr(p2_p2, p3_p3, p2_p3, &u);
    
    mppp::complex res = 1.0;

    *out = res / (s*t*t*u);
}