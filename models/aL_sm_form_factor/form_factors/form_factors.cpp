#include "form_factors.h"

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
    
    *out = 1;
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
    
    *out = 1;
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
    
    *out = 1;
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
    
    *out = 1;
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
    
    *out = 1;
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
    
    *out = 1;
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
    
    *out = 1;
}