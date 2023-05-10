#include <complex>
#include <mp++/complex128.hpp>
#include <mp++/complex.hpp>

using namespace std;
using complex128 = mppp::complex128;

extern void FFS_f64(complex<double>, complex<double>, complex<double>, complex<double>*);
extern void FFS_f128(complex128, complex128, complex128, complex128*);
extern void FFS_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void FFT_f64(complex<double>, complex<double>, complex<double>, complex<double>*);
extern void FFT_f128(complex128, complex128, complex128, complex128*);
extern void FFT_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void FFU_f64(complex<double>, complex<double>, complex<double>, complex<double>*);
extern void FFU_f128(complex128, complex128, complex128, complex128*);
extern void FFU_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void FFSINV_f64(complex<double>, complex<double>, complex<double>, complex<double>*);
extern void FFSINV_f128(complex128, complex128, complex128, complex128*);
extern void FFSINV_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void FFTINV_f64(complex<double>, complex<double>, complex<double>, complex<double>*);
extern void FFTINV_f128(complex128, complex128, complex128, complex128*);
extern void FFTINV_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void FFUINV_f64(complex<double>, complex<double>, complex<double>, complex<double>*);
extern void FFUINV_f128(complex128, complex128, complex128, complex128*);
extern void FFUINV_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex*);


extern void APHOAMPFFSTU_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void APHOAMPFFSTU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void APHOAMPFFSTU_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void APHOAMPFFTSU_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void APHOAMPFFTSU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void APHOAMPFFTSU_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void APHOAMPFFUST_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void APHOAMPFFUST_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void APHOAMPFFUST_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void BPHOAMPFFSTU_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void BPHOAMPFFSTU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void BPHOAMPFFSTU_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void BPHOAMPFFTSU_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void BPHOAMPFFTSU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void BPHOAMPFFTSU_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void BPHOAMPFFUST_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void BPHOAMPFFUST_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void BPHOAMPFFUST_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);

extern void CPHOAMPFFSTU_f64(complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>, complex<double>*);
extern void CPHOAMPFFSTU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128*);
extern void CPHOAMPFFSTU_mpfr(mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex, mppp::complex*);
