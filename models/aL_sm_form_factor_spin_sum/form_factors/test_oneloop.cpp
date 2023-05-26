#include <complex>
#include <mp++/complex.hpp>
using namespace std;
#include "cavh_olo.h"
#include "form_factors.h"

complex<double> lorentz_dot(vector<complex<double>> p1,
                            vector<complex<double>> p2) {
  return p1[0] * p2[0] - p1[1] * p2[1] - p1[2] * p2[2] - p1[3] * p2[3];
}

complex<double> conjugate(complex<double> c) {
  complex<double> cstar(c.real(), -c.imag());
  return cstar;
}

int main() {
  vector<complex<double>> p1 = {complex<double>(5, 0), complex<double>(0, 0),
                                complex<double>(0, 0), complex<double>(5, 0)};

  vector<complex<double>> p2 = {complex<double>(5, 0), complex<double>(0, 0),
                                complex<double>(0, 0), complex<double>(-5, 0)};

  vector<complex<double>> p3 = {
      complex<double>(-13, 0), complex<double>(12.0, 0),
      complex<double>(0.0, 0), complex<double>(-5, 0)};

  complex<double> p1_p1 = lorentz_dot(p1, p1);
  complex<double> p2_p2 = lorentz_dot(p2, p2);
  complex<double> p1_p2 = lorentz_dot(p1, p2);
  complex<double> p3_p3 = lorentz_dot(p3, p3);
  complex<double> p1_p3 = lorentz_dot(p1, p3);
  complex<double> p2_p3 = lorentz_dot(p2, p3);

  complex128 p1_p1_f128 = complex128(p1_p1);
  complex128 p2_p2_f128 = complex128(p2_p2);
  complex128 p3_p3_f128 = complex128(p3_p3);
  complex128 p1_p2_f128 = complex128(p1_p2);
  complex128 p1_p3_f128 = complex128(p1_p3);
  complex128 p2_p3_f128 = complex128(p2_p3);

  mppp::complex p1_p1_mpfr = mppp::complex(p1_p1);
  mppp::complex p2_p2_mpfr = mppp::complex(p2_p2);
  mppp::complex p3_p3_mpfr = mppp::complex(p3_p3);
  mppp::complex p1_p2_mpfr = mppp::complex(p1_p2);
  mppp::complex p1_p3_mpfr = mppp::complex(p1_p3);
  mppp::complex p2_p3_mpfr = mppp::complex(p2_p3);

  mppp::complex s_mpfr;
  FFS_mpfr(p1_p1_mpfr, p2_p2_mpfr, p1_p2_mpfr, &s_mpfr);

  complex<double> s;
  complex<double> t;
  complex<double> u;

  complex128 s_f128;
  FFS_f128(p1_p1_f128, p2_p2_f128, p1_p2_f128, &s_f128);

  FFS_f64(p1_p1, p2_p2, p1_p2, &s);
  printf("s: %f\n", s.real());

  FFT_f64(p3_p3, p2_p2, p2_p3, &t);
  printf("t: %f\n", t.real());

  FFU_f64(p1_p1, p3_p3, p1_p3, &u);
  printf("u: %f\n", u.real());

  complex<double> m(1.0, 0.0);

  complex<double> zero(0.0, 0.0);
  double mu = 10.0;

  complex<double> res[3];
  OLO_SCALE(&mu);

  OLO_A0c(&res, &m);
  printf("A0c: %f + i%f\n", res[0].real(), res[0].imag());

  OLO_B0cc(&res, &s, &m, &m);
  printf("B0(s): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_B0cc(&res, &t, &m, &m);
  printf("B0(t): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_B0cc(&res, &u, &m, &m);
  printf("B0(u): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_B0cc(&res, &zero, &m, &m);
  printf("B0(0): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_C0cc(&res, &zero, &zero, &s, &m, &m, &m);
  printf("C0(s): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_C0cc(&res, &zero, &zero, &t, &m, &m, &m);
  printf("C0(t): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_C0cc(&res, &zero, &zero, &u, &m, &m, &m);
  printf("C0(u): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_D0cc(&res, &zero, &zero, &zero, &zero, &s, &t, &m, &m, &m, &m);
  printf("D0(s,t): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_D0cc(&res, &zero, &zero, &zero, &zero, &s, &u, &m, &m, &m, &m);
  printf("D0(s,u): %f + i%f\n", res[0].real(), res[0].imag());

  OLO_D0cc(&res, &zero, &zero, &zero, &zero, &t, &u, &m, &m, &m, &m);
  printf("D0(t,u): %f + i%f\n", res[0].real(), res[0].imag());

  complex<double> out;

  complex<double> astu;
  complex<double> atsu;
  complex<double> aust;
  complex<double> bstu;
  complex<double> btsu;
  complex<double> bust;
  complex<double> cstu;

  APHOAMPFFSTU_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  astu = out;
  printf("astu: %.16e+ i%.16e\n", out.real(), out.imag());

  APHOAMPFFTSU_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  atsu = out;
  printf("atsu: %.16e + i%.16e\n", out.real(), out.imag());

  APHOAMPFFUST_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  aust = out;
  printf("aust: %.16e + i%.16e\n", out.real(), out.imag());

  BPHOAMPFFSTU_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  bstu = out;
  printf("bstu: %.16e + i%.16e\n", out.real(), out.imag());

  BPHOAMPFFTSU_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  btsu = out;
  printf("btsu: %.16e + i%.16e\n", out.real(), out.imag());

  BPHOAMPFFUST_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  bust = out;
  printf("bust: %.16e + i%.16e\n", out.real(), out.imag());

  CPHOAMPFFSTU_f64(p1_p1, p2_p2, p3_p3, p1_p2, p1_p3, p2_p3, &out);
  out = out;
  cstu = out;
  printf("cstu: %.16e + i%.16e\n", out.real(), out.imag());

  complex<double> castu = (astu);
  complex<double> catsu = (atsu);
  complex<double> caust = (aust);
  complex<double> cbstu = (bstu);
  complex<double> cbtsu = (btsu);
  complex<double> cbust = (bust);
  complex<double> ccstu = (cstu);

  complex<double> ans =
      (1 / 8.0) *
      (2.0 * aust * castu + 2.0 * aust * catsu + 4.0 * aust * caust +
       2.0 * aust * cbust * t + 2.0 * caust * bust * t +
       2.0 * cbust * bust * t * t + 2.0 * castu * bstu * u +
       2.0 * catsu * btsu * u - aust * ccstu * s * u +
       2.0 * cbstu * bstu * u * u + 2.0 * cbtsu * btsu * u * u -
       cstu * s * u * (castu + catsu + caust - ccstu * s * u) +
       astu *
           (4.0 * castu - ccstu * s * u + 2.0 * (catsu + caust + cbstu * u)) +
       atsu * (-(ccstu * s * u) +
               2.0 * (castu + 2.0 * catsu + caust + cbtsu * u)));

  printf("ans: %.16e+ i%.16e\n", ans.real(), ans.imag());
  return 0;
}