#include <complex>
#include <iostream>
#include <stdio.h>
#include <vector>

using namespace std;

#include "form_factors_f64.h"

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
      complex<double>(5, 0), complex<double>(-1.109243, 0),
      complex<double>(-4.448308, 0), complex<double>(1.995529299, 0)};

  complex<double> s = 100.0;
  complex<double> t = 10.0;
  complex<double> u = -s - t;

  complex<double> p1_p2 = s / 2.0;
  complex<double> p1_p3 = u / 2.0;
  complex<double> p2_p3 = t / 2.0;

  cout << p1_p2 << endl;
  cout << p2_p3 << endl;
  cout << p2_p3 << endl;

  complex<double> E1 = p1[0];
  complex<double> E2 = p2[0];
  complex<double> E3 = p3[0];
  complex<double> E4 = -E1 - E2 - E3;

  complex<double> E12 = E1 * E1;
  complex<double> E22 = E2 * E2;
  complex<double> E32 = E3 * E3;
  complex<double> E42 = E4 * E4;

  complex<double> astu;
  complex<double> atsu;
  complex<double> aust;
  complex<double> bstu;
  complex<double> btsu;
  complex<double> bust;
  complex<double> cstu;

  APHOAMPFFSTU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &astu);
  astu = astu * E1 * E2 * E3 * E42 * s * s;
  APHOAMPFFTSU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &atsu);
  atsu = atsu * E1 * E2 * E3 * E42 * t * t;
  APHOAMPFFUST_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &aust);
  aust = aust * E1 * E2 * E3 * E42 * u * u;
  BPHOAMPFFSTU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &bstu);
  bstu = bstu * E1 * E2 * E32 * E42 * s * s * t * t * u;
  BPHOAMPFFTSU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &btsu);
  btsu = btsu * E12 * E2 * E3 * E42 * s * s * t * t * u;
  BPHOAMPFFUST_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &bust);
  bust = bust * E1 * E22 * E3 * E42 * s * s * t * u * u;
  CPHOAMPFFSTU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &cstu);
  cstu = cstu * E1 * E2 * E3 * E4 * s * t * t * u;
  cout << "this thing: " << (E1 * E2 * E3 * E4 * s * t * t * u) << endl;

  printf("astu: %.16e + i %.16e\n", astu.real(), astu.imag());
  printf("atsu: %.16e + i %.16e\n", atsu.real(), atsu.imag());
  printf("aust: %.16e + i %.16e\n", aust.real(), aust.imag());
  printf("bstu: %.16e + i %.16e\n", bstu.real(), bstu.imag());
  printf("btsu: %.16e + i %.16e\n", btsu.real(), btsu.imag());
  printf("bust: %.16e + i %.16e\n", bust.real(), bust.imag());
  printf("cstu: %.16e + i %.16e\n", cstu.real(), cstu.imag());

  return 0;
}