#include <complex>
#include <iostream>
#include <vector>

using namespace std;

#include "form_factor_benchmark.h"

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

  complex<double> p1_p2 = lorentz_dot(p1, p2);
  complex<double> p1_p3 = lorentz_dot(p1, p3);
  complex<double> p2_p3 = lorentz_dot(p2, p3);

  cout << p1_p2 << endl;
  cout << p2_p3 << endl;
  cout << p2_p3 << endl;

  complex<double> E1 = p1[0];
  complex<double> E2 = p2[0];
  complex<double> E3 = p3[0];

  complex<double> astu;
  complex<double> atsu;
  complex<double> aust;
  complex<double> bstu;
  complex<double> btsu;
  complex<double> bust;
  complex<double> cstu;

  APHOAMPFFSTU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &astu);
  APHOAMPFFTSU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &atsu);
  APHOAMPFFUST_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &aust);
  BPHOAMPFFSTU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &bstu);
  BPHOAMPFFTSU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &btsu);
  BPHOAMPFFUST_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &bust);
  CPHOAMPFFSTU_f64(E1, E2, E3, p1_p2, p1_p3, p2_p3, &cstu);

  printf("astu: %.16e \n", astu.real());
  printf("atsu: %.16e \n", atsu.real());
  printf("aust: %.16e \n", aust.real());
  printf("bstu: %.16e \n", bstu.real());
  printf("btsu: %.16e \n", btsu.real());
  printf("bust: %.16e \n", bust.real());
  printf("cstu: %.16e \n", cstu.real());

  return 0;
}