#include "amp.hpp"
#include <complex>
#include <iostream>

using namespace std;

int main() {
  LorentzVector<complex<double>> p1(13.0, 0.0, 0.0, -13.0);
  LorentzVector<complex<double>> p2(13.0, 0.0, 0.0, 13.0);

  LorentzVector<complex<double>> p3(-13, 12.0, 0.0, -5.0);
  LorentzVector<complex<double>> p4(-13, -12.0, 0, 5.0);

  LorentzVector<complex<double>> n1(1.0, 2.0, 3.0, 4.0);
  LorentzVector<complex<double>> n2(2.0, 4.0, 8.0, 16.0);

  LorentzVector<complex<double>> n3(3.0, 2.0, 1.0, 0.0);
  LorentzVector<complex<double>> n4(9.0, 32.0, 2.0, 1.0);

  cout << "k1p2" << n1.dot(p2) << endl;
  cout << "k1p3" << n1.dot(p3) << endl;
  cout << "k1p1" << n1.dot(p1) << endl;
  cout << "k2p2" << n2.dot(p2) << endl;
  cout << "k2p1" << n2.dot(p1) << endl;
  cout << "k2p3" << n2.dot(p3) << endl;
  cout << "k3p3" << n3.dot(p3) << endl;
  cout << "k3p1" << n3.dot(p1) << endl;
  cout << "k3p2" << n3.dot(p2) << endl;
  cout << "k4p4" << n4.dot(p4) << endl;
  cout << "k4p2" << n4.dot(p2) << endl;
  cout << "k4p3" << n4.dot(p3) << endl;

  cout << "s: " << 2.0 * p1.dot(p2) << endl;
  cout << "t: " << 2.0 * p2.dot(p3) << endl;
  cout << "u: " << 2.0 * p1.dot(p3) << endl;

  complex<double> normal = amp(n1, n2, n3, n4, p1, p2, p3, p4, 1);

  cout << "normal " << normal << endl;

  complex<double> ward1 = amp(p1, n2, n3, n4, p1, p2, p3, p4);
  complex<double> ward2 = amp(n1, p2, n3, n4, p1, p2, p3, p4);
  complex<double> ward3 = amp(n1, n2, p3, n4, p1, p2, p3, p4);
  complex<double> ward4 = amp(n1, n2, n3, p4, p1, p2, p3, p4);

  cout << "ward1 " << ward1 << endl;
  cout << "ward2 " << ward2 << endl;
  cout << "ward3 " << ward3 << endl;
  cout << "ward4 " << ward4 << endl;

  return 0;
}