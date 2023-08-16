#include "form_factors_f128.hpp"
#include "amp.hpp"
#include "form_factors_f64.hpp"
#include <complex>
#include <mp++/complex128.hpp>
#include <mp++/integer.hpp>

typedef mppp::complex128 complex128;

using namespace std;
void APHOAMPFFSTU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  APHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}

void APHOAMPFFTSU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  APHOAMPFFTSU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void APHOAMPFFUST_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  APHOAMPFFUST_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void BPHOAMPFFSTU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  BPHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void BPHOAMPFFTSU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  BPHOAMPFFTSU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void BPHOAMPFFUST_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  BPHOAMPFFUST_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}
void CPHOAMPFFSTU_f128(complex128 E1, complex128 E2, complex128 E3,
                       complex128 p1_p2, complex128 p1_p3, complex128 p2_p3,
                       complex128 *out) {
  complex<double> out_f64;
  CPHOAMPFFSTU_f64((complex<double>)E1, (complex<double>)E2,
                   (complex<double>)E3, (complex<double>)p1_p2,
                   (complex<double>)p1_p3, (complex<double>)p2_p3, &out_f64);
  *out = (complex128)out_f64;
}

AmpTensor_f128::AmpTensor_f128(LorentzVector<complex128> mom_1,
                               LorentzVector<complex128> mom_2,
                               LorentzVector<complex128> mom_3,
                               LorentzVector<complex128> mom_4) {
  complex128 f_astu;
  complex128 f_atsu;
  complex128 f_aust;
  complex128 f_bstu;
  complex128 f_btsu;
  complex128 f_bust;
  complex128 f_cstu;

  complex128 p1_p2 = mom_1.dot(mom_2);
  complex128 p2_p3 = mom_2.dot(mom_3);
  complex128 p1_p3 = mom_3.dot(mom_1);

  APHOAMPFFSTU_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_astu);
  APHOAMPFFTSU_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_atsu);
  APHOAMPFFUST_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_aust);
  BPHOAMPFFSTU_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_bstu);
  BPHOAMPFFTSU_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_btsu);
  BPHOAMPFFUST_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_bust);
  CPHOAMPFFSTU_f128(mom_1.t, mom_2.t, mom_3.t, p1_p2, p1_p3, p2_p3, &f_cstu);
  amp_tensor_f128 =
      AmpTensor<complex128>(mom_1, mom_2, mom_3, mom_4, f_astu, f_atsu, f_aust,
                            f_bstu, f_btsu, f_bust, f_cstu);
}

void AmpTensor_f128::amp_f128(LorentzVector<complex128> mom_1,
                              LorentzVector<complex128> mom_2,
                              LorentzVector<complex128> mom_3,
                              LorentzVector<complex128> mom_4,
                              complex128 *out) {
  *out = amp_tensor_f128.amp(
      SpatialVector<complex128>(mom_1), SpatialVector<complex128>(mom_2),
      SpatialVector<complex128>(mom_3), SpatialVector<complex128>(mom_4));
}

void AmpTensor_f128::amp12_f128(LorentzVector<complex128> mom_3,
                                LorentzVector<complex128> mom_4,
                                complex128 *out) {
  *out = amp_tensor_f128.amp12(SpatialVector<complex128>(mom_3),
                               SpatialVector<complex128>(mom_4));
}

void AmpTensor_f128::amp13_f128(LorentzVector<complex128> mom_2,
                                LorentzVector<complex128> mom_4,
                                complex128 *out) {
  *out = amp_tensor_f128.amp13(SpatialVector<complex128>(mom_2),
                               SpatialVector<complex128>(mom_4));
}

void AmpTensor_f128::amp14_f128(LorentzVector<complex128> mom_2,
                                LorentzVector<complex128> mom_3,
                                complex128 *out) {
  *out = amp_tensor_f128.amp14(SpatialVector<complex128>(mom_2),
                               SpatialVector<complex128>(mom_3));
}

void AmpTensor_f128::amp23_f128(LorentzVector<complex128> mom_1,
                                LorentzVector<complex128> mom_4,
                                complex128 *out) {
  *out = amp_tensor_f128.amp23(SpatialVector<complex128>(mom_1),
                               SpatialVector<complex128>(mom_4));
}

void AmpTensor_f128::amp24_f128(LorentzVector<complex128> mom_1,
                                LorentzVector<complex128> mom_3,
                                complex128 *out) {
  *out = amp_tensor_f128.amp23(SpatialVector<complex128>(mom_1),
                               SpatialVector<complex128>(mom_3));
}

void AmpTensor_f128::amp34_f128(LorentzVector<complex128> mom_1,
                                LorentzVector<complex128> mom_2,
                                complex128 *out) {
  *out = amp_tensor_f128.amp34(SpatialVector<complex128>(mom_1),
                               SpatialVector<complex128>(mom_2));
}

void AmpTensor_f128::amp1122_f128(complex128 *out) {
  *out = amp_tensor_f128.entry1122;
}

void AmpTensor_f128::amp1212_f128(complex128 *out) {
  *out = amp_tensor_f128.entry1212;
}

void AmpTensor_f128::amp1221_f128(complex128 *out) {
  *out = amp_tensor_f128.entry1221;
}
