#include <mp++/complex128.hpp>
#include "amp.hpp"

typedef mppp::complex128 complex128;

extern void APHOAMPFFSTU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);
extern void APHOAMPFFTSU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);
extern void APHOAMPFFUST_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);
extern void BPHOAMPFFSTU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);
extern void BPHOAMPFFTSU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);
extern void BPHOAMPFFUST_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);
extern void CPHOAMPFFSTU_f128(complex128, complex128, complex128, complex128, complex128, complex128, complex128 *);

extern void amp_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *); 
extern void amp12_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp13_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp23_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp14_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp24_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp34_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp1122_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp1221_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);
extern void amp1212_f128(LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, LorentzVector<complex128>, complex128 *);