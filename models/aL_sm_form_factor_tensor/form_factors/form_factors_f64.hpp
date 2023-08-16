#include "amp.hpp"
#include <complex>

extern void APHOAMPFFSTU_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);
extern void APHOAMPFFTSU_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);
extern void APHOAMPFFUST_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);
extern void BPHOAMPFFSTU_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);
extern void BPHOAMPFFTSU_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);
extern void BPHOAMPFFUST_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);
extern void CPHOAMPFFSTU_f64(complex<double>, complex<double>, complex<double>,
                             complex<double>, complex<double>, complex<double>,
                             complex<double> *);

extern void
amp_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
        LorentzVector<complex<double>>, LorentzVector<complex<double>>,
        LorentzVector<complex<double>>, LorentzVector<complex<double>>,
        LorentzVector<complex<double>>, LorentzVector<complex<double>>,
        complex<double> *);
extern void amp12_f64(LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>, complex<double> *);
extern void amp13_f64(LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>, complex<double> *);
extern void amp23_f64(LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>, complex<double> *);
extern void amp14_f64(LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>, complex<double> *);
extern void amp24_f64(LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>, complex<double> *);
extern void amp34_f64(LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>,
                      LorentzVector<complex<double>>, complex<double> *);
extern void amp1122_f64(LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>, complex<double> *);
extern void amp1221_f64(LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>, complex<double> *);
extern void amp1212_f64(LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>,
                        LorentzVector<complex<double>>, complex<double> *);

class AmpTensor_f64 {
public:
  AmpTensor<complex<double>> amp_tensor_f64;
  AmpTensor_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                LorentzVector<complex<double>>, LorentzVector<complex<double>>);
  void amp_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
               LorentzVector<complex<double>>, LorentzVector<complex<double>>,
               complex<double> *);
  void amp12_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                 complex<double> *);
  void amp13_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                 complex<double> *);
  void amp14_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                 complex<double> *);
  void amp23_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                 complex<double> *);
  void amp24_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                 complex<double> *);
  void amp34_f64(LorentzVector<complex<double>>, LorentzVector<complex<double>>,
                 complex<double> *);
  void amp1122_f64(complex<double> *);
  void amp1212_f64(complex<double> *);
  void amp1221_f64(complex<double> *);
};