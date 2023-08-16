#ifndef TENSOR
#define TENSOR

#include <array>
#include <bits/types/struct_sched_param.h>
#include <complex>
#include <iostream>
#include <mp++/complex128.hpp>

using namespace std;

template <class F> using array2d = array<array<F, 3>, 3>;
template <class F> using array4d = array<array<array<array<F, 3>, 3>, 3>, 3>;

typedef mppp::complex128 complex128;

template <class F> class LorentzVector {
public:
  F t;
  F x;
  F y;
  F z;

  LorentzVector() {
    t = 0.;
    x = 0.;
    y = 0.;
    z = 0.;
  }

  LorentzVector(F t_new, F x_new, F y_new, F z_new) {
    t = t_new;
    x = x_new;
    y = y_new;
    z = z_new;
  }

  LorentzVector(const F *mom) {
    t = mom[0];
    x = mom[1];
    y = mom[2];
    z = mom[3];
  }

  LorentzVector<F> operator+(LorentzVector<F> const &rhs) {
    LorentzVector<F> res(t + rhs.t, x + rhs.x, y + rhs.y, z + rhs.z);
    return res;
  }

  LorentzVector<F> operator-(LorentzVector<F> const &rhs) {
    LorentzVector<F> res(t - rhs.t, x - rhs.x, y - rhs.y, z - rhs.z);
    return res;
  }

  LorentzVector<F> operator-() {
    LorentzVector<F> res(-t, -x, -y, -z);
    return res;
  }

  LorentzVector<F> operator*(F const &scalar) {
    LorentzVector<F> res(scalar * t, scalar * x, scalar * y, scalar * z);
    return res;
  }

  F dot(LorentzVector<F> const &rhs) {
    return t * rhs.t - x * rhs.x - y * rhs.y - z * rhs.z;
  }

  F get(int mu) {
    switch (mu) {
    case 0:
      return t;
    case 1:
      return x;
    case 2:
      return y;
    case 3:
      return z;
    }
  }

  void set(int mu, F val) {
    switch (mu) {
    case 0:
      t = val;
    case 1:
      x = val;
    case 2:
      y = val;
    case 3:
      z = val;
    }
  }
};

template <class F>
LorentzVector<F> operator*(F const &scalar, LorentzVector<F> const &vector) {
  LorentzVector<F> res(scalar * vector.t, scalar * vector.x, scalar * vector.y,
                       scalar * vector.z);
  return res;
}

template <class F> class SpatialVector {
public:
  array<F, 3> entries;

  SpatialVector() {}

  SpatialVector(array<F, 3> new_entries) { entries = new_entries; }
  SpatialVector(LorentzVector<F> mom) {
    entries[0] = mom.x;
    entries[1] = mom.y;
    entries[2] = mom.z;
  }

  SpatialVector(F x, F y, F z) {
    entries[0] = x;
    entries[1] = y;
    entries[2] = z;
  }

  F dot(SpatialVector<F> rhs) {
    F res = entries[0] * rhs.entries[0] + entries[1] * rhs.entries[1] +
            entries[2] * rhs.entries[2];
    return res;
  }

  F energy() {
    return sqrt(entries[0] * entries[0] + entries[1] * entries[1] +
                entries[2] * entries[2]);
  }

  F operator[](int i) { return entries[i]; }

  SpatialVector<F> operator+(SpatialVector<F> const &rhs) {
    array<F, 3> new_entries;
    for (int i = 0; i < 3; i++) {
      new_entries[i] = entries[i] + rhs.entries[i];
    }
    SpatialVector<F> res(new_entries);
    return res;
  }

  SpatialVector<F> operator-(SpatialVector<F> const &rhs) {
    array<F, 3> new_entries;
    for (int i = 0; i < 3; i++) {
      new_entries[i] = entries[i] - rhs.entries[i];
    }
    SpatialVector<F> res(new_entries);
    return res;
  }
};

template <class F>
SpatialVector<F> operator*(F const &lhs, SpatialVector<F> const &rhs) {
  array<F, 3> new_entries;

  for (int i = 0; i < 3; i++) {
    new_entries[i] = lhs * rhs.entries[i];
  }

  SpatialVector<F> res(new_entries);
  return res;
}

template <class F> class Tensor2 {
public:
  array<array<F, 3>, 3> entries;
  Tensor2(array<array<F, 3>, 3> new_entries) { entries = new_entries; }

  Tensor2(SpatialVector<F> lhs, SpatialVector<F> rhs) {

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        entries[i][j] = lhs.entries[i] * rhs.entries[j];
      }
    }
  }

  void debug() {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        cout << entries[i][j];
      }
      cout << endl;
    }
  }

  Tensor2<F> spin_sum_project(Tensor2<F> spin_sum_left,
                              Tensor2<F> spin_sum_right) {
    array<array<F, 3>, 3> new_entries;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        new_entries[i][j] = 0.;
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            new_entries[i][j] += spin_sum_left.entries[i][k] * entries[k][l] *
                                 spin_sum_right.entries[j][l];
          }
        }
      }
    }

    Tensor2<F> res(new_entries);
    return res;
  }

  void add_kronecker(F scalar) {
    for (int i = 0; i < 3; i++) {
      entries[i][i] += scalar;
    }
  }
};

template <class F> Tensor2<F> spin_sum(SpatialVector<F> mom) {
  SpatialVector<F> rescaled_mom = (1. / mom.energy()) * mom;
  array<array<F, 3>, 3> new_entries;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      new_entries[i][j] = -rescaled_mom.entries[i] * rescaled_mom.entries[j];
    }
  }

  for (int i = 0; i < 3; i++) {
    new_entries[i][i] += 1;
  }

  Tensor2<F> res(new_entries);
  return res;
}

template <class F> class Tensor4 {
public:
  array<array<array<array<F, 3>, 3>, 3>, 3> entries;
  Tensor4() {}
  Tensor4(array<array<array<array<F, 3>, 3>, 3>, 3> new_entries) {
    entries = new_entries;
  }

  F full_contract(SpatialVector<F> mom_1, SpatialVector<F> mom_2,
                  SpatialVector<F> mom_3, SpatialVector<F> mom_4) {
    F res = 0.;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            res += mom_1.entries[i] * mom_2.entries[j] * mom_3.entries[k] *
                   mom_4.entries[l] * entries[i][j][k][l];
          }
        }
      }
    }

    return res;
  }
};

template <class F> inline F delta(int i, int j) {
  if (i == j) {
    return 1.;
  } else {
    return 0.;
  }
}

template <class F> class AmpTensor {
public:
  array4d<F> total_entries;

  array2d<F> entries12;
  array2d<F> entries13;
  array2d<F> entries14;
  array2d<F> entries23;
  array2d<F> entries24;
  array2d<F> entries34;

  F entry1122;
  F entry1212;
  F entry1221;

  SpatialVector<F> unit_mom_1;
  SpatialVector<F> unit_mom_2;
  SpatialVector<F> unit_mom_3;
  SpatialVector<F> unit_mom_4;

  AmpTensor(LorentzVector<F> l_mom_1, LorentzVector<F> l_mom_2,
            LorentzVector<F> l_mom_3, LorentzVector<F> l_mom_4, F f_astu,
            F f_atsu, F f_aust, F f_bstu, F f_btsu, F f_bust, F f_cstu) {

    array2d<F> spin_sum_1_entries;
    array2d<F> spin_sum_2_entries;
    array2d<F> spin_sum_3_entries;
    array2d<F> spin_sum_4_entries;

    F E1 = l_mom_1.t;
    F E2 = l_mom_2.t;
    F E3 = l_mom_3.t;
    F E4 = l_mom_4.t;

    SpatialVector<F> mom_1(l_mom_1);
    SpatialVector<F> mom_2(l_mom_2);
    SpatialVector<F> mom_3(l_mom_3);
    SpatialVector<F> mom_4(l_mom_4);

    F s_half = E1 * E2 - mom_1.dot(mom_2);
    F t_half = E2 * E3 - mom_2.dot(mom_3);
    F u_half = E3 * E1 - mom_3.dot(mom_1);

    unit_mom_1 = (1. / E1) * mom_1;
    unit_mom_2 = (1. / E2) * mom_2;
    unit_mom_3 = (1. / E3) * mom_3;
    unit_mom_4 = (1. / E4) * mom_4;

    for (int i1 = 0; i1 < 3; i1++) {
      for (int j1 = 0; j1 < 3; j1++) {
        spin_sum_1_entries[i1][j1] =
            delta<F>(i1, j1) - unit_mom_1[i1] * unit_mom_1[j1];
        spin_sum_2_entries[i1][j1] =
            delta<F>(i1, j1) - unit_mom_2[i1] * unit_mom_2[j1];
        spin_sum_3_entries[i1][j1] =
            delta<F>(i1, j1) - unit_mom_3[i1] * unit_mom_3[j1];
        spin_sum_4_entries[i1][j1] =
            delta<F>(i1, j1) - unit_mom_4[i1] * unit_mom_4[j1];
      }
    }

    array4d<F> total_entries_no_spin_sum;

    for (int i1 = 0; i1 < 3; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i4 = 0; i4 < 3; i4++) {

            // astu
            total_entries_no_spin_sum[i1][i2][i3][i4] =
                f_astu *
                (-delta<F>(i1, i2) - (1. / s_half) * mom_2[i1] * mom_1[i2]) *
                (-delta<F>(i3, i4) +
                 (1. / s_half) * (mom_1[i3] + mom_2[i3]) * mom_3[i4]);

            // //   // atsu
            total_entries_no_spin_sum[i1][i2][i3][i4] +=
                f_atsu *
                (-delta<F>(i3, i2) - (1. / t_half) * mom_2[i3] * mom_3[i2]) *
                (-delta<F>(i1, i4) +
                 (1. / t_half) * (mom_3[i1] + mom_2[i1]) * mom_1[i4]);

            // aust
            total_entries_no_spin_sum[i1][i2][i3][i4] +=
                f_aust *
                (-delta<F>(i3, i1) - (1. / u_half) * mom_1[i3] * mom_3[i1]) *
                (-delta<F>(i2, i4) +
                 (1. / u_half) * (mom_3[i2] + mom_1[i2]) * mom_2[i4]);

            // bstu
            total_entries_no_spin_sum[i1][i2][i3][i4] +=
                f_bstu *
                ((delta<F>(i1, i2) + (1. / s_half) * mom_2[i1] * mom_1[i2]) *
                     (mom_1[i3] - u_half / t_half * mom_2[i3]) *
                     (mom_2[i4] - u_half / s_half * mom_3[i4]) +
                 (mom_3[i1] - u_half / s_half * mom_2[i1]) *
                     (mom_1[i2] - s_half / t_half * mom_3[i2]) *
                     (-delta<F>(i3, i4) + 1. / t_half * mom_2[i3] * mom_1[i4] +
                      1. / u_half * mom_1[i3] * mom_2[i4]));

            // btsu
            total_entries_no_spin_sum[i1][i2][i3][i4] +=
                f_btsu *
                ((delta<F>(i3, i2) + (1. / t_half) * mom_3[i2] * mom_2[i3]) *
                     (mom_3[i1] - u_half / s_half * mom_2[i1]) *
                     (mom_2[i4] - u_half / t_half * mom_1[i4]) +
                 (mom_1[i3] - u_half / t_half * mom_2[i3]) *
                     (mom_3[i2] - t_half / s_half * mom_1[i2]) *
                     (-delta<F>(i1, i4) +
                      (1. / t_half) * mom_2[i1] * mom_3[i4] +
                      1. / u_half * mom_3[i1] * mom_2[i4]));

            // bust
            total_entries_no_spin_sum[i1][i2][i3][i4] +=
                f_bust *
                ((delta<F>(i3, i1) + 1. / u_half * mom_3[i1] * mom_1[i3]) *
                     (mom_3[i2] - t_half / s_half * mom_1[i2]) *
                     (mom_1[i4] - t_half / u_half * mom_2[i4]) +
                 (mom_2[i3] - t_half / u_half * mom_1[i3]) *
                     (mom_3[i1] - u_half / s_half * mom_2[i1]) *
                     (-delta<F>(i2, i4) + 1. / s_half * mom_1[i2] * mom_3[i4] +
                      1. / t_half * mom_3[i2] * mom_1[i4]));

            // cstu
            total_entries_no_spin_sum[i1][i2][i3][i4] +=
                f_cstu / (s_half * t_half * t_half * u_half) *
                (s_half * mom_3[i1] - u_half * mom_2[i1]) *
                (s_half * mom_3[i2] - t_half * mom_1[i2]) *
                (u_half * mom_2[i3] - t_half * mom_1[i3]) *
                (s_half * mom_2[i4] - u_half * mom_3[i4]);
          }
        }
      }
    }

    for (int i1 = 0; i1 < 3; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i4 = 0; i4 < 3; i4++) {
            total_entries[i1][i2][i3][i4] = 0.;

            for (int j1 = 0; j1 < 3; j1++) {
              for (int j2 = 0; j2 < 3; j2++) {
                for (int j3 = 0; j3 < 3; j3++) {
                  for (int j4 = 0; j4 < 3; j4++) {
                    total_entries[i1][i2][i3][i4] +=
                        spin_sum_1_entries[i1][j1] *
                        spin_sum_2_entries[i2][j2] *
                        spin_sum_3_entries[i3][j3] *
                        spin_sum_4_entries[i4][j4] *
                        total_entries_no_spin_sum[j1][j2][j3][j4];
                  }
                }
              }
            }
          }
        }
      }
    }

    entry1122 = 0.;
    entry1212 = 0.;
    entry1221 = 0.;

    for (int i1 = 0; i1 < 3; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {

        entries12[i1][i2] = 0;
        entries13[i1][i2] = 0;
        entries14[i1][i2] = 0;
        entries23[i1][i2] = 0;
        entries24[i1][i2] = 0;
        entries34[i1][i2] = 0;

        for (int i3 = 0; i3 < 3; i3++) {
          entries12[i1][i2] += -total_entries[i3][i3][i1][i2];
          entries13[i1][i2] += -total_entries[i3][i1][i3][i2];
          entries14[i1][i2] += -total_entries[i3][i1][i2][i3];
          entries23[i1][i2] += -total_entries[i1][i3][i3][i2];
          entries24[i1][i2] += -total_entries[i1][i3][i2][i3];
          entries34[i1][i2] += -total_entries[i1][i2][i3][i3];
        }

        entry1122 += total_entries[i1][i1][i2][i2];
        entry1212 += total_entries[i1][i2][i1][i2];
        entry1221 += total_entries[i1][i2][i2][i1];
      }
    }
  }

  F amp(SpatialVector<F> ex_mom_1, SpatialVector<F> ex_mom_2,
        SpatialVector<F> ex_mom_3, SpatialVector<F> ex_mom_4) {
    F res = 0;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            res += ex_mom_1.entries[i] * ex_mom_2.entries[j] *
                   ex_mom_3.entries[k] * ex_mom_4.entries[l] *
                   total_entries[i][j][k][l];
          }
        }
      }
    }

    return res;
  }

  F amp12(SpatialVector<F> ex_mom_3, SpatialVector<F> ex_mom_4) {
    F res = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res += ex_mom_3[i] * ex_mom_4[j] * entries12[i][j];
      }
    }

    return res;
  }

  F amp13(SpatialVector<F> ex_mom_2, SpatialVector<F> ex_mom_4) {
    F res = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res += ex_mom_2[i] * ex_mom_4[j] * entries13[i][j];
      }
    }

    return res;
  }

  F amp14(SpatialVector<F> ex_mom_2, SpatialVector<F> ex_mom_3) {
    F res = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res += ex_mom_2[i] * ex_mom_3[j] * entries14[i][j];
      }
    }

    return res;
  }

  F amp23(SpatialVector<F> ex_mom_1, SpatialVector<F> ex_mom_4) {
    F res = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res += ex_mom_1[i] * ex_mom_4[j] * entries23[i][j];
      }
    }

    return res;
  }

  F amp24(SpatialVector<F> ex_mom_1, SpatialVector<F> ex_mom_3) {
    F res = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res += ex_mom_1[i] * ex_mom_3[j] * entries24[i][j];
      }
    }

    return res;
  }

  F amp34(SpatialVector<F> ex_mom_1, SpatialVector<F> ex_mom_2) {
    F res = 0.;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        res += ex_mom_1[i] * ex_mom_2[j] * entries34[i][j];
      }
    }

    return res;
  }
  AmpTensor() {}
  F amp1122() { return entry1122; }
  F amp1212() { return entry1212; }
  F amp1221() { return entry1221; }
};

#endif