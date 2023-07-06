#include <tuple>

using namespace std;

template <class T> class LorentzVector {
public:
  T t;
  T x;
  T y;
  T z;

  LorentzVector(T t_new, T x_new, T y_new, T z_new) {
    t = t_new;
    x = x_new;
    y = y_new;
    z = z_new;
  }

  LorentzVector(const T *mom) {
    t = mom[0];
    x = mom[1];
    y = mom[2];
    z = mom[3];
  }

  LorentzVector<T> operator+(LorentzVector<T> const &rhs) {
    LorentzVector<T> res(t + rhs.t, x + rhs.x, y + rhs.y, z + rhs.z);
    return res;
  }

  LorentzVector<T> operator-(LorentzVector<T> const &rhs) {
    LorentzVector<T> res(t - rhs.t, x - rhs.x, y - rhs.y, z - rhs.z);
    return res;
  }

  LorentzVector<T> operator-() {
    LorentzVector<T> res(-t, -x, -y, -z);
    return res;
  }

  LorentzVector<T> operator*(T const &scalar) {
    LorentzVector<T> res(scalar * t, scalar * x, scalar * y, scalar * z);
    return res;
  }

  T dot(LorentzVector<T> const &rhs) {
    return t * rhs.t - x * rhs.x - y * rhs.y - z * rhs.z;
  }
};

template <class T>
tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
c_tensor_common(LorentzVector<T> internal_mom_1,
                LorentzVector<T> internal_mom_2,
                LorentzVector<T> internal_mom_3,
                LorentzVector<T> internal_mom_4) {

  T E1 = internal_mom_1.t;
  T E2 = internal_mom_2.t;
  T E3 = internal_mom_3.t;
  T E4 = internal_mom_4.t;

  T s = 2.0 * internal_mom_1.dot(internal_mom_2);
  T t = 2.0 * internal_mom_2.dot(internal_mom_3);
  T u = 2.0 * internal_mom_1.dot(internal_mom_3);

  // factor 1

  LorentzVector<T> factor_1 = internal_mom_2 * u * E1 -
                              E1 * s * internal_mom_3 -
                              E2 * u * internal_mom_1 + E3 * s * internal_mom_1;

  // factor 2

  LorentzVector<T> factor_2 = -E2 * t * internal_mom_2 +
                              E2 * t * internal_mom_1 -
                              E2 * s * internal_mom_3 + E3 * s * internal_mom_2;

  // factor 3

  LorentzVector<T> factor_3 = -E1 * t * internal_mom_3 +
                              E2 * u * internal_mom_3 +
                              E3 * t * internal_mom_1 - E3 * u * internal_mom_2;

  // factor 4

  LorentzVector<T> factor_4 = -E2 * s * internal_mom_4 +
                              E3 * u * internal_mom_4 +
                              E4 * s * internal_mom_2 - E4 * u * internal_mom_3;

  return make_tuple(factor_1, factor_2, factor_3, factor_4);
}

template <class T>
void amp_f64(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_2,
             LorentzVector<T> external_mom_3, LorentzVector<T> external_mom_4,
             LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
             LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
             T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = external_mom_1.dot(factor1) * external_mom_2.dot(factor2) *
         external_mom_3.dot(factor3) * external_mom_4.dot(factor4);
}

template <class T>
void amp12_f64(LorentzVector<T> external_mom_3, LorentzVector<T> external_mom_4,
               LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
               LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
               T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor1.dot(factor2) * external_mom_3.dot(factor3) *
         external_mom_4.dot(factor4);
}

template <class T>
void amp13_f64(LorentzVector<T> external_mom_2, LorentzVector<T> external_mom_4,
               LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
               LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
               T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor1.dot(factor3) * external_mom_2.dot(factor2) *
         external_mom_4.dot(factor4);
}

template <class T>
void amp23_f64(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_4,
               LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
               LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
               T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor2.dot(factor3) * external_mom_1.dot(factor1) *
         external_mom_4.dot(factor4);
}

template <class T>
void amp14_f64(LorentzVector<T> external_mom_2, LorentzVector<T> external_mom_3,
               LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
               LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
               T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor1.dot(factor4) * external_mom_2.dot(factor2) *
         external_mom_3.dot(factor3);
}

template <class T>
void amp24_f64(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_3,
               LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
               LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
               T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor2.dot(factor4) * external_mom_1.dot(factor1) *
         external_mom_3.dot(factor3);
}

template <class T>
void amp34_f64(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_2,
               LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
               LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
               T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor3.dot(factor4) * external_mom_1.dot(factor1) *
         external_mom_2.dot(factor2);
}

template <class T>
void amp1122_f64(LorentzVector<T> internal_mom_1,
                 LorentzVector<T> internal_mom_2,
                 LorentzVector<T> internal_mom_3,
                 LorentzVector<T> internal_mom_4, T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor1.dot(factor2) * factor3.dot(factor4);
}

template <class T>
void amp1212_f64(LorentzVector<T> internal_mom_1,
                 LorentzVector<T> internal_mom_2,
                 LorentzVector<T> internal_mom_3,
                 LorentzVector<T> internal_mom_4, T *out) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor1.dot(factor3) * factor2.dot(factor4);
}

template <class T>
void amp1221_f64(LorentzVector<T> internal_mom_1,
                 LorentzVector<T> internal_mom_2,
                 LorentzVector<T> internal_mom_3,
                 LorentzVector<T> internal_mom_4, T *out) {
  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  *out = factor1.dot(factor4) * factor2.dot(factor3);
}