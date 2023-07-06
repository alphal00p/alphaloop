#ifndef AMP
#define AMP

#include <iostream>
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

  void debug() {
    cout << "( t = " << t << ", x = " << x << ", y = " << y << ", z = " << z
         << ")" << endl;
  }
};

template <class T>
LorentzVector<T> operator*(T const &scalar, LorentzVector<T> const &vector) {
  LorentzVector<T> res(scalar * vector.t, scalar * vector.x, scalar * vector.y,
                       scalar * vector.z);
  return res;
}

template <class T>
tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
c_tensor_common(LorentzVector<T> internal_mom_1,
                LorentzVector<T> internal_mom_2,
                LorentzVector<T> internal_mom_3,
                LorentzVector<T> internal_mom_4, int debug = 0) {

  T E1 = internal_mom_1.t;
  T E2 = internal_mom_2.t;
  T E3 = internal_mom_3.t;
  T E4 = internal_mom_4.t;

  string red = "\033[31m";
  string green = "\033[32m";
  string yellow = "\033[33m";
  string blue = "\033[34m";
  string white = "\033[37m";
  {}

  if (debug > 0) {
    cout << blue << "start evaluation c_tensor_common ------------------"
         << endl;
    cout << white << "(E1 = " << E1;
    cout << ", E2 = " << E2;
    cout << ", E3 = " << E3;
    cout << ", E4 = " << E4;
    cout << ")" << endl;
  }
  T s = 2.0 * internal_mom_1.dot(internal_mom_2);
  T t = 2.0 * internal_mom_2.dot(internal_mom_3);
  T u = 2.0 * internal_mom_1.dot(internal_mom_3);

  if (debug > 0) {
    cout << "(s = " << s;
    cout << ", t = " << t;
    cout << ", u = " << u;
    cout << ")" << endl;
  }

  // factor 1

  LorentzVector<T> factor_1 = internal_mom_2 * u * E1 -
                              E1 * s * internal_mom_3 -
                              E2 * u * internal_mom_1 + E3 * s * internal_mom_1;
  if (debug > 0) {
    if (debug > 1) {
      cout << yellow << "evaluation of factor 1 -------" << endl;
      cout << white << "term1 = ";
      (internal_mom_2 * u * E1).debug();
      cout << "term2 = ";
      (-E1 * s * internal_mom_3).debug();
      cout << "term3 = ";
      (-E2 * u * internal_mom_1).debug();
      cout << "term4 = ";
      (E3 * s * internal_mom_1).debug();
    }
    cout << "factor_1 = ";
    factor_1.debug();
  }

  // factor 2

  LorentzVector<T> factor_2 = -E1 * t * internal_mom_2 +
                              E2 * t * internal_mom_1 -
                              E2 * s * internal_mom_3 + E3 * s * internal_mom_2;

  if (debug > 0) {
    if (debug > 1) {
      cout << yellow << "evaluation of factor 2 -------" << endl;
      cout << white << "term1 = ";
      (-E2 * t * internal_mom_2).debug();
      cout << "term2 = ";
      (E2 * t * internal_mom_1).debug();
      cout << "term3 = ";
      (-E2 * s * internal_mom_3).debug();
      cout << "term4 = ";
      (E3 * s * internal_mom_2).debug();
    }
    cout << "factor_2 = ";
    factor_2.debug();
  }

  // factor 3

  LorentzVector<T> factor_3 = -E1 * t * internal_mom_3 +
                              E2 * u * internal_mom_3 +
                              E3 * t * internal_mom_1 - E3 * u * internal_mom_2;

  if (debug > 0) {
    if (debug > 1) {
      cout << yellow << "evaluation of factor 3 -------" << endl;
      cout << white << "term1 = ";
      (-E1 * t * internal_mom_3).debug();
      cout << "term2 = ";
    }
    cout << "factor_3 = ";
    factor_3.debug();
  }

  // factor 4

  LorentzVector<T> factor_4 = -E2 * s * internal_mom_4 +
                              E3 * u * internal_mom_4 +
                              E4 * s * internal_mom_2 - E4 * u * internal_mom_3;
  if (debug > 0) {
    if (debug > 1) {
      cout << yellow << "evaluation of factor 4 -------" << endl;
      cout << white << "term1 = ";
      (-E2 * s * internal_mom_4).debug();
      cout << "term2 = ";
      (E3 * u * internal_mom_4).debug();
      cout << "term3 = ";
      (E4 * s * internal_mom_2).debug();
      cout << "term4 = ";
      (-E4 * u * internal_mom_3).debug();
    }
    cout << "factor_4 = ";
    factor_4.debug();
  }

  if (debug > 0) {

    cout << "factor_4 = ";
    factor_4.debug();
  }

  return make_tuple(factor_1, factor_2, factor_3, factor_4);
}

template <class T>
T amp(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_2,
      LorentzVector<T> external_mom_3, LorentzVector<T> external_mom_4,
      LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
      LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4,
      int debug = 0) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4, debug);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  if (debug > 0) {
    cout << "factor1 = " << external_mom_1.dot(factor1) << endl;
    cout << "factor2 = " << external_mom_2.dot(factor2) << endl;
    cout << "factor3 = " << external_mom_3.dot(factor3) << endl;
    cout << "factor4 = " << external_mom_4.dot(factor4) << endl;
  }

  return external_mom_1.dot(factor1) * external_mom_2.dot(factor2) *
         external_mom_3.dot(factor3) * external_mom_4.dot(factor4);
}

template <class T>
T amp12(LorentzVector<T> external_mom_3, LorentzVector<T> external_mom_4,
        LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
        LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor1.dot(factor2) * external_mom_3.dot(factor3) *
         external_mom_4.dot(factor4);
}

template <class T>
T amp13(LorentzVector<T> external_mom_2, LorentzVector<T> external_mom_4,
        LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
        LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor1.dot(factor3) * external_mom_2.dot(factor2) *
         external_mom_4.dot(factor4);
}

template <class T>
T amp23(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_4,
        LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
        LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor2.dot(factor3) * external_mom_1.dot(factor1) *
         external_mom_4.dot(factor4);
}

template <class T>
T amp14(LorentzVector<T> external_mom_2, LorentzVector<T> external_mom_3,
        LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
        LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor1.dot(factor4) * external_mom_2.dot(factor2) *
         external_mom_3.dot(factor3);
}

template <class T>
T amp24(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_3,
        LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
        LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor2.dot(factor4) * external_mom_1.dot(factor1) *
         external_mom_3.dot(factor3);
}

template <class T>
T amp34(LorentzVector<T> external_mom_1, LorentzVector<T> external_mom_2,
        LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
        LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor3.dot(factor4) * external_mom_1.dot(factor1) *
         external_mom_2.dot(factor2);
}

template <class T>
T amp1122(LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
          LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor1.dot(factor2) * factor3.dot(factor4);
}

template <class T>
T amp1212(LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
          LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {

  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor1.dot(factor3) * factor2.dot(factor4);
}

template <class T>
T amp1221(LorentzVector<T> internal_mom_1, LorentzVector<T> internal_mom_2,
          LorentzVector<T> internal_mom_3, LorentzVector<T> internal_mom_4) {
  tuple<LorentzVector<T>, LorentzVector<T>, LorentzVector<T>, LorentzVector<T>>
      factors = c_tensor_common(internal_mom_1, internal_mom_2, internal_mom_3,
                                internal_mom_4);

  LorentzVector<T> factor1 = get<0>(factors);
  LorentzVector<T> factor2 = get<1>(factors);
  LorentzVector<T> factor3 = get<2>(factors);
  LorentzVector<T> factor4 = get<3>(factors);

  return factor1.dot(factor4) * factor2.dot(factor3);
}

#endif