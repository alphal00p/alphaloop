//===-- duals/dualt2 - Dual number class --------------------------*- C++ -*-===//
//
// Extension of cppduals project.
// https://tesch1.gitlab.io/cppduals
//
// (c)2019 Michael Tesch. tesch1@gmail.com
//
// See https://gitlab.com/tesch1/cppdualt2s/blob/master/LICENSE.txt for
// license information.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DUALT2
#define DUALT2

#ifndef PARSED_BY_DOXYGEN
#include <cmath>
#include <ctgmath>
#include <limits>
#include <type_traits>
#include <complex>
#include <random>
#include <iostream>
#include <cassert>
#endif

#if !defined(CPPDUALST2_IGNORE_COMPILER_VERSION) && !defined(_WIN32)
#if __cplusplus < 201103L
  #error CPPDUALST2 needs at least a C++11 compliant compiler
#endif
#endif

namespace dualst2 {

#ifndef PARSED_BY_DOXYGEN
template<class T> class dualt2;
#endif

/// Check if T is dualt2<>, match non-dualt2s.
template <class T> struct is_dualt2 : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is dualt2<>, match dualt2<>.
template <class T> struct is_dualt2<dualt2<T> > : std::true_type {};

#endif

/// Check if T is std::complex<>, match non- std::complex<>.
template <class T> struct is_complex : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is std::complex<>, match std::complex<>.
template <class T> struct is_complex<std::complex<T> > : std::true_type {};

#endif

/// dualt2_traits helper class.
template <class T> struct dualt2_traits
{
  /// Depth of T - for T=scalar this is 0. for dualt2_traits<double> it
  /// is 1.
  enum { depth = 0 }; // -Wenum-compare

  /// The real storage type.
  typedef T real_type;
};

#ifndef PARSED_BY_DOXYGEN

/// dualt2_traits for dualt2<> types
template <class T> struct dualt2_traits<dualt2<T>>
{
  /// Depth to which this dualt2<> type is nested.  One (1) is a
  /// first-level dualt2, whereas non-dualt2s have a depth of 0.
  enum { depth = dualt2_traits<T>::depth + 1 };

  /// The real storage type.
  typedef typename dualt2_traits<T>::real_type real_type;
};

template <class T> struct dualt2_traits<std::complex<dualt2<T>>>
{
  /// complex<dualt2<T>> have the same 'depth' as their dual.
  enum { depth = dualt2_traits<T>::depth };

  /// The real storage type.
  typedef typename dualt2_traits<T>::real_type real_type;
};

namespace detail {

template<class T>
struct Void { typedef void type; };
template<class T, class U = void>
struct has_member_type : std::false_type {};
template<class T>
struct has_member_type<T, typename Void<typename T::type>::type > : std::true_type {
  struct wrap {
    typedef typename T::type type;
    typedef typename T::type ReturnType;
  };
};

} // namespace detail

/// Promote two types - default according to common_type
template <class T, class U, class V = void> struct promote : std::common_type<T,U> {};

/// Can types A and B be promoted to a common type?
template<class A, class B> using can_promote = detail::has_member_type<promote<A,B>>;

template <class T, class U>
struct promote<dualt2<T>, dualt2<U>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualt2_traits<T>::depth == (int)dualt2_traits<U>::depth)>::type>
{
  typedef dualt2<typename promote<U,T>::type> type;
};
template <class T, class U>
struct promote<dualt2<T>, U,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualt2_traits<T>::depth >= (int)dualt2_traits<U>::depth
                                        && !is_complex<U>::value)>::type>
{
  typedef dualt2<typename promote<U,T>::type> type;
};
template <class T, class U>
struct promote<U, dualt2<T>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualt2_traits<T>::depth >= (int)dualt2_traits<U>::depth
                                        && !is_complex<U>::value)>::type>
{
  typedef dualt2<typename promote<U,T>::type> type;
};
// /////////////////////////////////////////////////
template <class T, class U>
struct promote<std::complex<T>, std::complex<U>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualt2<T>::value || is_dualt2<U>::value))>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
template <class T, class U>
struct promote<std::complex<T>, U,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualt2<T>::value || is_dualt2<U>::value)
                                        && !is_complex<U>::value)>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
template <class T, class U>
struct promote<U, std::complex<T>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualt2<T>::value || is_dualt2<U>::value)
                                        && !is_complex<U>::value)>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
// /////////////////////////////////////////////////

#endif // PARSED_BY_DOXYGEN

} // namespace dualst2

#define NOMACRO

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifdef CPPDUALST2_ENABLE_STD_IS_ARITHMETIC

/// Duals are as arithmetic as their value_type is arithmetic.
template <class T>
struct is_arithmetic<dualst2::dualt2<T>> : is_arithmetic<T> {};

#endif // CPPDUALST2_ENABLE_IS_ARITHMETIC

/// Duals are compound types.
template <class T>
struct is_compound<dualst2::dualt2<T>> : true_type {};

// Modification of std::numeric_limits<> per
// C++03 17.4.3.1/1, and C++11 18.3.2.3/1.
template <class T>
struct numeric_limits<dualst2::dualt2<T>> : numeric_limits<T> {
  static constexpr bool is_specialized = true;
  static constexpr dualst2::dualt2<T> min NOMACRO ()  { return numeric_limits<T>::min NOMACRO (); }
  static constexpr dualst2::dualt2<T> lowest()        { return numeric_limits<T>::lowest(); }
  static constexpr dualst2::dualt2<T> max NOMACRO ()  { return numeric_limits<T>::max NOMACRO (); }
  static constexpr dualst2::dualt2<T> epsilon()       { return numeric_limits<T>::epsilon(); }
  static constexpr dualst2::dualt2<T> round_error()   { return numeric_limits<T>::round_error(); }
  static constexpr dualst2::dualt2<T> infinity()      { return numeric_limits<T>::infinity(); }
  static constexpr dualst2::dualt2<T> quiet_NaN()     { return numeric_limits<T>::quiet_NaN(); }
  static constexpr dualst2::dualt2<T> signaling_NaN() { return numeric_limits<T>::signaling_NaN(); }
  static constexpr dualst2::dualt2<T> denorm_min()    { return numeric_limits<T>::denorm_min(); }
};

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_END_NAMESPACE_STD
#else
} // namespace std
#endif

namespace dualst2 {

#ifndef PARSED_BY_DOXYGEN

// T and X are wrapped in a dualt2<>
#define CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X)                              \
  typename std::enable_if<(int)dualst2::dualt2_traits<X>::depth ==          \
                          (int)dualst2::dualt2_traits<T>::depth, int>::type = 0,   \
    typename std::enable_if<can_promote<T,X>::value,int>::type = 0

// Both T and U are wrapped in a dualt2<>
#define CPPDUALST2_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)                    \
  typename std::enable_if<(int)dualst2::dualt2_traits<T>::depth ==          \
                          (int)dualst2::dualt2_traits<U>::depth, int>::type = 0, \
    typename std::enable_if<can_promote<T,U>::value,int>::type = 0,     \
    typename common_t = dualt2<typename dualst2::promote<T,U>::type>

// T is wrapped in a dualt2<>, U may or may not be.
#define CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)                     \
  typename std::enable_if<((int)dualst2::dualt2_traits<T>::depth >=         \
                           (int)dualst2::dualt2_traits<U>::depth), int>::type = 0, \
    typename std::enable_if<can_promote<dualt2<T>,U>::value,int>::type = 0, \
    typename common_t = typename dualst2::promote<dualt2<T>, U>::type

// T is wrapped in complex<dualt2<>>
#define CPPDUALST2_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U)                  \
  typename std::enable_if<((int)dualst2::dualt2_traits<T>::depth >=         \
                           (int)dualst2::dualt2_traits<U>::depth), int>::type = 0, \
  typename std::enable_if<can_promote<std::complex<dualt2<T>>,U>::value,int>::type = 0, \
    typename common_t = typename dualst2::promote<std::complex<dualt2<T> >,U>::type

#define CPPDUALST2_ENABLE_IF(...) typename std::enable_if< (__VA_ARGS__) , int>::type = 0

#endif

template<class T>
class dualt2
{
public:
  typedef T value_type;

  value_type _real;
  value_type _ep_t;
  value_type _ep_t2;

  /// Construct dualt2 from optional real and dualt2 parts.
  constexpr
  dualt2(const value_type re = value_type(), const value_type t = value_type(), const value_type t2 = value_type())
    : _real(re), _ep_t(t), _ep_t2(t2) {}

  /// Copy construct from a dualt2 of equal depth.
  template<class X, CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X),
           CPPDUALST2_ENABLE_IF(!is_complex<X>::value)>
  dualt2(const dualt2<X> & x)
    : _real((T)x._real), _ep_t((T)x._ep_t), _ep_t2((T)x._ep_t2) {}

  /// Explicit cast to an arithmetic type retains the rpart()
  template <class X,
            CPPDUALST2_ENABLE_IF(std::is_arithmetic<X>::value && !is_dualt2<X>::value)>
  explicit operator X() const { return X(_real); }

  /// Get the real part.
  T rpart() const { return _real; }

  /// Set the real part.
  void rpart(value_type re) { _real = re; }

  /// Unary negation
  dualt2<T> operator-() const { return dualt2<T>(-_real, -_ep_t, -_ep_t2); }

  /// Unary nothing
  dualt2<T> operator+() const { return *this; }

  /// Assignment of `value_type` assigns the real part and zeros the dualt2 part.
  dualt2<T> & operator= (const T & x) { _real =  x; _ep_t = value_type(); _ep_t2 = value_type(); return *this; }

  /// Add a relatively-scalar to this dualt2.
  dualt2<T> & operator+=(const T & x) { _real += x; return *this; }

  /// Subtract a relatively-scalar from this dualt2.
  dualt2<T> & operator-=(const T & x) { _real -= x; return *this; }

  /// Multiply a relatively-scalar with this dualt2.
  dualt2<T> & operator*=(const T & x) { _real *= x; _ep_t *= x; _ep_t2 *= x; return *this; }

  /// Divide this dualt2 by relatively-scalar.
  dualt2<T> & operator/=(const T & x) { _real /= x; _ep_t /= x; _ep_t2 /= x; return *this; }


  // NOTE: Added by Ben Ruijl
  dualt2<T> operator+(const T & x) const { return dualt2<T>(_real + x, _ep_t, _ep_t2); }
  dualt2<T> operator-(const T & x) const { return dualt2<T>(_real - x, _ep_t, _ep_t2); }
  dualt2<T> operator*(const T & x) const { return dualt2<T>(_real * x, _ep_t * x, _ep_t2 * x); }
  dualt2<T> operator/(const T & x) const { return dualt2<T>(_real / x, _ep_t / x, _ep_t2 / x); }

  template<class X, CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualt2<T> & operator= (const dualt2<X> & x) { _real =  x._real; _ep_t =  x._ep_t; _ep_t2 =  x._ep_t2; return *this; }

  /// Add a dualt2 of the same depth to this dualt2.
  template<class X, CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualt2<T> & operator+=(const dualt2<X> & x) { _real += x._real; _ep_t +=  x._ep_t; _ep_t2 +=  x._ep_t2; return *this; }

  /// Subtract a dualt2 of the same depth from this dualt2.
  template<class X, CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualt2<T> & operator-=(const dualt2<X> & x) { _real -= x._real; _ep_t -=  x._ep_t; _ep_t2 -=  x._ep_t2; return *this; }

  /// Multiply this dualt2 with a dualt2 of same depth.
  template<class X, CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualt2<T> & operator*=(const dualt2<X> & rhs) {
    T real = _real * rhs._real;
    T ep_t = _real * rhs._ep_t + _ep_t * rhs._real;
    T ep_t2 = _ep_t * rhs._ep_t + _real * rhs._ep_t2 + _ep_t2 * rhs._real;
    _real = real;
    _ep_t = ep_t;
    _ep_t2 = ep_t2;

    return *this;
  }

  /// Divide this dualt2 by another dualt2 of the same or lower depth.
  template<class X, CPPDUALST2_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualt2<T> & operator/=(const dualt2<X> & rhs) {
    T r1 = T(1) / rhs._real;
    T r2 = r1 * r1;
    T r3 = r2 * r1;

    T real = _real * r1;
    T ep_t = _ep_t * r1 - _real * rhs._ep_t * r2;
    T ep_t2 = _ep_t2 * r1 - _ep_t * rhs._ep_t * r2 - _real * rhs._ep_t2 * r2
        + _real * rhs._ep_t * rhs._ep_t * r3;

    _real = real;
    _ep_t = ep_t;
    _ep_t2 = ep_t2;
    return *this;
  }


};

/// Get the dualt2's real part.
template <class T> T rpart(const dualt2<T> & x) { return x.rpart(); }

#ifndef PARSED_BY_DOXYGEN

/// Dual +-*/ ops with another entity
#define CPPDUALST2_BINARY_OP(op)                                          \
  template<class T, class U, CPPDUALST2_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)> \
  common_t                                                              \
  operator op(const dualt2<T> & z, const dualt2<U> & w) {                   \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
    CPPDUALST2_ENABLE_IF(!std::is_same<U,std::complex<dualt2<T>>>::value)>  \
  common_t                                                              \
  operator op(const dualt2<T> & z, const U & w) {                         \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
    CPPDUALST2_ENABLE_IF(!std::is_same<U,std::complex<dualt2<T>>>::value)>  \
  common_t                                                              \
  operator op(const U & z, const dualt2<T> & w) {                         \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U), \
    CPPDUALST2_ENABLE_IF(!std::is_same<U,std::complex<dualt2<T>>>::value)>  \
  common_t                                                              \
  operator op(const std::complex<dualt2<T>> & z, const U & w) {           \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U), \
    CPPDUALST2_ENABLE_IF(!std::is_same<U,std::complex<dualt2<T>>>::value)>  \
  common_t                                                              \
  operator op(const U & z, const std::complex<dualt2<T>> & w) {           \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                \

CPPDUALST2_BINARY_OP(+)
CPPDUALST2_BINARY_OP(-)
CPPDUALST2_BINARY_OP(*)
CPPDUALST2_BINARY_OP(/)

/// Dual compared to a non-complex lower rank thing
#define CPPDUALST2_LHS_COMPARISON(op)                                     \
  template<class T, class U, CPPDUALST2_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)> \
  bool                                                                  \
  operator op(const dualt2<T> & a, const dualt2<U> & b) {                   \
    return a.rpart() op b.rpart();                                      \
  }                                                                     \
  template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
           CPPDUALST2_ENABLE_IF(!is_complex<U>::value)>                   \
  bool                                                                  \
  operator op(const U & a, const dualt2<T> & b) {                         \
    return a op b.rpart();                                              \
  }                                                                     \
  template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
           CPPDUALST2_ENABLE_IF(!is_complex<U>::value)>                   \
  bool                                                                  \
  operator op(const dualt2<T> & a, const U & b) {                         \
    return a.rpart() op b;                                              \
  }

CPPDUALST2_LHS_COMPARISON(<)
CPPDUALST2_LHS_COMPARISON(>)
CPPDUALST2_LHS_COMPARISON(<=)
CPPDUALST2_LHS_COMPARISON(>=)
CPPDUALST2_LHS_COMPARISON(==)
CPPDUALST2_LHS_COMPARISON(!=)

#endif // PARSED_BY_DOXYGEN

// NOTE: Added by Ben Ruijl
template<class T> dualt2<T> operator*(const T & lhs, const dualt2<T> rhs) { return rhs * lhs; }
template<class T> dualt2<T> operator/(const T & lhs, const dualt2<T> rhs) { return rhs.inv() * lhs; }
template<class T> dualt2<T> operator+(const T & lhs, const dualt2<T> rhs) { return rhs + lhs; }
template<class T> dualt2<T> operator-(const T & lhs, const dualt2<T> rhs) { return -rhs + lhs; }

/// Exponential e^x
template<class T> dualt2<T> exp(const dualt2<T> & x) {
  using std::exp;
  T r = exp(x._real);
  return dualt2<T>(
      r,
      r * x._ep_t,
      r * x._ep_t2 + r * x._ep_t * x._ep_t / T(2)
  );
}

/// Natural log ln(x)
template<class T> dualt2<T> log(const dualt2<T> & x) {
  using std::log;
  T v = log(x.rpart());
  assert(false);
}

template<class T> dualt2<T> log10(const dualt2<T> & x) {
  using std::log;
  return log(x) / log(static_cast<T>(10));
}

template<class T> dualt2<T> log2(const dualt2<T> & x) {
  using std::log;
  return log(x) / log(static_cast<T>(2));
}

template<class T, class U, CPPDUALST2_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)>
common_t
pow(const dualt2<T> & f, const dualt2<U> & g) {
  using std::pow;
  using std::log;
  /*T v = pow(f.rpart(), g.rpart());
  common_t(v,
                  pow(f.rpart(), g.rpart() - T(1)) *
                  (g.rpart() * f.dpart()
                   + f.rpart() * log(f.rpart()) * g.dpart()));*/
  assert(false);
}

template<class T> dualt2<T> pow(const dualt2<T> & x, int y) {
  using std::pow;
  T rmm = pow(x._real, y - 2);
  T rm = rmm * x._real;
  T r = rm * x._real;
  T dr = T(y) * rm;
  T ddr = T(y) * rmm;
  return dualt2<T>(
      r,
      x._ep_t * dr,
      x._ep_t2 * dr
          + x._ep_t * x._ep_t / T(2) * (T(y) * ddr - ddr)
  );
}


template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
pow(const U & x, const dualt2<T> & y) {
  return pow(common_t(x), y);
}

namespace utils {
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
}

template<class T> dualt2<T> abs(const dualt2<T> & x) {
  using std::abs;
  T s = utils::sgn(x.rpart());
  return dualt2<T>(abs(x._real), x._ep_t * s, x._ep_t2 * s);
}

template<class T> dualt2<T> fabs(const dualt2<T> & x) {
  using std::fabs;
  T s = utils::sgn(x.rpart());
  return dualt2<T>(fabs(x._real), x._ep_t * s, x._ep_t2 * s);
}

template<class T> dualst2::dualt2<T> copysign(const dualst2::dualt2<T> & x, const dualst2::dualt2<T> & y) {
  using std::copysign;
  T r = copysign(x.rpart(), y.rpart());
  //return dualst2::dualt2<T>(r, r == x.rpart() ? x.dpart() : -x.dpart());
  assert(false);
}

template<class T> dualst2::dualt2<T> hypot(const dualst2::dualt2<T> & x, const dualst2::dualt2<T> & y) {
  return sqrt(x*x + y*y);
}

template<class T> dualst2::dualt2<T> scalbn(const dualst2::dualt2<T> & arg, int ex) {
  return arg * std::pow((T)2, ex);
}

template<class T> dualst2::dualt2<T> (fmax)(const dualst2::dualt2<T> & x, const dualst2::dualt2<T> & y) {
  return x.rpart() > y.rpart() ? x : y;
}

template<class T> dualst2::dualt2<T> (fmin)(const dualst2::dualt2<T> & x, const dualst2::dualt2<T> & y) {
  return x.rpart() <= y.rpart() ? x : y;
}

template<class T> dualst2::dualt2<T> logb(const dualst2::dualt2<T> & x) {
  return dualst2::log2(x);
}

template<class T> int (fpclassify)(const dualst2::dualt2<T> & d) { using std::fpclassify; return (fpclassify)(d.rpart()); }
template<class T> bool (isfinite)(const dualst2::dualt2<T> & d) { using std::isfinite; return (isfinite)(d.rpart()); }
template<class T> bool (isnormal)(const dualst2::dualt2<T> & d) { using std::isnormal; return (isnormal)(d.rpart()); }
template<class T> bool (isinf)(const dualst2::dualt2<T> & d) { using std::isinf; return (isinf)(d.rpart()); }
template<class T> bool (isnan)(const dualst2::dualt2<T> & d) { using std::isnan; return (isnan)(d.rpart()); }
template<class T> bool (signbit)(const dualst2::dualt2<T> & d) { using std::signbit; return (signbit)(d.rpart()); }

template<class T> dualt2<T> sqrt(const dualt2<T> & x) {
  using std::sqrt;
  T half = T(1) / T(2);
  T r = sqrt(x._real);
  T ir = T(1) / r * half;
  T irr = ir / x._real * half;
  
  return dualt2<T>(r, x._ep_t * ir, x._ep_t2 * ir - x._ep_t * x._ep_t * irr * half);
}

template<class T> dualt2<T> cbrt(const dualt2<T> & x) {
  assert(false);
}

template<class T> dualt2<T> sin(const dualt2<T> & x) {
  using std::sin;
  using std::cos;
  T c = cos(x._real);
  T s = sin(x._real);
  
  return dualt2(
      s,
      x._ep_t * c,
      x._ep_t2 * c - x._ep_t * x._ep_t * s / T(2)
  );
}

template<class T> dualt2<T> cos(const dualt2<T> & x) {
  using std::cos;
  using std::sin;
  T c = cos(x._real);
  T s = sin(x._real);
  return dualt2(
      c,
      -x._ep_t * s,
      -x._ep_t2 * s - x._ep_t * x._ep_t * c / T(2)
  );
}

template<class T> dualt2<T> tan(const dualt2<T> & x) {
  using std::tan;
  T v = tan(x.rpart());
  assert(false);
}

template<class T> dualt2<T> asin(const dualt2<T> & x) {
  using std::asin;
  using std::sqrt;
  T v = asin(x.rpart());
  assert(false);
}

template<class T> dualt2<T> acos(const dualt2<T> & x) {
  using std::acos;
  using std::sqrt;
  assert(false);
}

template<class T> dualt2<T> atan(const dualt2<T> & x) {
  using std::atan;
  assert(false);
}

template<class T> dualt2<T> atan2(const dualt2<T> & y, const dualt2<T> & x) {
  using std::atan2;
  T v = atan2(y.rpart(), x.rpart());
  assert(false);
}

// TODO
template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
atan2(const dualt2<T> & y, const U & x) {
  using std::atan2;
  T v = atan2(y.rpart(), x);
  assert(false);
}

// TODO
template<class T, class U, CPPDUALST2_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
atan2(const U & y, const dualt2<T> & x) {
  using std::atan2;
  T v = atan2(y, x.rpart());
  assert(false);
}

// TODO
template<class T> dualt2<T> sinh(const dualt2<T> & x);
template<class T> dualt2<T> cosh(const dualt2<T> & x);
template<class T> dualt2<T> tanh(const dualt2<T> & x);
template<class T> dualt2<T> asinh(const dualt2<T> & x);
template<class T> dualt2<T> acosh(const dualt2<T> & x);
template<class T> dualt2<T> atanh(const dualt2<T> & x);
template<class T> dualt2<T> log1p(const dualt2<T> & x);
template<class T> dualt2<T> expm1(const dualt2<T> & x);

/// The error function.  Make sure to `#include <math.h>` before
/// `#include <dualt2s/dualt2>` to use this function.
template<class T> dualt2<T> erf(const dualt2<T> & x) {
  using std::erf;
  using std::sqrt;
  using std::pow;
  using std::exp;
  assert(false);
}

/// Error function complement (1 - erf()).
template<class T> dualt2<T> erfc(const dualt2<T> & x) {
  using std::erfc;
  using std::sqrt;
  using std::pow;
  using std::exp;
  assert(false);
}

/// Gamma function.  Approximation of the dualt2 part.
// TODO specialize for integers
template<class T> dualt2<T> tgamma(const dualt2<T> & x) {
  using std::tgamma;
  assert(false);
}

/// Log of absolute value of gamma function.  Approximation of the dualt2 part.
template<class T> dualt2<T> lgamma(const dualt2<T> & x) {
  using std::lgamma;
  T v = lgamma(x.rpart());
  assert(false);
}

/// Putto operator
template<class T, class _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits> &
operator<<(std::basic_ostream<_CharT, _Traits> & os, const dualt2<T> & x)
{
  std::basic_ostringstream<_CharT, _Traits> s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision());
  s << '(' << x._real
    << " + "
    << x._ep_t
    << "_e_t"
    << " + "
    << x._ep_t2
    << "_e_t2"
    << ")";
  return os << s.str();
}

} // namespace dualst2

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifndef PARSED_BY_DOXYGEN

#define make_math(T)                                                    \
  inline T (frexp)(const dualst2::dualt2<T> & arg, int* exp ) { return (frexp)(arg.rpart(), exp); } \
  inline dualst2::dualt2<T> (ldexp)(const dualst2::dualt2<T> & arg, int exp ) { return arg * std::pow((T)2,exp); } \
  inline T (trunc)(const dualst2::dualt2<T> & d) { return (trunc)(d.rpart()); } \
  inline T (floor)(const dualst2::dualt2<T> & d) { return (floor)(d.rpart()); } \
  inline T (ceil)(const dualst2::dualt2<T> & d)  { return (ceil)(d.rpart()); } \
  inline T (round)(const dualst2::dualt2<T> & d) { return (round)(d.rpart()); } \
  inline int (fpclassify)(const dualst2::dualt2<T> & d) { return (fpclassify)(d.rpart()); } \
  inline bool (isfinite)(const dualst2::dualt2<T> & d) { return (isfinite)(d.rpart()); } \
  inline bool (isnormal)(const dualst2::dualt2<T> & d) { return (isnormal)(d.rpart()); } \
  inline bool (isinf)(const dualst2::dualt2<T> & d) { return (isinf)(d.rpart()); } \
  inline bool (isnan)(const dualst2::dualt2<T> & d) { return (isnan)(d.rpart()); } \

//make_math(float)
make_math(double)
make_math(long double)

#undef make_math

#endif // PARSED_BY_DOXYGEN

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_END_NAMESPACE_STD
#else
} // namespace std
#endif

#endif // DUALT2
