//===-- duals/dualkt3 - Dual number class --------------------------*- C++ -*-===//
//
// Extension of cppduals project.
// https://tesch1.gitlab.io/cppduals
//
// (c)2019 Michael Tesch. tesch1@gmail.com
//
// See https://gitlab.com/tesch1/cppdualkt3s/blob/master/LICENSE.txt for
// license information.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DUALKT3
#define DUALKT3

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

#if !defined(CPPDUALSKT3_IGNORE_COMPILER_VERSION) && !defined(_WIN32)
#if __cplusplus < 201103L
  #error CPPDUALSKT3 needs at least a C++11 compliant compiler
#endif
#endif

namespace dualskt3 {

#ifndef PARSED_BY_DOXYGEN
template<class T> class dualkt3;
#endif

/// Check if T is dualkt3<>, match non-dualkt3s.
template <class T> struct is_dualkt3 : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is dualkt3<>, match dualkt3<>.
template <class T> struct is_dualkt3<dualkt3<T> > : std::true_type {};

#endif

/// Check if T is std::complex<>, match non- std::complex<>.
template <class T> struct is_complex : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is std::complex<>, match std::complex<>.
template <class T> struct is_complex<std::complex<T> > : std::true_type {};

#endif

/// dualkt3_traits helper class.
template <class T> struct dualkt3_traits
{
  /// Depth of T - for T=scalar this is 0. for dualkt3_traits<double> it
  /// is 1.
  enum { depth = 0 }; // -Wenum-compare

  /// The real storage type.
  typedef T real_type;
};

#ifndef PARSED_BY_DOXYGEN

/// dualkt3_traits for dualkt3<> types
template <class T> struct dualkt3_traits<dualkt3<T>>
{
  /// Depth to which this dualkt3<> type is nested.  One (1) is a
  /// first-level dualkt3, whereas non-dualkt3s have a depth of 0.
  enum { depth = dualkt3_traits<T>::depth + 1 };

  /// The real storage type.
  typedef typename dualkt3_traits<T>::real_type real_type;
};

template <class T> struct dualkt3_traits<std::complex<dualkt3<T>>>
{
  /// complex<dualkt3<T>> have the same 'depth' as their dual.
  enum { depth = dualkt3_traits<T>::depth };

  /// The real storage type.
  typedef typename dualkt3_traits<T>::real_type real_type;
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
struct promote<dualkt3<T>, dualkt3<U>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualkt3_traits<T>::depth == (int)dualkt3_traits<U>::depth)>::type>
{
  typedef dualkt3<typename promote<U,T>::type> type;
};
template <class T, class U>
struct promote<dualkt3<T>, U,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualkt3_traits<T>::depth >= (int)dualkt3_traits<U>::depth
                                        && !is_complex<U>::value)>::type>
{
  typedef dualkt3<typename promote<U,T>::type> type;
};
template <class T, class U>
struct promote<U, dualkt3<T>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualkt3_traits<T>::depth >= (int)dualkt3_traits<U>::depth
                                        && !is_complex<U>::value)>::type>
{
  typedef dualkt3<typename promote<U,T>::type> type;
};
// /////////////////////////////////////////////////
template <class T, class U>
struct promote<std::complex<T>, std::complex<U>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualkt3<T>::value || is_dualkt3<U>::value))>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
template <class T, class U>
struct promote<std::complex<T>, U,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualkt3<T>::value || is_dualkt3<U>::value)
                                        && !is_complex<U>::value)>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
template <class T, class U>
struct promote<U, std::complex<T>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualkt3<T>::value || is_dualkt3<U>::value)
                                        && !is_complex<U>::value)>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
// /////////////////////////////////////////////////

#endif // PARSED_BY_DOXYGEN

} // namespace dualskt3

#define NOMACRO

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifdef CPPDUALSKT3_ENABLE_STD_IS_ARITHMETIC

/// Duals are as arithmetic as their value_type is arithmetic.
template <class T>
struct is_arithmetic<dualskt3::dualkt3<T>> : is_arithmetic<T> {};

#endif // CPPDUALSKT3_ENABLE_IS_ARITHMETIC

/// Duals are compound types.
template <class T>
struct is_compound<dualskt3::dualkt3<T>> : true_type {};

// Modification of std::numeric_limits<> per
// C++03 17.4.3.1/1, and C++11 18.3.2.3/1.
template <class T>
struct numeric_limits<dualskt3::dualkt3<T>> : numeric_limits<T> {
  static constexpr bool is_specialized = true;
  static constexpr dualskt3::dualkt3<T> min NOMACRO ()  { return numeric_limits<T>::min NOMACRO (); }
  static constexpr dualskt3::dualkt3<T> lowest()        { return numeric_limits<T>::lowest(); }
  static constexpr dualskt3::dualkt3<T> max NOMACRO ()  { return numeric_limits<T>::max NOMACRO (); }
  static constexpr dualskt3::dualkt3<T> epsilon()       { return numeric_limits<T>::epsilon(); }
  static constexpr dualskt3::dualkt3<T> round_error()   { return numeric_limits<T>::round_error(); }
  static constexpr dualskt3::dualkt3<T> infinity()      { return numeric_limits<T>::infinity(); }
  static constexpr dualskt3::dualkt3<T> quiet_NaN()     { return numeric_limits<T>::quiet_NaN(); }
  static constexpr dualskt3::dualkt3<T> signaling_NaN() { return numeric_limits<T>::signaling_NaN(); }
  static constexpr dualskt3::dualkt3<T> denorm_min()    { return numeric_limits<T>::denorm_min(); }
};

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_END_NAMESPACE_STD
#else
} // namespace std
#endif

namespace dualskt3 {

#ifndef PARSED_BY_DOXYGEN

// T and X are wrapped in a dualkt3<>
#define CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X)                              \
  typename std::enable_if<(int)dualskt3::dualkt3_traits<X>::depth ==          \
                          (int)dualskt3::dualkt3_traits<T>::depth, int>::type = 0,   \
    typename std::enable_if<can_promote<T,X>::value,int>::type = 0

// Both T and U are wrapped in a dualkt3<>
#define CPPDUALSKT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)                    \
  typename std::enable_if<(int)dualskt3::dualkt3_traits<T>::depth ==          \
                          (int)dualskt3::dualkt3_traits<U>::depth, int>::type = 0, \
    typename std::enable_if<can_promote<T,U>::value,int>::type = 0,     \
    typename common_t = dualkt3<typename dualskt3::promote<T,U>::type>

// T is wrapped in a dualkt3<>, U may or may not be.
#define CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)                     \
  typename std::enable_if<((int)dualskt3::dualkt3_traits<T>::depth >=         \
                           (int)dualskt3::dualkt3_traits<U>::depth), int>::type = 0, \
    typename std::enable_if<can_promote<dualkt3<T>,U>::value,int>::type = 0, \
    typename common_t = typename dualskt3::promote<dualkt3<T>, U>::type

// T is wrapped in complex<dualkt3<>>
#define CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U)                  \
  typename std::enable_if<((int)dualskt3::dualkt3_traits<T>::depth >=         \
                           (int)dualskt3::dualkt3_traits<U>::depth), int>::type = 0, \
  typename std::enable_if<can_promote<std::complex<dualkt3<T>>,U>::value,int>::type = 0, \
    typename common_t = typename dualskt3::promote<std::complex<dualkt3<T> >,U>::type

#define CPPDUALSKT3_ENABLE_IF(...) typename std::enable_if< (__VA_ARGS__) , int>::type = 0

#endif

template<class T>
class dualkt3
{
public:
  typedef T value_type;

  value_type _real;
  value_type _ep_k0;
  value_type _ep_t;
  value_type _ep_k0_t;
  value_type _ep_t2;
  value_type _ep_k0_t2;
  value_type _ep_t3;

  /// Construct dualkt3 from optional real and dualkt3 parts.
  constexpr
  dualkt3(const value_type re = value_type(), const value_type k0 = value_type(), const value_type t = value_type(), const value_type k0_t = value_type(), const value_type t2 = value_type(), const value_type k0_t2 = value_type(), const value_type t3 = value_type())
    : _real(re), _ep_k0(k0), _ep_t(t), _ep_k0_t(k0_t), _ep_t2(t2), _ep_k0_t2(k0_t2), _ep_t3(t3) {}

  /// Copy construct from a dualkt3 of equal depth.
  template<class X, CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X),
           CPPDUALSKT3_ENABLE_IF(!is_complex<X>::value)>
  dualkt3(const dualkt3<X> & x)
    : _real((T)x._real), _ep_k0((T)x._ep_k0), _ep_t((T)x._ep_t), _ep_k0_t((T)x._ep_k0_t), _ep_t2((T)x._ep_t2), _ep_k0_t2((T)x._ep_k0_t2), _ep_t3((T)x._ep_t3) {}

  /// Explicit cast to an arithmetic type retains the rpart()
  template <class X,
            CPPDUALSKT3_ENABLE_IF(std::is_arithmetic<X>::value && !is_dualkt3<X>::value)>
  explicit operator X() const { return X(_real); }

  /// Get the real part.
  T rpart() const { return _real; }

  /// Set the real part.
  void rpart(value_type re) { _real = re; }

  /// Unary negation
  dualkt3<T> operator-() const { return dualkt3<T>(-_real, -_ep_k0, -_ep_t, -_ep_k0_t, -_ep_t2, -_ep_k0_t2, -_ep_t3); }

  /// Unary nothing
  dualkt3<T> operator+() const { return *this; }

  /// Assignment of `value_type` assigns the real part and zeros the dualkt3 part.
  dualkt3<T> & operator= (const T & x) { _real =  x; _ep_k0 = value_type(); _ep_t = value_type(); _ep_k0_t = value_type(); _ep_t2 = value_type(); _ep_k0_t2 = value_type(); _ep_t3 = value_type(); return *this; }

  /// Add a relatively-scalar to this dualkt3.
  dualkt3<T> & operator+=(const T & x) { _real += x; return *this; }

  /// Subtract a relatively-scalar from this dualkt3.
  dualkt3<T> & operator-=(const T & x) { _real -= x; return *this; }

  /// Multiply a relatively-scalar with this dualkt3.
  dualkt3<T> & operator*=(const T & x) { _real *= x; _ep_k0 *= x; _ep_t *= x; _ep_k0_t *= x; _ep_t2 *= x; _ep_k0_t2 *= x; _ep_t3 *= x; return *this; }

  /// Divide this dualkt3 by relatively-scalar.
  dualkt3<T> & operator/=(const T & x) { _real /= x; _ep_k0 /= x; _ep_t /= x; _ep_k0_t /= x; _ep_t2 /= x; _ep_k0_t2 /= x; _ep_t3 /= x; return *this; }


  // NOTE: Added by Ben Ruijl
  dualkt3<T> operator+(const T & x) const { return dualkt3<T>(_real + x, _ep_k0, _ep_t, _ep_k0_t, _ep_t2, _ep_k0_t2, _ep_t3); }
  dualkt3<T> operator-(const T & x) const { return dualkt3<T>(_real - x, _ep_k0, _ep_t, _ep_k0_t, _ep_t2, _ep_k0_t2, _ep_t3); }
  dualkt3<T> operator*(const T & x) const { return dualkt3<T>(_real * x, _ep_k0 * x, _ep_t * x, _ep_k0_t * x, _ep_t2 * x, _ep_k0_t2 * x, _ep_t3 * x); }
  dualkt3<T> operator/(const T & x) const { return dualkt3<T>(_real / x, _ep_k0 / x, _ep_t / x, _ep_k0_t / x, _ep_t2 / x, _ep_k0_t2 / x, _ep_t3 / x); }

  template<class X, CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualkt3<T> & operator= (const dualkt3<X> & x) { _real =  x._real; _ep_k0 =  x._ep_k0;_ep_t =  x._ep_t; _ep_k0_t =  x._ep_k0_t; _ep_t2 =  x._ep_t2; _ep_k0_t2 =  x._ep_k0_t2; _ep_t3 =  x._ep_t3; return *this; }

  /// Add a dualkt3 of the same depth to this dualkt3.
  template<class X, CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualkt3<T> & operator+=(const dualkt3<X> & x) { _real += x._real; _ep_k0 += x._ep_k0; _ep_t += x._ep_t; _ep_k0_t += x._ep_k0_t; _ep_t2 += x._ep_t2; _ep_k0_t2 += x._ep_k0_t2; _ep_t3 += x._ep_t3; return *this; }

  /// Subtract a dualkt3 of the same depth from this dualkt3.
  template<class X, CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualkt3<T> & operator-=(const dualkt3<X> & x) { _real -= x._real; _ep_k0 -= x._ep_k0; _ep_t -= x._ep_t; _ep_k0_t -= x._ep_k0_t; _ep_t2 -= x._ep_t2; _ep_k0_t2 -= x._ep_k0_t2; _ep_t3 -= x._ep_t3; return *this; }

  /// Multiply this dualkt3 with a dualkt3 of same depth.
  template<class X, CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualkt3<T> & operator*=(const dualkt3<X> & rhs) {
    T real = _real * rhs._real;
    T ep_k0 = _real * rhs._ep_k0 + _ep_k0 * rhs._real;
    T ep_t = _real * rhs._ep_t + _ep_t * rhs._real;
    T ep_k0_t = _real * rhs._ep_k0_t
        + _ep_k0_t * rhs._real
        + _ep_k0 * rhs._ep_t
        + _ep_t * rhs._ep_k0;
    T ep_t2 = _ep_t * rhs._ep_t + _real * rhs._ep_t2 + _ep_t2 * rhs._real;
    T ep_k0_t2 = _real * rhs._ep_k0_t2
        + rhs._real * _ep_k0_t2
        + _ep_t * rhs._ep_k0_t
        + _ep_k0_t * rhs._ep_t
        + _ep_t2 * rhs._ep_k0
        + _ep_k0 * rhs._ep_t2;
    T ep_t3 = _ep_t * rhs._ep_t2
        + _ep_t2 * rhs._ep_t
        + _real * rhs._ep_t3
        + rhs._real * _ep_t3;

    _real = real;
    _ep_k0 = ep_k0;
    _ep_t = ep_t;
    _ep_k0_t = ep_k0_t;
    _ep_t2 = ep_t2;
    _ep_k0_t2 = ep_k0_t2;
    _ep_t3 = ep_t3;
    return *this;
}

  /// Divide this dualkt3 by another dualkt3 of the same or lower depth.
  template<class X, CPPDUALSKT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualkt3<T> & operator/=(const dualkt3<X> & rhs) {
    T r1 = T(1) / rhs._real;
    T r2 = r1 * r1;
    T r3 = r2 * r1;
    T r4 = r3 * r1;
    T real = _real * r1;
    T ep_k0 = _ep_k0 * r1 - _real * rhs._ep_k0 * r2;
    T ep_t = _ep_t * r1 - _real * rhs._ep_t * r2;
    T ep_k0_t = _ep_k0_t * r1
        - (rhs._ep_k0_t * _real + _ep_k0 * rhs._ep_t + _ep_t * rhs._ep_k0) * r2
        + T(2) * _real * rhs._ep_k0 * rhs._ep_t * r3;
    T ep_t2 = _ep_t2 * r1 - (_ep_t * rhs._ep_t + _real * rhs._ep_t2) * r2
        + _real * rhs._ep_t * rhs._ep_t * r3;
    T ep_k0_t2 = _ep_k0_t2 * r1
        - r2 * (_real * rhs._ep_k0_t2
            + rhs._ep_k0_t * _ep_t
            + _ep_k0_t * rhs._ep_t
            + _ep_t2 * rhs._ep_k0
            + _ep_k0 * rhs._ep_t2)
        + r3 * (T(2)
            * _real
            * (rhs._ep_k0_t * rhs._ep_t + rhs._ep_t2 * rhs._ep_k0)
            + T(2) * _ep_t * rhs._ep_k0 * rhs._ep_t
            + _ep_k0 * rhs._ep_t * rhs._ep_t)
        - r4 * T(3)
            * _real
            * rhs._ep_k0
            * rhs._ep_t
            * rhs._ep_t;
    T ep_t3 = _ep_t3 * r1
        - (rhs._ep_t * _ep_t2 + _ep_t * rhs._ep_t2 + _real * rhs._ep_t3) * r2
        + (_ep_t * rhs._ep_t * rhs._ep_t
            + T(2) * _real * rhs._ep_t * rhs._ep_t2)
            * r3
        - _real * rhs._ep_t * rhs._ep_t * rhs._ep_t * r4;

    _real = real;
    _ep_k0 = ep_k0;
    _ep_t = ep_t;
    _ep_k0_t = ep_k0_t;
    _ep_t2 = ep_t2;
    _ep_k0_t2 = ep_k0_t2;
    _ep_t3 = ep_t3;
    return *this;
}

};

/// Get the dualkt3's real part.
template <class T> T rpart(const dualkt3<T> & x) { return x.rpart(); }

#ifndef PARSED_BY_DOXYGEN

/// Dual +-*/ ops with another entity
#define CPPDUALSKT3_BINARY_OP(op)                                          \
  template<class T, class U, CPPDUALSKT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)> \
  common_t                                                              \
  operator op(const dualkt3<T> & z, const dualkt3<U> & w) {                   \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
    CPPDUALSKT3_ENABLE_IF(!std::is_same<U,std::complex<dualkt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const dualkt3<T> & z, const U & w) {                         \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
    CPPDUALSKT3_ENABLE_IF(!std::is_same<U,std::complex<dualkt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const U & z, const dualkt3<T> & w) {                         \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U), \
    CPPDUALSKT3_ENABLE_IF(!std::is_same<U,std::complex<dualkt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const std::complex<dualkt3<T>> & z, const U & w) {           \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U), \
    CPPDUALSKT3_ENABLE_IF(!std::is_same<U,std::complex<dualkt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const U & z, const std::complex<dualkt3<T>> & w) {           \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                \

CPPDUALSKT3_BINARY_OP(+)
CPPDUALSKT3_BINARY_OP(-)
CPPDUALSKT3_BINARY_OP(*)
CPPDUALSKT3_BINARY_OP(/)

/// Dual compared to a non-complex lower rank thing
#define CPPDUALSKT3_LHS_COMPARISON(op)                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)> \
  bool                                                                  \
  operator op(const dualkt3<T> & a, const dualkt3<U> & b) {                   \
    return a.rpart() op b.rpart();                                      \
  }                                                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
           CPPDUALSKT3_ENABLE_IF(!is_complex<U>::value)>                   \
  bool                                                                  \
  operator op(const U & a, const dualkt3<T> & b) {                         \
    return a op b.rpart();                                              \
  }                                                                     \
  template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
           CPPDUALSKT3_ENABLE_IF(!is_complex<U>::value)>                   \
  bool                                                                  \
  operator op(const dualkt3<T> & a, const U & b) {                         \
    return a.rpart() op b;                                              \
  }

CPPDUALSKT3_LHS_COMPARISON(<)
CPPDUALSKT3_LHS_COMPARISON(>)
CPPDUALSKT3_LHS_COMPARISON(<=)
CPPDUALSKT3_LHS_COMPARISON(>=)
CPPDUALSKT3_LHS_COMPARISON(==)
CPPDUALSKT3_LHS_COMPARISON(!=)

#endif // PARSED_BY_DOXYGEN

// NOTE: Added by Ben Ruijl
template<class T> dualkt3<T> operator*(const T & lhs, const dualkt3<T> rhs) { return rhs * lhs; }
template<class T> dualkt3<T> operator/(const T & lhs, const dualkt3<T> rhs) { return rhs.inv() * lhs; }
template<class T> dualkt3<T> operator+(const T & lhs, const dualkt3<T> rhs) { return rhs + lhs; }
template<class T> dualkt3<T> operator-(const T & lhs, const dualkt3<T> rhs) { return -rhs + lhs; }

/// Exponential e^x
template<class T> dualkt3<T> exp(const dualkt3<T> & x) {
  using std::exp;
  T r = exp(x._real);
  return dualkt3<T>(
      r,
      r * x._ep_k0,
      r * x._ep_t,
      r * (x._ep_k0_t + x._ep_k0 * x._ep_t),
      r * (x._ep_t * x._ep_t / T(2) + x._ep_t2),
      r * (x._ep_k0_t2
              + x._ep_k0_t * x._ep_t
              + x._ep_k0
                  * (x._ep_t2 + x._ep_t * x._ep_t / T(2))),
      r * (x._ep_t * x._ep_t * x._ep_t / T(6)
              + x._ep_t3
              + x._ep_t2 * x._ep_t)
  );
}

/// Natural log ln(x)
template<class T> dualkt3<T> log(const dualkt3<T> & x) {
  using std::log;
  T v = log(x.rpart());
  assert(false);
}

template<class T> dualkt3<T> log10(const dualkt3<T> & x) {
  using std::log;
  return log(x) / log(static_cast<T>(10));
}

template<class T> dualkt3<T> log2(const dualkt3<T> & x) {
  using std::log;
  return log(x) / log(static_cast<T>(2));
}

template<class T, class U, CPPDUALSKT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)>
common_t
pow(const dualkt3<T> & f, const dualkt3<U> & g) {
  using std::pow;
  using std::log;
  /*T v = pow(f.rpart(), g.rpart());
  common_t(v,
                  pow(f.rpart(), g.rpart() - T(1)) *
                  (g.rpart() * f.dpart()
                   + f.rpart() * log(f.rpart()) * g.dpart()));*/
  assert(false);
}

template<class T> dualkt3<T> pow(const dualkt3<T> & x, int n) {
  using std::pow;
  T rmmm = pow(x._real, n - 3);
  T rmm = rmmm * x._real;
  T rm = rmm * x._real;
  T r = rm * x._real;
  T dr = T(n) * rm;
  T ddr = T(n) * rmm;
  T dddr = T(n) * rmmm;
  return dualkt3<T>(
    r,
    x._ep_k0 * dr,
    x._ep_t * dr,
    x._ep_k0_t * dr
        + x._ep_k0 * x._ep_t * T(n - 1) * ddr,
    x._ep_t2 * dr
        + x._ep_t * x._ep_t / T(2)
            * T(n-1)
            * ddr,
    x._ep_k0_t2 * dr
        + ddr
            * (x._ep_k0 * x._ep_t2 + x._ep_k0_t * x._ep_t)
            * T(n - 1)
        + dddr / T(2)
            * x._ep_k0
            * x._ep_t
            * x._ep_t
            * T((n - 1) * (n - 2)),
    x._ep_t3 * dr
        + x._ep_t * x._ep_t2 * T(n - 1) * ddr
        + x._ep_t * x._ep_t * x._ep_t / T(6)
            * T((n - 1) * (n - 2))
            * dddr
  );
}


template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
pow(const U & x, const dualkt3<T> & y) {
  return pow(common_t(x), y);
}

namespace utils {
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
}

template<class T> dualkt3<T> abs(const dualkt3<T> & x) {
  using std::abs;
  T s = utils::sgn(x.rpart());
  return dualkt3<T>(abs(x._real), x._ep_k0 * s, x._ep_t * s, x._ep_k0_t * s, x._ep_t2 * s);
}

template<class T> dualkt3<T> fabs(const dualkt3<T> & x) {
  using std::fabs;
  T s = utils::sgn(x.rpart());
  return dualkt3<T>(fabs(x._real), x._ep_k0 * s, x._ep_t * s, x._ep_k0_t * s, x._ep_t2 * s);
}

template<class T> dualskt3::dualkt3<T> copysign(const dualskt3::dualkt3<T> & x, const dualskt3::dualkt3<T> & y) {
  using std::copysign;
  T r = copysign(x.rpart(), y.rpart());
  //return dualskt3::dualkt3<T>(r, r == x.rpart() ? x.dpart() : -x.dpart());
  assert(false);
}

template<class T> dualskt3::dualkt3<T> hypot(const dualskt3::dualkt3<T> & x, const dualskt3::dualkt3<T> & y) {
  return sqrt(x*x + y*y);
}

template<class T> dualskt3::dualkt3<T> scalbn(const dualskt3::dualkt3<T> & arg, int ex) {
  return arg * std::pow((T)2, ex);
}

template<class T> dualskt3::dualkt3<T> (fmax)(const dualskt3::dualkt3<T> & x, const dualskt3::dualkt3<T> & y) {
  return x.rpart() > y.rpart() ? x : y;
}

template<class T> dualskt3::dualkt3<T> (fmin)(const dualskt3::dualkt3<T> & x, const dualskt3::dualkt3<T> & y) {
  return x.rpart() <= y.rpart() ? x : y;
}

template<class T> dualskt3::dualkt3<T> logb(const dualskt3::dualkt3<T> & x) {
  return dualskt3::log2(x);
}

template<class T> int (fpclassify)(const dualskt3::dualkt3<T> & d) { using std::fpclassify; return (fpclassify)(d.rpart()); }
template<class T> bool (isfinite)(const dualskt3::dualkt3<T> & d) { using std::isfinite; return (isfinite)(d.rpart()); }
template<class T> bool (isnormal)(const dualskt3::dualkt3<T> & d) { using std::isnormal; return (isnormal)(d.rpart()); }
template<class T> bool (isinf)(const dualskt3::dualkt3<T> & d) { using std::isinf; return (isinf)(d.rpart()); }
template<class T> bool (isnan)(const dualskt3::dualkt3<T> & d) { using std::isnan; return (isnan)(d.rpart()); }
template<class T> bool (signbit)(const dualskt3::dualkt3<T> & d) { using std::signbit; return (signbit)(d.rpart()); }

template<class T> dualkt3<T> sqrt(const dualkt3<T> & x) {
  using std::sqrt;
  T half = T(0.5);
  T r = sqrt(x._real);
  T ir = T(1) / r * half;
  T irr = ir / x._real * half;
  T irrr = irr / x._real * half;

  return dualkt3<T>(
    r,
    ir * x._ep_k0,
    ir * x._ep_t,
    ir * x._ep_k0_t - irr * x._ep_k0 * x._ep_t,
    ir * x._ep_t2 - irr * half * x._ep_t * x._ep_t,
    irrr
        * T(3)
        * half
        * x._ep_k0
        * x._ep_t
        * x._ep_t
        + ir * x._ep_k0_t2
        - irr * (x._ep_k0 * x._ep_t2 + x._ep_k0_t * x._ep_t),
    irrr * half * x._ep_t * x._ep_t * x._ep_t - irr * x._ep_t2 * x._ep_t
        + ir * x._ep_t3
  );
}

template<class T> dualkt3<T> cbrt(const dualkt3<T> & x) {
  assert(false);
}

template<class T> dualkt3<T> sin(const dualkt3<T> & x) {
  using std::sin;
  using std::cos;
  T c = cos(x._real);
  T s = sin(x._real);
  
  return dualkt3(
    s,
    c * x._ep_k0,
    c * x._ep_t,
    c * x._ep_k0_t - s * x._ep_k0 * x._ep_t,
    c * x._ep_t2 - s * x._ep_t * x._ep_t / T(2),
    c
        * (x._ep_k0_t - x._ep_k0 * x._ep_t * x._ep_t / T(2))
        - s * (x._ep_k0_t * x._ep_t + x._ep_k0 * x._ep_t2),
    c
        * (x._ep_t3
            - x._ep_t * x._ep_t * x._ep_t / T(6))
        - s * (x._ep_t * x._ep_t2)
  );
}

template<class T> dualkt3<T> cos(const dualkt3<T> & x) {
  using std::cos;
  using std::sin;
  T c = cos(x._real);
  T s = sin(x._real);
  return dualkt3(
    c,
    -s * x._ep_k0,
    -s * x._ep_t,
    -s * x._ep_k0_t - c * x._ep_k0 * x._ep_t,
    -s * x._ep_t2 - c * x._ep_t * x._ep_t / T(2),
    -s
        * (x._ep_k0_t - x._ep_k0 * x._ep_t * x._ep_t / T(2))
        - c * (x._ep_k0_t * x._ep_t + x._ep_k0 * x._ep_t2),
    -s
        * (x._ep_t3
            - x._ep_t * x._ep_t * x._ep_t / T(6))
        - c * (x._ep_t * x._ep_t2)
  );
}

template<class T> dualkt3<T> tan(const dualkt3<T> & x) {
  using std::tan;
  T v = tan(x.rpart());
  assert(false);
}

template<class T> dualkt3<T> asin(const dualkt3<T> & x) {
  using std::asin;
  using std::sqrt;
  T v = asin(x.rpart());
  assert(false);
}

template<class T> dualkt3<T> acos(const dualkt3<T> & x) {
  using std::acos;
  using std::sqrt;
  assert(false);
}

template<class T> dualkt3<T> atan(const dualkt3<T> & x) {
  using std::atan;
  assert(false);
}

template<class T> dualkt3<T> atan2(const dualkt3<T> & y, const dualkt3<T> & x) {
  using std::atan2;
  T v = atan2(y.rpart(), x.rpart());
  assert(false);
}

// TODO
template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
atan2(const dualkt3<T> & y, const U & x) {
  using std::atan2;
  T v = atan2(y.rpart(), x);
  assert(false);
}

// TODO
template<class T, class U, CPPDUALSKT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
atan2(const U & y, const dualkt3<T> & x) {
  using std::atan2;
  T v = atan2(y, x.rpart());
  assert(false);
}

// TODO
template<class T> dualkt3<T> sinh(const dualkt3<T> & x);
template<class T> dualkt3<T> cosh(const dualkt3<T> & x);
template<class T> dualkt3<T> tanh(const dualkt3<T> & x);
template<class T> dualkt3<T> asinh(const dualkt3<T> & x);
template<class T> dualkt3<T> acosh(const dualkt3<T> & x);
template<class T> dualkt3<T> atanh(const dualkt3<T> & x);
template<class T> dualkt3<T> log1p(const dualkt3<T> & x);
template<class T> dualkt3<T> expm1(const dualkt3<T> & x);

/// The error function.  Make sure to `#include <math.h>` before
/// `#include <dualkt3s/dualkt3>` to use this function.
template<class T> dualkt3<T> erf(const dualkt3<T> & x) {
  using std::erf;
  using std::sqrt;
  using std::pow;
  using std::exp;
  assert(false);
}

/// Error function complement (1 - erf()).
template<class T> dualkt3<T> erfc(const dualkt3<T> & x) {
  using std::erfc;
  using std::sqrt;
  using std::pow;
  using std::exp;
  assert(false);
}

/// Gamma function.  Approximation of the dualkt3 part.
// TODO specialize for integers
template<class T> dualkt3<T> tgamma(const dualkt3<T> & x) {
  using std::tgamma;
  assert(false);
}

/// Log of absolute value of gamma function.  Approximation of the dualkt3 part.
template<class T> dualkt3<T> lgamma(const dualkt3<T> & x) {
  using std::lgamma;
  T v = lgamma(x.rpart());
  assert(false);
}

/// Putto operator
template<class T, class _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits> &
operator<<(std::basic_ostream<_CharT, _Traits> & os, const dualkt3<T> & x)
{
  std::basic_ostringstream<_CharT, _Traits> s;
  s.flags(os.flags());
  s.imbue(os.getloc());
  s.precision(os.precision());
  s << '(' << x._real
    << " + "
    << x._ep_k0
    << "_e_k"
    << " + "
    << x._ep_t
    << "_e_t"
    << " + "
    << x._ep_k0_t
    << "_e_kt"
    << " + "
    << x._ep_t2
    << "_e_t2"
    << " + "
    << x._ep_k0_t2
    << "_e_kt2"
    << " + "
    << x._ep_t3
    << "_e_t3"
    << ")";
  return os << s.str();
}

} // namespace dualskt3

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifndef PARSED_BY_DOXYGEN

#define make_math(T)                                                    \
  inline T (frexp)(const dualskt3::dualkt3<T> & arg, int* exp ) { return (frexp)(arg.rpart(), exp); } \
  inline dualskt3::dualkt3<T> (ldexp)(const dualskt3::dualkt3<T> & arg, int exp ) { return arg * std::pow((T)2,exp); } \
  inline T (trunc)(const dualskt3::dualkt3<T> & d) { return (trunc)(d.rpart()); } \
  inline T (floor)(const dualskt3::dualkt3<T> & d) { return (floor)(d.rpart()); } \
  inline T (ceil)(const dualskt3::dualkt3<T> & d)  { return (ceil)(d.rpart()); } \
  inline T (round)(const dualskt3::dualkt3<T> & d) { return (round)(d.rpart()); } \
  inline int (fpclassify)(const dualskt3::dualkt3<T> & d) { return (fpclassify)(d.rpart()); } \
  inline bool (isfinite)(const dualskt3::dualkt3<T> & d) { return (isfinite)(d.rpart()); } \
  inline bool (isnormal)(const dualskt3::dualkt3<T> & d) { return (isnormal)(d.rpart()); } \
  inline bool (isinf)(const dualskt3::dualkt3<T> & d) { return (isinf)(d.rpart()); } \
  inline bool (isnan)(const dualskt3::dualkt3<T> & d) { return (isnan)(d.rpart()); } \

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

#endif // DUALKT3
