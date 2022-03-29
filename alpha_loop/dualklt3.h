//===-- duals/dualklt3 - Dual number class --------------------------*- C++ -*-===//
//
// Extension of cppduals project.
// https://tesch1.gitlab.io/cppduals
//
// (c)2019 Michael Tesch. tesch1@gmail.com
//
// See https://gitlab.com/tesch1/cppdualklt3s/blob/master/LICENSE.txt for
// license information.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DUALKLT3
#define DUALKLT3

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

#if !defined(CPPDUALSKLT3_IGNORE_COMPILER_VERSION) && !defined(_WIN32)
#if __cplusplus < 201103L
  #error CPPDUALSKLT3 needs at least a C++11 compliant compiler
#endif
#endif

namespace dualsklt3 {

#ifndef PARSED_BY_DOXYGEN
template<class T> class dualklt3;
#endif

/// Check if T is dualklt3<>, match non-dualklt3s.
template <class T> struct is_dualklt3 : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is dualklt3<>, match dualklt3<>.
template <class T> struct is_dualklt3<dualklt3<T> > : std::true_type {};

#endif

/// Check if T is std::complex<>, match non- std::complex<>.
template <class T> struct is_complex : std::false_type {};

#ifndef PARSED_BY_DOXYGEN

/// Check if T is std::complex<>, match std::complex<>.
template <class T> struct is_complex<std::complex<T> > : std::true_type {};

#endif

/// dualklt3_traits helper class.
template <class T> struct dualklt3_traits
{
  /// Depth of T - for T=scalar this is 0. for dualklt3_traits<double> it
  /// is 1.
  enum { depth = 0 }; // -Wenum-compare

  /// The real storage type.
  typedef T real_type;
};

#ifndef PARSED_BY_DOXYGEN

/// dualklt3_traits for dualklt3<> types
template <class T> struct dualklt3_traits<dualklt3<T>>
{
  /// Depth to which this dualklt3<> type is nested.  One (1) is a
  /// first-level dualklt3, whereas non-dualklt3s have a depth of 0.
  enum { depth = dualklt3_traits<T>::depth + 1 };

  /// The real storage type.
  typedef typename dualklt3_traits<T>::real_type real_type;
};

template <class T> struct dualklt3_traits<std::complex<dualklt3<T>>>
{
  /// complex<dualklt3<T>> have the same 'depth' as their dual.
  enum { depth = dualklt3_traits<T>::depth };

  /// The real storage type.
  typedef typename dualklt3_traits<T>::real_type real_type;
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
struct promote<dualklt3<T>, dualklt3<U>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualklt3_traits<T>::depth == (int)dualklt3_traits<U>::depth)>::type>
{
  typedef dualklt3<typename promote<U,T>::type> type;
};
template <class T, class U>
struct promote<dualklt3<T>, U,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualklt3_traits<T>::depth >= (int)dualklt3_traits<U>::depth
                                        && !is_complex<U>::value)>::type>
{
  typedef dualklt3<typename promote<U,T>::type> type;
};
template <class T, class U>
struct promote<U, dualklt3<T>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (int)dualklt3_traits<T>::depth >= (int)dualklt3_traits<U>::depth
                                        && !is_complex<U>::value)>::type>
{
  typedef dualklt3<typename promote<U,T>::type> type;
};
// /////////////////////////////////////////////////
template <class T, class U>
struct promote<std::complex<T>, std::complex<U>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualklt3<T>::value || is_dualklt3<U>::value))>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
template <class T, class U>
struct promote<std::complex<T>, U,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualklt3<T>::value || is_dualklt3<U>::value)
                                        && !is_complex<U>::value)>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
template <class T, class U>
struct promote<U, std::complex<T>,
               typename std::enable_if<(can_promote<T,U>::value
                                        && (is_dualklt3<T>::value || is_dualklt3<U>::value)
                                        && !is_complex<U>::value)>::type>
{
  typedef std::complex<typename promote<T,U>::type> type;
};
// /////////////////////////////////////////////////

#endif // PARSED_BY_DOXYGEN

} // namespace dualsklt3

#define NOMACRO

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifdef CPPDUALSKLT3_ENABLE_STD_IS_ARITHMETIC

/// Duals are as arithmetic as their value_type is arithmetic.
template <class T>
struct is_arithmetic<dualsklt3::dualklt3<T>> : is_arithmetic<T> {};

#endif // CPPDUALSKLT3_ENABLE_IS_ARITHMETIC

/// Duals are compound types.
template <class T>
struct is_compound<dualsklt3::dualklt3<T>> : true_type {};

// Modification of std::numeric_limits<> per
// C++03 17.4.3.1/1, and C++11 18.3.2.3/1.
template <class T>
struct numeric_limits<dualsklt3::dualklt3<T>> : numeric_limits<T> {
  static constexpr bool is_specialized = true;
  static constexpr dualsklt3::dualklt3<T> min NOMACRO ()  { return numeric_limits<T>::min NOMACRO (); }
  static constexpr dualsklt3::dualklt3<T> lowest()        { return numeric_limits<T>::lowest(); }
  static constexpr dualsklt3::dualklt3<T> max NOMACRO ()  { return numeric_limits<T>::max NOMACRO (); }
  static constexpr dualsklt3::dualklt3<T> epsilon()       { return numeric_limits<T>::epsilon(); }
  static constexpr dualsklt3::dualklt3<T> round_error()   { return numeric_limits<T>::round_error(); }
  static constexpr dualsklt3::dualklt3<T> infinity()      { return numeric_limits<T>::infinity(); }
  static constexpr dualsklt3::dualklt3<T> quiet_NaN()     { return numeric_limits<T>::quiet_NaN(); }
  static constexpr dualsklt3::dualklt3<T> signaling_NaN() { return numeric_limits<T>::signaling_NaN(); }
  static constexpr dualsklt3::dualklt3<T> denorm_min()    { return numeric_limits<T>::denorm_min(); }
};

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_END_NAMESPACE_STD
#else
} // namespace std
#endif

namespace dualsklt3 {

#ifndef PARSED_BY_DOXYGEN

// T and X are wrapped in a dualklt3<>
#define CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X)                              \
  typename std::enable_if<(int)dualsklt3::dualklt3_traits<X>::depth ==          \
                          (int)dualsklt3::dualklt3_traits<T>::depth, int>::type = 0,   \
    typename std::enable_if<can_promote<T,X>::value,int>::type = 0

// Both T and U are wrapped in a dualklt3<>
#define CPPDUALSKLT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)                    \
  typename std::enable_if<(int)dualsklt3::dualklt3_traits<T>::depth ==          \
                          (int)dualsklt3::dualklt3_traits<U>::depth, int>::type = 0, \
    typename std::enable_if<can_promote<T,U>::value,int>::type = 0,     \
    typename common_t = dualklt3<typename dualsklt3::promote<T,U>::type>

// T is wrapped in a dualklt3<>, U may or may not be.
#define CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)                     \
  typename std::enable_if<((int)dualsklt3::dualklt3_traits<T>::depth >=         \
                           (int)dualsklt3::dualklt3_traits<U>::depth), int>::type = 0, \
    typename std::enable_if<can_promote<dualklt3<T>,U>::value,int>::type = 0, \
    typename common_t = typename dualsklt3::promote<dualklt3<T>, U>::type

// T is wrapped in complex<dualklt3<>>
#define CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U)                  \
  typename std::enable_if<((int)dualsklt3::dualklt3_traits<T>::depth >=         \
                           (int)dualsklt3::dualklt3_traits<U>::depth), int>::type = 0, \
  typename std::enable_if<can_promote<std::complex<dualklt3<T>>,U>::value,int>::type = 0, \
    typename common_t = typename dualsklt3::promote<std::complex<dualklt3<T> >,U>::type

#define CPPDUALSKLT3_ENABLE_IF(...) typename std::enable_if< (__VA_ARGS__) , int>::type = 0

#endif

template<class T>
class dualklt3
{
public:
  typedef T value_type;

  value_type _real;
  value_type _ep_k0;
  value_type _ep_l0;
  value_type _ep_t;
  value_type _ep_k0_l0;
  value_type _ep_k0_t;
  value_type _ep_l0_t;
  value_type _ep_t2;
  value_type _ep_k0_l0_t;
  value_type _ep_k0_t2;
  value_type _ep_l0_t2;
  value_type _ep_t3;

  /// Construct dualklt3 from optional real and dualklt3 parts.
  constexpr
  dualklt3(const value_type re = value_type(), const value_type k0 = value_type(), const value_type l0 = value_type(), const value_type t = value_type(), const value_type k0_l0 = value_type(), const value_type k0_t = value_type(), const value_type l0_t = value_type(), const value_type t2 = value_type(), const value_type k0_l0_t = value_type(), const value_type k0_t2 = value_type(), const value_type l0_t2 = value_type(), const value_type t3 = value_type())
    : _real(re), _ep_k0(k0), _ep_l0(l0), _ep_t(t), _ep_k0_l0(k0_l0), _ep_k0_t(k0_t), _ep_l0_t(l0_t), _ep_t2(t2), _ep_k0_l0_t(k0_l0_t), _ep_k0_t2(k0_t2), _ep_l0_t2(l0_t2), _ep_t3(t3) {}

  /// Copy construct from a dualklt3 of equal depth.
  template<class X, CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X),
           CPPDUALSKLT3_ENABLE_IF(!is_complex<X>::value)>
  dualklt3(const dualklt3<X> & x)
    : _real((T)x._real), _ep_k0((T)x._ep_k0), _ep_l0((T)x._ep_l0), _ep_t((T)x._ep_t), _ep_k0_l0((T)x._ep_k0_l0), _ep_k0_t((T)x._ep_k0_t), _ep_l0_t((T)x._ep_l0_t), _ep_t2((T)x._ep_t2), _ep_k0_l0_t((T)x._ep_k0_l0_t), _ep_k0_t2((T)x._ep_k0_t2), _ep_l0_t2((T)x._ep_l0_t2), _ep_t3((T)x._ep_t3) {}

  /// Explicit cast to an arithmetic type retains the rpart()
  template <class X,
            CPPDUALSKLT3_ENABLE_IF(std::is_arithmetic<X>::value && !is_dualklt3<X>::value)>
  explicit operator X() const { return X(_real); }

  /// Get the real part.
  T rpart() const { return _real; }

  /// Set the real part.
  void rpart(value_type re) { _real = re; }

  /// Unary negation
  dualklt3<T> operator-() const { return dualklt3<T>(-_real, -_ep_k0, -_ep_l0, -_ep_t, -_ep_k0_l0, -_ep_k0_t, -_ep_l0_t, -_ep_t2, -_ep_k0_l0_t, -_ep_k0_t2, -_ep_l0_t2, -_ep_t3); }

  /// Unary nothing
  dualklt3<T> operator+() const { return *this; }

  /// Assignment of `value_type` assigns the real part and zeros the dualklt3 part.
  dualklt3<T> & operator= (const T & x) { _real =  x; _ep_k0 = value_type(); _ep_l0 = value_type(); _ep_t = value_type(); _ep_k0_l0 = value_type(); _ep_k0_t = value_type(); _ep_l0_t = value_type(); _ep_t2 = value_type(); _ep_k0_l0_t = value_type(); _ep_k0_t2 = value_type(); _ep_l0_t2 = value_type(); _ep_t3 = value_type(); return *this; }

  /// Add a relatively-scalar to this dualklt3.
  dualklt3<T> & operator+=(const T & x) { _real += x; return *this; }

  /// Subtract a relatively-scalar from this dualklt3.
  dualklt3<T> & operator-=(const T & x) { _real -= x; return *this; }

  /// Multiply a relatively-scalar with this dualklt3.
  dualklt3<T> & operator*=(const T & x) { _real *= x; _ep_k0 *= x; _ep_l0 *= x; _ep_t *= x; _ep_k0_l0 *= x; _ep_k0_t *= x; _ep_l0_t *= x; _ep_t2 *= x; _ep_k0_l0_t *= x; _ep_k0_t2 *= x; _ep_l0_t2 *= x; _ep_t3 *= x; return *this; }

  /// Divide this dualklt3 by relatively-scalar.
  dualklt3<T> & operator/=(const T & x) { _real /= x; _ep_k0 /= x; _ep_l0 /= x; _ep_t /= x; _ep_k0_l0 /= x; _ep_k0_t /= x; _ep_l0_t /= x; _ep_t2 /= x; _ep_k0_l0_t /= x; _ep_k0_t2 /= x; _ep_l0_t2 /= x; _ep_t3 /= x; return *this; }


  // NOTE: Added by Ben Ruijl
  dualklt3<T> operator+(const T & x) const { return dualklt3<T>(_real + x, _ep_k0, _ep_l0, _ep_t, _ep_k0_l0, _ep_k0_t, _ep_l0_t, _ep_t2, _ep_k0_l0_t, _ep_k0_t2, _ep_l0_t2, _ep_t3); }
  dualklt3<T> operator-(const T & x) const { return dualklt3<T>(_real - x, _ep_k0, _ep_l0, _ep_t, _ep_k0_l0, _ep_k0_t, _ep_l0_t, _ep_t2, _ep_k0_l0_t, _ep_k0_t2, _ep_l0_t2, _ep_t3); }
  dualklt3<T> operator*(const T & x) const { return dualklt3<T>(_real * x, _ep_k0 * x, _ep_l0 * x, _ep_t * x, _ep_k0_l0 * x, _ep_k0_t * x, _ep_l0_t * x, _ep_t2 * x, _ep_k0_l0_t * x, _ep_k0_t2 * x, _ep_l0_t2 * x, _ep_t3 * x); }
  dualklt3<T> operator/(const T & x) const { return dualklt3<T>(_real / x, _ep_k0 / x, _ep_l0 / x, _ep_t / x, _ep_k0_l0 / x, _ep_k0_t / x, _ep_l0_t / x, _ep_t2 / x, _ep_k0_l0_t / x, _ep_k0_t2 / x, _ep_l0_t2 / x, _ep_t3 / x); }

  template<class X, CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualklt3<T> & operator= (const dualklt3<X> & x) { _real =  x._real; _ep_k0 =  x._ep_k0; _ep_l0 =  x._ep_l0; _ep_t =  x._ep_t; _ep_k0_l0 =  x._ep_k0_l0; _ep_k0_t =  x._ep_k0_t; _ep_l0_t =  x._ep_l0_t; _ep_t2 =  x._ep_t2; _ep_k0_l0_t =  x._ep_k0_l0_t; _ep_k0_t2 =  x._ep_k0_t2; _ep_l0_t2 =  x._ep_l0_t2; _ep_t3 =  x._ep_t3; return *this; }

  /// Add a dualklt3 of the same depth to this dualklt3.
  template<class X, CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualklt3<T> & operator+=(const dualklt3<X> & x) { _real += x._real; _ep_k0 += x._ep_k0; _ep_l0 += x._ep_l0; _ep_t += x._ep_t; _ep_k0_l0 += x._ep_k0_l0; _ep_k0_t += x._ep_k0_t; _ep_l0_t += x._ep_l0_t; _ep_t2 += x._ep_t2; _ep_k0_l0_t += x._ep_k0_l0_t; _ep_k0_t2 += x._ep_k0_t2; _ep_l0_t2 += x._ep_l0_t2; _ep_t3 += x._ep_t3; return *this; }

  /// Subtract a dualklt3 of the same depth from this dualklt3.
  template<class X, CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualklt3<T> & operator-=(const dualklt3<X> & x) { _real -= x._real; _ep_k0 -= x._ep_k0; _ep_l0 -= x._ep_l0; _ep_t -= x._ep_t; _ep_k0_l0 -= x._ep_k0_l0; _ep_k0_t -= x._ep_k0_t; _ep_l0_t -= x._ep_l0_t; _ep_t2 -= x._ep_t2; _ep_k0_l0_t -= x._ep_k0_l0_t; _ep_k0_t2 -= x._ep_k0_t2; _ep_l0_t2 -= x._ep_l0_t2; _ep_t3 -= x._ep_t3; return *this; }

  /// Multiply this dualklt3 with a dualklt3 of same depth.
  template<class X, CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualklt3<T> & operator*=(const dualklt3<X> & rhs) {
    T real = _real * rhs._real;
    T ep_k0 = _real * rhs._ep_k0 + _ep_k0 * rhs._real;
    T ep_l0 = _real * rhs._ep_l0 + _ep_l0 * rhs._real;
    T ep_t = _real * rhs._ep_t + _ep_t * rhs._real;
    T ep_k0_l0 = _real * rhs._ep_k0_l0
        + _ep_k0_l0 * rhs._real
        + _ep_l0 * rhs._ep_k0
        + _ep_k0 * rhs._ep_l0;
    T ep_k0_t = _real * rhs._ep_k0_t
        + _ep_k0_t * rhs._real
        + _ep_k0 * rhs._ep_t
        + _ep_t * rhs._ep_k0;
    T ep_l0_t = _real * rhs._ep_l0_t
        + _ep_l0_t * rhs._real
        + _ep_l0 * rhs._ep_t
        + _ep_t * rhs._ep_l0;
    T ep_t2 = _ep_t * rhs._ep_t + _real * rhs._ep_t2 + _ep_t2 * rhs._real;
    T ep_k0_l0_t = rhs._ep_k0_t * _ep_l0
        + _ep_k0_t * rhs._ep_l0
        + rhs._ep_k0 * _ep_l0_t
        + _ep_k0 * rhs._ep_l0_t
        + rhs._ep_k0_l0_t * _real
        + _ep_k0_l0_t * rhs._real
        + rhs._ep_k0_l0 * _ep_t
        + _ep_k0_l0 * rhs._ep_t;
    T ep_k0_t2 = _real * rhs._ep_k0_t2
        + rhs._real * _ep_k0_t2
        + _ep_t * rhs._ep_k0_t
        + _ep_k0_t * rhs._ep_t
        + _ep_t2 * rhs._ep_k0
        + _ep_k0 * rhs._ep_t2;
    T ep_l0_t2 = _real * rhs._ep_l0_t2
        + rhs._real * _ep_l0_t2
        + _ep_t * rhs._ep_l0_t
        + _ep_l0_t * rhs._ep_t
        + _ep_t2 * rhs._ep_l0
        + _ep_l0 * rhs._ep_t2;
    T ep_t3 = _ep_t * rhs._ep_t2
        + _ep_t2 * rhs._ep_t
        + _real * rhs._ep_t3
        + rhs._real * _ep_t3;

    _real = real;
    _ep_k0 = ep_k0;
    _ep_l0 = ep_l0;
    _ep_t = ep_t;
    _ep_k0_l0 = ep_k0_l0;
    _ep_k0_t = ep_k0_t;
    _ep_l0_t = ep_l0_t;
    _ep_t2 = ep_t2;
    _ep_k0_l0_t = ep_k0_l0_t;
    _ep_k0_t2 = ep_k0_t2;
    _ep_l0_t2 = ep_l0_t2;
    _ep_t3 = ep_t3;
    return *this;
}

  /// Divide this dualklt3 by another dualklt3 of the same or lower depth.
  template<class X, CPPDUALSKLT3_ONLY_SAME_DEPTH_AS_T(T,X)>
  dualklt3<T> & operator/=(const dualklt3<X> & rhs) {
    T r1 = T(1) / rhs._real;
    T r2 = r1 * r1;
    T r3 = r2 * r1;
    T r4 = r3 * r1;
    T real = _real * r1;
    T ep_k0 = _ep_k0 * r1 - _real * rhs._ep_k0 * r2;
    T ep_l0 = _ep_l0 * r1 - _real * rhs._ep_l0 * r2;
    T ep_t = _ep_t * r1 - _real * rhs._ep_t * r2;
    T ep_k0_l0 = _ep_k0_l0 * r1
        - (rhs._ep_k0_l0 * _real + _ep_l0 * rhs._ep_k0 + _ep_k0 * rhs._ep_l0) * r2
        + T(2) * _real * rhs._ep_k0 * rhs._ep_l0 * r3;
    T ep_k0_t = _ep_k0_t * r1
        - (rhs._ep_k0_t * _real + _ep_k0 * rhs._ep_t + _ep_t * rhs._ep_k0) * r2
        + T(2) * _real * rhs._ep_k0 * rhs._ep_t * r3;
    T ep_l0_t = _ep_l0_t * r1
        - (rhs._ep_l0_t * _real + _ep_l0 * rhs._ep_t + _ep_t * rhs._ep_l0) * r2
        + T(2) * _real * rhs._ep_l0 * rhs._ep_t * r3;
    T ep_t2 = _ep_t2 * r1 - (_ep_t * rhs._ep_t + _real * rhs._ep_t2) * r2
        + _real * rhs._ep_t * rhs._ep_t * r3;
    T ep_k0_l0_t = _ep_k0_l0_t * r1
        - r2 * (_ep_k0 * rhs._ep_l0_t
            + rhs._ep_k0 * _ep_l0_t
            + _ep_l0 * rhs._ep_k0_t
            + rhs._ep_l0 * _ep_k0_t
            + _ep_t * rhs._ep_k0_l0
            + rhs._ep_t * _ep_k0_l0
            + _real * rhs._ep_k0_l0_t)
        + r3 * T(2)
            * (_real
                * (rhs._ep_k0_t * rhs._ep_l0
                    + rhs._ep_l0_t * rhs._ep_k0
                    + rhs._ep_k0_l0 * rhs._ep_t)
                + _ep_t * rhs._ep_k0 * rhs._ep_l0
                + _ep_k0 * rhs._ep_l0 * rhs._ep_t
                + _ep_l0 * rhs._ep_k0 * rhs._ep_t)
        - r4 * T(6)
            * _real
            * rhs._ep_k0
            * rhs._ep_l0
            * rhs._ep_t;
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
    T ep_l0_t2 = _ep_l0_t2 * r1
        - r2 * (_real * rhs._ep_l0_t2
            + rhs._ep_l0_t * _ep_t
            + _ep_l0_t * rhs._ep_t
            + _ep_t2 * rhs._ep_l0
            + _ep_l0 * rhs._ep_t2)
        + r3 * (T(2)
            * _real
            * (rhs._ep_l0_t * rhs._ep_t + rhs._ep_t2 * rhs._ep_l0)
            + T(2) * _ep_t * rhs._ep_l0 * rhs._ep_t
            + _ep_l0 * rhs._ep_t * rhs._ep_t)
        - r4 * T(3)
            * _real
            * rhs._ep_l0
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
    _ep_l0 = ep_l0;
    _ep_t = ep_t;
    _ep_k0_l0 = ep_k0_l0;
    _ep_k0_t = ep_k0_t;
    _ep_l0_t = ep_l0_t;
    _ep_t2 = ep_t2;
    _ep_k0_l0_t = ep_k0_l0_t;
    _ep_k0_t2 = ep_k0_t2;
    _ep_l0_t2 = ep_l0_t2;
    _ep_t3 = ep_t3;
    return *this;
}

};

/// Get the dualklt3's real part.
template <class T> T rpart(const dualklt3<T> & x) { return x.rpart(); }

#ifndef PARSED_BY_DOXYGEN

/// Dual +-*/ ops with another entity
#define CPPDUALSKLT3_BINARY_OP(op)                                          \
  template<class T, class U, CPPDUALSKLT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)> \
  common_t                                                              \
  operator op(const dualklt3<T> & z, const dualklt3<U> & w) {                   \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
    CPPDUALSKLT3_ENABLE_IF(!std::is_same<U,std::complex<dualklt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const dualklt3<T> & z, const U & w) {                         \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
    CPPDUALSKLT3_ENABLE_IF(!std::is_same<U,std::complex<dualklt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const U & z, const dualklt3<T> & w) {                         \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U), \
    CPPDUALSKLT3_ENABLE_IF(!std::is_same<U,std::complex<dualklt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const std::complex<dualklt3<T>> & z, const U & w) {           \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_CX_COMMON_T(T,U), \
    CPPDUALSKLT3_ENABLE_IF(!std::is_same<U,std::complex<dualklt3<T>>>::value)>  \
  common_t                                                              \
  operator op(const U & z, const std::complex<dualklt3<T>> & w) {           \
    common_t x(z);                                                      \
    return x op##= w;                                                   \
  }                                                                \

CPPDUALSKLT3_BINARY_OP(+)
CPPDUALSKLT3_BINARY_OP(-)
CPPDUALSKLT3_BINARY_OP(*)
CPPDUALSKLT3_BINARY_OP(/)

/// Dual compared to a non-complex lower rank thing
#define CPPDUALSKLT3_LHS_COMPARISON(op)                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)> \
  bool                                                                  \
  operator op(const dualklt3<T> & a, const dualklt3<U> & b) {                   \
    return a.rpart() op b.rpart();                                      \
  }                                                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
           CPPDUALSKLT3_ENABLE_IF(!is_complex<U>::value)>                   \
  bool                                                                  \
  operator op(const U & a, const dualklt3<T> & b) {                         \
    return a op b.rpart();                                              \
  }                                                                     \
  template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U), \
           CPPDUALSKLT3_ENABLE_IF(!is_complex<U>::value)>                   \
  bool                                                                  \
  operator op(const dualklt3<T> & a, const U & b) {                         \
    return a.rpart() op b;                                              \
  }

CPPDUALSKLT3_LHS_COMPARISON(<)
CPPDUALSKLT3_LHS_COMPARISON(>)
CPPDUALSKLT3_LHS_COMPARISON(<=)
CPPDUALSKLT3_LHS_COMPARISON(>=)
CPPDUALSKLT3_LHS_COMPARISON(==)
CPPDUALSKLT3_LHS_COMPARISON(!=)

#endif // PARSED_BY_DOXYGEN

// NOTE: Added by Ben Ruijl
template<class T> dualklt3<T> operator*(const T & lhs, const dualklt3<T> rhs) { return rhs * lhs; }
template<class T> dualklt3<T> operator/(const T & lhs, const dualklt3<T> rhs) { return rhs.inv() * lhs; }
template<class T> dualklt3<T> operator+(const T & lhs, const dualklt3<T> rhs) { return rhs + lhs; }
template<class T> dualklt3<T> operator-(const T & lhs, const dualklt3<T> rhs) { return -rhs + lhs; }

/// Exponential e^x
template<class T> dualklt3<T> exp(const dualklt3<T> & x) {
  using std::exp;
  T r = exp(x._real);
  return dualklt3<T>(
      r,
      r * x._ep_k0,
      r * x._ep_l0,
      r * x._ep_t,
      r * (x._ep_k0_l0 + x._ep_k0 * x._ep_l0),
      r * (x._ep_k0_t + x._ep_k0 * x._ep_t),
      r * (x._ep_l0_t + x._ep_l0 * x._ep_t),
      r * (x._ep_t * x._ep_t / T(2) + x._ep_t2),
      r * (x._ep_k0_l0_t
              + x._ep_k0_l0 * x._ep_t
              + x._ep_k0_t * x._ep_l0
              + x._ep_l0_t * x._ep_k0
              + x._ep_k0 * x._ep_l0 * x._ep_t),
      r * (x._ep_k0_t2
              + x._ep_k0_t * x._ep_t
              + x._ep_k0
                  * (x._ep_t2 + x._ep_t * x._ep_t / T(2))),
      r * (x._ep_l0_t2
              + x._ep_l0_t * x._ep_t
              + x._ep_l0
                  * (x._ep_t2 + x._ep_t * x._ep_t / T(2))),
      r * (x._ep_t * x._ep_t * x._ep_t / T(6)
              + x._ep_t3
              + x._ep_t2 * x._ep_t)
  );
}

/// Natural log ln(x)
template<class T> dualklt3<T> log(const dualklt3<T> & x) {
  using std::log;
  T v = log(x.rpart());
  assert(false);
}

template<class T> dualklt3<T> log10(const dualklt3<T> & x) {
  using std::log;
  return log(x) / log(static_cast<T>(10));
}

template<class T> dualklt3<T> log2(const dualklt3<T> & x) {
  using std::log;
  return log(x) / log(static_cast<T>(2));
}

template<class T, class U, CPPDUALSKLT3_ENABLE_SAME_DEPTH_AND_COMMON_T(T,U)>
common_t
pow(const dualklt3<T> & f, const dualklt3<U> & g) {
  using std::pow;
  using std::log;
  /*T v = pow(f.rpart(), g.rpart());
  common_t(v,
                  pow(f.rpart(), g.rpart() - T(1)) *
                  (g.rpart() * f.dpart()
                   + f.rpart() * log(f.rpart()) * g.dpart()));*/
  assert(false);
}

template<class T> dualklt3<T> pow(const dualklt3<T> & x, int n) {
  using std::pow;
  T rmmm = pow(x._real, n - 3);
  T rmm = rmmm * x._real;
  T rm = rmm * x._real;
  T r = rm * x._real;
  T dr = T(n) * rm;
  T ddr = T(n) * rmm;
  T dddr = T(n) * rmmm;
  return dualklt3<T>(
    r,
    x._ep_k0 * dr,
    x._ep_l0 * dr,
    x._ep_t * dr,
    x._ep_k0_l0 * dr
        + x._ep_k0 * x._ep_l0 * T(n - 1) * ddr,
    x._ep_k0_t * dr
        + x._ep_k0 * x._ep_t * T(n - 1) * ddr,
    x._ep_l0_t * dr
        + x._ep_l0 * x._ep_t * T(n - 1) * ddr,
    x._ep_t2 * dr
        + x._ep_t * x._ep_t / T(2)
            * T(n-1)
            * ddr,
    x._ep_k0_l0_t * dr
        + ddr
            * (x._ep_k0_l0 * x._ep_t
                + x._ep_k0_t * x._ep_l0
                + x._ep_l0_t * x._ep_k0)
            * T(n - 1)
        + dddr
            * x._ep_k0
            * x._ep_l0
            * x._ep_t
            * T((n - 1) * (n - 2)),
    x._ep_k0_t2 * dr
        + ddr
            * (x._ep_k0 * x._ep_t2 + x._ep_k0_t * x._ep_t)
            * T(n - 1)
        + dddr / T(2)
            * x._ep_k0
            * x._ep_t
            * x._ep_t
            * T((n - 1) * (n - 2)),
    x._ep_l0_t2 * dr
        + ddr
            * (x._ep_l0 * x._ep_t2 + x._ep_l0_t * x._ep_t)
            * T(n - 1)
        + dddr / T(2)
            * x._ep_l0
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


template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
pow(const U & x, const dualklt3<T> & y) {
  return pow(common_t(x), y);
}

namespace utils {
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
}

template<class T> dualklt3<T> abs(const dualklt3<T> & x) {
  using std::abs;
  T s = utils::sgn(x.rpart());
  return dualklt3<T>(abs(x._real), x._ep_k0 * s, x._ep_t * s, x._ep_k0_t * s, x._ep_t2 * s);
}

template<class T> dualklt3<T> fabs(const dualklt3<T> & x) {
  using std::fabs;
  T s = utils::sgn(x.rpart());
  return dualklt3<T>(fabs(x._real), x._ep_k0 * s, x._ep_t * s, x._ep_k0_t * s, x._ep_t2 * s);
}

template<class T> dualsklt3::dualklt3<T> copysign(const dualsklt3::dualklt3<T> & x, const dualsklt3::dualklt3<T> & y) {
  using std::copysign;
  T r = copysign(x.rpart(), y.rpart());
  //return dualsklt3::dualklt3<T>(r, r == x.rpart() ? x.dpart() : -x.dpart());
  assert(false);
}

template<class T> dualsklt3::dualklt3<T> hypot(const dualsklt3::dualklt3<T> & x, const dualsklt3::dualklt3<T> & y) {
  return sqrt(x*x + y*y);
}

template<class T> dualsklt3::dualklt3<T> scalbn(const dualsklt3::dualklt3<T> & arg, int ex) {
  return arg * std::pow((T)2, ex);
}

template<class T> dualsklt3::dualklt3<T> (fmax)(const dualsklt3::dualklt3<T> & x, const dualsklt3::dualklt3<T> & y) {
  return x.rpart() > y.rpart() ? x : y;
}

template<class T> dualsklt3::dualklt3<T> (fmin)(const dualsklt3::dualklt3<T> & x, const dualsklt3::dualklt3<T> & y) {
  return x.rpart() <= y.rpart() ? x : y;
}

template<class T> dualsklt3::dualklt3<T> logb(const dualsklt3::dualklt3<T> & x) {
  return dualsklt3::log2(x);
}

template<class T> int (fpclassify)(const dualsklt3::dualklt3<T> & d) { using std::fpclassify; return (fpclassify)(d.rpart()); }
template<class T> bool (isfinite)(const dualsklt3::dualklt3<T> & d) { using std::isfinite; return (isfinite)(d.rpart()); }
template<class T> bool (isnormal)(const dualsklt3::dualklt3<T> & d) { using std::isnormal; return (isnormal)(d.rpart()); }
template<class T> bool (isinf)(const dualsklt3::dualklt3<T> & d) { using std::isinf; return (isinf)(d.rpart()); }
template<class T> bool (isnan)(const dualsklt3::dualklt3<T> & d) { using std::isnan; return (isnan)(d.rpart()); }
template<class T> bool (signbit)(const dualsklt3::dualklt3<T> & d) { using std::signbit; return (signbit)(d.rpart()); }

template<class T> dualklt3<T> sqrt(const dualklt3<T> & x) {
  using std::sqrt;
  T half = T(0.5);
  T r = sqrt(x._real);
  T ir = T(1) / r * half;
  T irr = ir / x._real * half;
  T irrr = irr / x._real * half;

  return dualklt3<T>(
    r,
    ir * x._ep_k0,
    ir * x._ep_l0,
    ir * x._ep_t,
    ir * x._ep_k0_l0 - irr * x._ep_k0 * x._ep_l0,
    ir * x._ep_k0_t - irr * x._ep_k0 * x._ep_t,
    ir * x._ep_l0_t - irr * x._ep_l0 * x._ep_t,
    ir * x._ep_t2 - irr * half * x._ep_t * x._ep_t,
    ir * x._ep_k0_l0_t
        + irrr * T(3) * x._ep_k0 * x._ep_l0 * x._ep_t
        - irr
            * (x._ep_k0_l0 * x._ep_t
                + x._ep_k0_t * x._ep_l0
                + x._ep_l0_t * x._ep_k0),
    irrr
        * T(3)
        * half
        * x._ep_k0
        * x._ep_t
        * x._ep_t
        + ir * x._ep_k0_t2
        - irr * (x._ep_k0 * x._ep_t2 + x._ep_k0_t * x._ep_t),
    irrr
        * T(3)
        * half
        * x._ep_l0
        * x._ep_t
        * x._ep_t
        + ir * x._ep_l0_t2
        - irr * (x._ep_l0 * x._ep_t2 + x._ep_l0_t * x._ep_t),
    irrr * half * x._ep_t * x._ep_t * x._ep_t - irr * x._ep_t2 * x._ep_t
        + ir * x._ep_t3
  );
}

template<class T> dualklt3<T> cbrt(const dualklt3<T> & x) {
  assert(false);
}

template<class T> dualklt3<T> sin(const dualklt3<T> & x) {
  using std::sin;
  using std::cos;
  T c = cos(x._real);
  T s = sin(x._real);
  
  return dualklt3(
    s,
    c * x._ep_k0,
    c * x._ep_l0,
    c * x._ep_t,
    c * x._ep_k0_l0 - s * x._ep_k0 * x._ep_l0,
    c * x._ep_k0_t - s * x._ep_k0 * x._ep_t,
    c * x._ep_l0_t - s * x._ep_l0 * x._ep_t,
    c * x._ep_t2 - s * x._ep_t * x._ep_t / T(2),
    c * (x._ep_k0_l0_t - x._ep_k0 * x._ep_l0 * x._ep_t)
        - s * (x._ep_k0_t * x._ep_l0
            + x._ep_k0 * x._ep_l0_t
            + x._ep_k0_l0 * x._ep_t),
    c
        * (x._ep_k0_t - x._ep_k0 * x._ep_t * x._ep_t / T(2))
        - s * (x._ep_k0_t * x._ep_t + x._ep_k0 * x._ep_t2),
    c
        * (x._ep_l0_t - x._ep_l0 * x._ep_t * x._ep_t / T(2))
        - s * (x._ep_l0_t * x._ep_t + x._ep_l0 * x._ep_t2),
    c
        * (x._ep_t3
            - x._ep_t * x._ep_t * x._ep_t / T(6))
        - s * (x._ep_t * x._ep_t2)
  );
}

template<class T> dualklt3<T> cos(const dualklt3<T> & x) {
  using std::cos;
  using std::sin;
  T c = cos(x._real);
  T s = sin(x._real);
  return dualklt3(
    c,
    -s * x._ep_k0,
    -s * x._ep_l0,
    -s * x._ep_t,
    -s * x._ep_k0_l0 - c * x._ep_k0 * x._ep_l0,
    -s * x._ep_k0_t - c * x._ep_k0 * x._ep_t,
    -s * x._ep_l0_t - c * x._ep_l0 * x._ep_t,
    -s * x._ep_t2 - c * x._ep_t * x._ep_t / T(2),
    -s * (x._ep_k0_l0_t - x._ep_k0 * x._ep_l0 * x._ep_t)
        - c * (x._ep_k0_t * x._ep_l0
            + x._ep_k0 * x._ep_l0_t
            + x._ep_k0_l0 * x._ep_t),
    -s
        * (x._ep_k0_t - x._ep_k0 * x._ep_t * x._ep_t / T(2))
        - c * (x._ep_k0_t * x._ep_t + x._ep_k0 * x._ep_t2),
    -s
        * (x._ep_l0_t - x._ep_l0 * x._ep_t * x._ep_t / T(2))
        - c * (x._ep_l0_t * x._ep_t + x._ep_l0 * x._ep_t2),
    -s
        * (x._ep_t3
            - x._ep_t * x._ep_t * x._ep_t / T(6))
        - c * (x._ep_t * x._ep_t2)
  );
}

template<class T> dualklt3<T> tan(const dualklt3<T> & x) {
  using std::tan;
  T v = tan(x.rpart());
  assert(false);
}

template<class T> dualklt3<T> asin(const dualklt3<T> & x) {
  using std::asin;
  using std::sqrt;
  T v = asin(x.rpart());
  assert(false);
}

template<class T> dualklt3<T> acos(const dualklt3<T> & x) {
  using std::acos;
  using std::sqrt;
  assert(false);
}

template<class T> dualklt3<T> atan(const dualklt3<T> & x) {
  using std::atan;
  assert(false);
}

template<class T> dualklt3<T> atan2(const dualklt3<T> & y, const dualklt3<T> & x) {
  using std::atan2;
  T v = atan2(y.rpart(), x.rpart());
  assert(false);
}

// TODO
template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
atan2(const dualklt3<T> & y, const U & x) {
  using std::atan2;
  T v = atan2(y.rpart(), x);
  assert(false);
}

// TODO
template<class T, class U, CPPDUALSKLT3_ENABLE_LEQ_DEPTH_AND_COMMON_T(T,U)>
common_t
atan2(const U & y, const dualklt3<T> & x) {
  using std::atan2;
  T v = atan2(y, x.rpart());
  assert(false);
}

// TODO
template<class T> dualklt3<T> sinh(const dualklt3<T> & x);
template<class T> dualklt3<T> cosh(const dualklt3<T> & x);
template<class T> dualklt3<T> tanh(const dualklt3<T> & x);
template<class T> dualklt3<T> asinh(const dualklt3<T> & x);
template<class T> dualklt3<T> acosh(const dualklt3<T> & x);
template<class T> dualklt3<T> atanh(const dualklt3<T> & x);
template<class T> dualklt3<T> log1p(const dualklt3<T> & x);
template<class T> dualklt3<T> expm1(const dualklt3<T> & x);

/// The error function.  Make sure to `#include <math.h>` before
/// `#include <dualklt3s/dualklt3>` to use this function.
template<class T> dualklt3<T> erf(const dualklt3<T> & x) {
  using std::erf;
  using std::sqrt;
  using std::pow;
  using std::exp;
  assert(false);
}

/// Error function complement (1 - erf()).
template<class T> dualklt3<T> erfc(const dualklt3<T> & x) {
  using std::erfc;
  using std::sqrt;
  using std::pow;
  using std::exp;
  assert(false);
}

/// Gamma function.  Approximation of the dualklt3 part.
// TODO specialize for integers
template<class T> dualklt3<T> tgamma(const dualklt3<T> & x) {
  using std::tgamma;
  assert(false);
}

/// Log of absolute value of gamma function.  Approximation of the dualklt3 part.
template<class T> dualklt3<T> lgamma(const dualklt3<T> & x) {
  using std::lgamma;
  T v = lgamma(x.rpart());
  assert(false);
}

/// Putto operator
template<class T, class _CharT, class _Traits>
std::basic_ostream<_CharT, _Traits> &
operator<<(std::basic_ostream<_CharT, _Traits> & os, const dualklt3<T> & x)
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
    << x._ep_l0
    << "_e_l"
    << " + "
    << x._ep_t
    << "_e_t"
    << " + "
    << x._ep_k0_l0
    << "_e_kl"
    << " + "
    << x._ep_k0_t
    << "_e_kt"
    << " + "
    << x._ep_l0_t
    << "_e_lt"
    << " + "
    << x._ep_t2
    << "_e_t2"
    << " + "
    << x._ep_k0_l0_t
    << "_e_klt"
    << " + "
    << x._ep_k0_t2
    << "_e_kt2"
    << " + "
    << x._ep_l0_t2
    << "_e_lt2"
    << " + "
    << x._ep_t3
    << "_e_t3"
    << ")";
  return os << s.str();
}

} // namespace dualsklt3

#ifdef _LIBCPP_BEGIN_NAMESPACE_STD
_LIBCPP_BEGIN_NAMESPACE_STD
#else
namespace std {
#endif

#ifndef PARSED_BY_DOXYGEN

#define make_math(T)                                                    \
  inline T (frexp)(const dualsklt3::dualklt3<T> & arg, int* exp ) { return (frexp)(arg.rpart(), exp); } \
  inline dualsklt3::dualklt3<T> (ldexp)(const dualsklt3::dualklt3<T> & arg, int exp ) { return arg * std::pow((T)2,exp); } \
  inline T (trunc)(const dualsklt3::dualklt3<T> & d) { return (trunc)(d.rpart()); } \
  inline T (floor)(const dualsklt3::dualklt3<T> & d) { return (floor)(d.rpart()); } \
  inline T (ceil)(const dualsklt3::dualklt3<T> & d)  { return (ceil)(d.rpart()); } \
  inline T (round)(const dualsklt3::dualklt3<T> & d) { return (round)(d.rpart()); } \
  inline int (fpclassify)(const dualsklt3::dualklt3<T> & d) { return (fpclassify)(d.rpart()); } \
  inline bool (isfinite)(const dualsklt3::dualklt3<T> & d) { return (isfinite)(d.rpart()); } \
  inline bool (isnormal)(const dualsklt3::dualklt3<T> & d) { return (isnormal)(d.rpart()); } \
  inline bool (isinf)(const dualsklt3::dualklt3<T> & d) { return (isinf)(d.rpart()); } \
  inline bool (isnan)(const dualsklt3::dualklt3<T> & d) { return (isnan)(d.rpart()); } \

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

#endif // DUALKLT3
