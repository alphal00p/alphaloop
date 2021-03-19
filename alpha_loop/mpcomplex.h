//
// Based on the GNU std::complex template library and turned into a specific
// class for handling multi-precision complex numbers.
//

#pragma once

#include "mpreal.h"
#include <sstream>
#include <complex.h>
#include <quadmath.h>

namespace mpfr
{
    // Forward declarations.
    class mpcomplex;

    mpreal abs(const mpcomplex&);
    mpreal arg(const mpcomplex&);
    mpreal norm(const mpcomplex&);

    mpcomplex conj(const mpcomplex&);
    mpcomplex polar(const mpfr::mpreal&, const mpfr::mpreal&);

    // Transcendentals:
    mpcomplex cos(const mpcomplex&);
    mpcomplex cosh(const mpcomplex&);
    mpcomplex exp(const mpcomplex&);
    mpcomplex log(const mpcomplex&);
    mpcomplex log10(const mpcomplex&);
    mpcomplex pow(const mpcomplex&, int);
    mpcomplex pow(const mpcomplex&, const mpreal&);
    mpcomplex pow(const mpcomplex&, const mpcomplex&);
    mpcomplex pow(const mpreal&, const mpcomplex&);
    mpcomplex sin(const mpcomplex&);
    mpcomplex sinh(const mpcomplex&);
    mpcomplex sqrt(const mpcomplex&);
    mpcomplex tan(const mpcomplex&);
    mpcomplex tanh(const mpcomplex&);  
   
    class mpcomplex
    {
    public:
        mpcomplex(const mpfr::mpreal& r= mpreal(), const mpreal& i = mpfr::mpreal())
        : m_real(r), m_imag(i) { }

        mpcomplex(const double& r)
        : m_real(mpfr::mpreal(r)), m_imag(mpfr::mpreal()) { }

        mpcomplex(const __complex128& z) : m_real(crealq(z)), m_imag(cimagq(z)) { }
        mpcomplex(const _Complex& z) : m_real(creal(z)), m_imag(cimag(z)) { }

        explicit operator __complex128 () const { return (__float128)m_real + (__float128)m_imag*I; }

        mpreal& real() { return m_real; }
        const mpreal& real() const { return m_real; }

        ///  Return imaginary part of complex number.
        mpreal& imag() { return m_imag; }
        const mpreal& imag() const { return m_imag; }

        /// Assign this complex number to scalar @a t.
        mpcomplex& operator=(const mpreal&);

        mpcomplex&
        operator+=(const mpreal& t)
        {
            m_real += t;
            return *this;
        }

        mpcomplex&
        operator-=(const mpreal& t)
        {
            m_real -= t;
            return *this;
        }

        mpcomplex& operator*=(const mpreal&);
        mpcomplex& operator/=(const mpreal&);

        // Lets the compiler synthesize the
        // copy and assignment operator
        // complex<_Tp>& operator= (const complex<_Tp>&);
        /// Assign this complex number to complex @a z.
        mpcomplex& operator=(const mpcomplex&);
        /// Add @a z to this complex number.
        mpcomplex& operator+=(const mpcomplex&);
        /// Subtract @a z from this complex number.
        mpcomplex& operator-=(const mpcomplex&);
        /// Multiply this complex number by @a z.
        mpcomplex& operator*=(const mpcomplex&);
        /// Divide this complex number by @a z.
        mpcomplex& operator/=(const mpcomplex&);

    private:
        mpreal m_real;
        mpreal m_imag;
    };

    
    inline mpcomplex operator+(const mpcomplex& x, const mpcomplex& y)
    {
        mpcomplex r = x;
        r += y;
        return r;
    }

    inline mpcomplex operator+(const mpcomplex& x, const mpreal& y)
    {
        mpcomplex r = x;
        r += y;
        return r;
    }

    inline mpcomplex operator+(const mpreal& x, const mpcomplex& y)
    {
        mpcomplex r = y;
        r += x;
        return r;
    }

    inline mpcomplex operator-(const mpcomplex& x, const mpcomplex& y)
    {
        mpcomplex r = x;
        r -= y;
        return r;
    }
    
    inline mpcomplex operator-(const mpcomplex& x, const mpreal& y)
    {
        mpcomplex r = x;
        r -= y;
        return r;
    }

    inline mpcomplex operator-(const mpreal& x, const mpcomplex& y)
    {
        mpcomplex r(x, - y.imag());
        r -= y.real();
        return r;
    }

    inline mpcomplex operator*(const mpcomplex& x, const mpcomplex& y)
    {
        mpcomplex r = x;
        r *= y;
        return r;
    }

    inline mpcomplex operator*(const mpcomplex& x, const mpreal& y)
    {
        mpcomplex r = x;
        r *= y;
        return r;
    }

    inline mpcomplex operator*(const mpreal& x, const mpcomplex& y)
    {
        mpcomplex r = y;
        r *= x;
        return r;
    }

    inline mpcomplex operator/(const mpcomplex& x, const mpcomplex& y)
    {
        mpcomplex r = x;
        r /= y;
        return r;
    }
    
    inline mpcomplex operator/(const mpcomplex& x, const mpreal& y)
    {
        mpcomplex r = x;
        r /= y;
        return r;
    }

    inline mpcomplex operator/(const mpreal& x, const mpcomplex& y)
    {
        mpcomplex r = x;
        r /= y;
        return r;
    }

    inline mpcomplex operator+(const mpcomplex& x)
    { return x; }

    inline mpcomplex operator-(const mpcomplex& x)
    { return mpcomplex(-x.real(), -x.imag()); }

    inline bool operator==(const mpcomplex& x, const mpcomplex& y)
    { return x.real() == y.real() &&x.imag() == y.imag(); }

    inline bool operator==(const mpcomplex& x, const mpreal& y)
    { return x.real() == y && x.imag() == mpreal(); }

    inline bool operator==(const mpreal& x, const mpcomplex& y)
    { return x == y.real() && mpreal(0, y.imag().getPrecision()) == y.imag(); }

    inline bool operator!=(const mpcomplex& x, const mpcomplex& y)
    { return x.real() != y.real() || x.imag() != y.imag(); }

    inline bool operator!=(const mpcomplex& x, const mpreal& y)
    { return x.real() != y || x.imag() != mpreal(0, x.imag().getPrecision()); }

    inline bool operator!=(const mpreal& x, const mpcomplex& y)
    { return x != y.real() || mpreal(0, y.imag().getPrecision()) != y.imag(); }

    // Values
    inline mpreal& real(mpcomplex& z) { return z.real(); }
    inline const mpreal& real(const mpcomplex& z) { return z.real(); }
    inline mpreal& imag(mpcomplex& z) { return z.imag(); }
    inline const mpreal& imag(const mpcomplex& z) { return z.imag(); }

    inline mpreal abs(const mpcomplex& z)
    {
        if (mpfr::isinf(mpfr::abs(z.real())) || mpfr::isinf(mpfr::abs(z.imag())))
        {
            auto inf = mpreal(0, z.real().getPrecision());
            inf.setInf(1);
            return inf;
        }
        else
        {
            mpreal x = z.real();
            mpreal y = z.imag();
            const mpreal s = mpfr::max(mpfr::abs(x), mpfr::abs(y));
            if (s == mpreal(0, x.getPrecision())) return s;
            x /= s; 
            y /= s;
            return s*mpfr::sqrt(x*x + y*y);
        }
    }

    inline mpreal arg(const mpcomplex& z) { return mpfr::atan2(z.imag(), z.real()); }

    inline mpreal norm(const mpcomplex& z)
    {
        if (mpfr::isinf(mpfr::abs(z.real())) || mpfr::isinf(mpfr::abs(z.imag())))
        {
            auto inf = mpreal(0, z.real().getPrecision());
            inf.setInf(1);
            return inf;
        }
        else
        {
            const mpreal x = z.real();
            const mpreal y = z.imag();
            return x*x + y*y;
        }
    }

    inline mpcomplex polar(const mpreal& rho, const mpreal& theta)
    {
        return mpcomplex(rho*cos(theta), rho*sin(theta));
    }

    inline mpcomplex conj(const mpcomplex& z)
    {
        return mpcomplex(z.real(), -z.imag());
    }
  
    inline mpcomplex cos(const mpcomplex& z)
    {
        const mpreal x = z.real();
        const mpreal y = z.imag();
        return mpcomplex(mpfr::cos(x)*mpfr::cosh(y), -mpfr::sin(x)*mpfr::sinh(y));
    }

    inline mpcomplex cosh(const mpcomplex& z)
    {
        const mpreal x = z.real();
        const mpreal y = z.imag();
        return mpcomplex(mpfr::cosh(x)*mpfr::cos(y), mpfr::sinh(x)*mpfr::sin(y));
    }

    inline mpcomplex exp(const mpcomplex& z)
    {
        return polar(mpfr::exp(z.real()), z.imag());
    }

    inline mpcomplex log(const mpcomplex& z)
    {
        return mpcomplex(mpfr::log(mpfr::abs(z)), arg(z));
    }

    inline mpcomplex log10(const mpcomplex& z)
    { 
        return log(z) / log(mpreal("10", z.real().getPrecision()));
    }

    inline mpcomplex sin(const mpcomplex& z)
    {
        const mpreal x = z.real();
        const mpreal y = z.imag();
        return mpcomplex(mpfr::sin(x)*mpfr::cosh(y), mpfr::cos(x)*sinh(y)); 
    }

    inline mpcomplex sinh(const mpcomplex& z)
    {
        const mpreal x = z.real();
        const mpreal y = z.imag();
        return mpcomplex(mpfr::sinh(x)*mpfr::cos(y), mpfr::cosh(x)*mpfr::sin(y));
    }

    inline mpcomplex tan(const mpcomplex& z)
    {
        return sin(z)/cos(z);
    }

    inline mpcomplex tanh(const mpcomplex& z)
    {
        return sinh(z)/cosh(z);
    }

    inline mpcomplex pow(const mpcomplex& z, int n)
    {
        mpcomplex x = z;
        mpcomplex one(mpreal(1, z.real().getPrecision()));
        bool inverse = n < 0;
        if (inverse) n = -n;
        mpcomplex y = n % 2 ? x : one;
        while (n >>= 1)
        {
            x *= x;
            if (n % 2)
            {
                y *= x;
            }
        }
        return inverse ? one/y : y;
    }

    inline mpcomplex pow(const mpcomplex& x, const mpreal& y)
    {
        auto zero = mpreal(0, x.real().getPrecision());
        if (x == mpcomplex(zero, zero))
        {
            if (y == zero)
            {
                auto a = zero;
                auto b = zero;
                a.setNan();
                b.setNan();
                return mpcomplex(a, b);
            }
            else
            {
                return mpcomplex(zero, zero);
            }
        }
        else if ((x.imag() == zero) && (y < zero))
        {
            return mpcomplex(mpfr::pow(x.real(),y), zero);
        }
        else if (x.imag() == zero)
        {
            auto a = zero;
            auto b = zero;
            a.setInf(1);
            b.setNan();
            return mpcomplex(a, b);
        }
        else if ((x.imag() == zero) && (y >= zero))
        {
            return mpcomplex(zero, zero);
        }
        else if ((x.imag() == zero) && (x.real() > zero))
        {
            return mpcomplex(mpfr::pow(x.real(),y), zero);
        }
        mpcomplex t = log(x);
        return polar(mpfr::exp(y*t.real()), y*t.imag());
    }

    inline mpcomplex pow(const mpcomplex& x, const mpcomplex& y)
    {
        auto zero = mpreal(0, x.real().getPrecision());
        if (x == mpcomplex(zero, zero))
        {
            if (y == mpcomplex(zero, zero) || ((y.real() == zero) && (y.imag() != zero)))
            {
                auto a = zero;
                auto b = zero;
                a.setNan();
                b.setNan();
                return mpcomplex(a, b);
            }
            else if (y.real() > zero)
            {
                return mpcomplex(zero, zero);
            }
            else if (y.real() < zero)
            {
                auto a = zero;
                auto b = zero;
                a.setInf(1);
                b.setNan();
                return mpcomplex(a, b);
            }
        }
        else if (y == mpcomplex(zero, zero))
        {
            if (mpfr::isinf(mpfr::abs(x.real())) || mpfr::isinf(mpfr::abs(x.imag())))
            {
                auto a = zero;
                auto b = zero;
                a.setNan();
                b.setNan();
                return mpcomplex(a, b);
            }
            return mpcomplex(zero, zero);
        }
        else if (mpfr::isinf(mpfr::abs(x.real())) || mpfr::isinf(mpfr::abs(x.imag())))
        {
            if (y.real() > zero)
            {
                return x;
            }
            else if (y.real() < zero)
            {
                return mpcomplex(zero, zero);
            }
            else
            {
                auto a = zero;
                auto b = zero;
                a.setNan();
                b.setNan();
                return mpcomplex(a, b);
            }
        }
        return exp(y*log(x));
    }

    inline mpcomplex pow(const mpreal& x, const mpcomplex& y)
    {
      return x > mpreal(0, x.getPrecision()) ? polar(pow(x, y.real()), y.imag()*log(x))
	                                     : pow(mpcomplex(x), y);
    }
    
    inline mpcomplex& mpcomplex::operator=(const mpreal& t)
    {
        m_real = t;
        m_imag = mpreal();
        return *this;
    } 

    inline mpcomplex& mpcomplex::operator*=(const mpreal& t)
    {
        m_real *= t;
        m_imag *= t;
        return *this;
    }

    inline mpcomplex& mpcomplex::operator/=(const mpreal& t)
    {
        m_real /= t;
        m_imag /= t;
        return *this;
    }

    inline mpcomplex& mpcomplex::operator=(const mpcomplex& z)
    {
        m_real = z.real();
        m_imag = z.imag();
        return *this;
    }

    inline mpcomplex& mpcomplex::operator+=(const mpcomplex& z)
    {
        m_real += z.real();
        m_imag += z.imag();
        return *this;
    }

    inline mpcomplex& mpcomplex::operator-=(const mpcomplex& z)
    {
        m_real -= z.real();
        m_imag -= z.imag();
        return *this;
    }

    inline mpcomplex& mpcomplex::operator*=(const mpcomplex& z)
    {
        auto zero = mpreal(0, z.real().getPrecision());
        if (mpfr::isinf(mpfr::abs(m_real)) || mpfr::isinf(mpfr::abs(m_imag)))
        {
            if (z == mpcomplex(zero, zero))
            {
                m_real.setNan();
                m_imag.setNan();
            }
            return *this;
        }
        else if (mpfr::isinf(mpfr::abs(z.real())) || mpfr::isinf(mpfr::abs(z.imag())))
        {
            if ((m_real == zero) && (m_imag == zero))
            {
                m_real.setNan();
                m_imag.setNan();
            }
            else
            {
                m_real = z.real();
                m_imag = z.imag();
            }
            return *this;
        }
        const mpreal r = m_real*z.real() - m_imag*z.imag();
        m_imag = m_real*z.imag() + m_imag*z.real();
        m_real = r;
        return *this;
    }

    inline mpcomplex& mpcomplex::operator/=(const mpcomplex& z)
    {
        auto zero = mpreal(0, z.real().getPrecision());
        if (z == mpcomplex(zero, zero))
        {
            if (m_real == zero)
            {
                m_real.setNan();
            }
            else if (m_real > zero)
            {
                m_real.setInf(1);
            }
            else
            {
                m_real.setInf(-1);
            }
            if (m_imag == zero)
            {
                m_imag.setNan();
            }
            else if (m_imag > zero)
            {
                m_imag.setInf(1);
            }
            else
            {
                m_imag.setInf(-1);
            }
        }
        else if (mpfr::isinf(mpfr::abs(z.real())) || mpfr::isinf(mpfr::abs(z.imag())))
        {
            if (mpfr::isinf(m_real) || mpfr::isinf(m_imag))
            {
                m_real.setNan();
                m_imag.setNan();
            }
            else if (mpfr::isnan(m_real) || mpfr::isnan(m_imag))
            {
                return *this;
            }
            else
            {
                m_real = zero;
                m_imag = zero;
            }
        }
        else
        {
            const mpreal r =  m_real*z.real() + m_imag*z.imag();
            const mpreal n = mpfr::norm(z);
            m_imag = (m_imag*z.real() - m_real*z.imag())/n;
            m_real = r/n;
        }
        return *this;
    }

    inline mpcomplex sqrt(const mpcomplex& z)
    {
        if (mpfr::isinf(mpfr::abs(z.real())) || mpfr::isinf(mpfr::abs(z.imag())))
        {
            return z;
        }
        mpreal x = z.real();
        mpreal y = z.imag();
        mpfr::mpreal two(2, z.real().getPrecision());

        if (x == mpfr::mpreal())
        {
            mpreal t = mpfr::sqrt(mpfr::abs(y) / two);
            return mpfr::mpcomplex(t, y < mpfr::mpreal() ? -t : t);
        }
        else
        {
            mpreal t = mpfr::sqrt(two * (mpfr::abs(z) + mpfr::abs(x)));
            mpreal u = t / two;
            return x > mpfr::mpreal()
                   ? mpcomplex(u, y/t)
                   : mpcomplex(mpfr::abs(y) / t, y < mpreal() ? -u : u);
        }
    }

    /*std::ostream& operator<<(std::ostream& os, const mpcomplex& z)
    {
      os << '(' << z.real() << ',' << z.imag() << ')';
      return os;
    }*/
}

