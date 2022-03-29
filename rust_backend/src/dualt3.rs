extern crate num_traits;

use lorentz_vector::{Field, RealNumberLike};
pub use num_traits::{Float, FloatConst, Num, One, Zero};
use num_traits::{
    FromPrimitive, Inv, MulAdd, MulAddAssign, NumCast, Pow, Signed, ToPrimitive, Unsigned,
};
use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter, LowerExp, Result as FmtResult};
use std::iter::Sum;
use std::num::FpCategory;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

#[repr(C)]
#[derive(Copy, Clone, Hash, Debug, Default)]
pub struct Dualt3<T> {
    pub real: T,
    pub ep_t: T,
    pub ep_t2: T,
    pub ep_t3: T,
}

impl<T> Dualt3<T> {
    pub fn new(real: T, ep_t: T, ep_t2: T, ep_t3: T) -> Dualt3<T> {
        Dualt3 {
            real,
            ep_t,
            ep_t2,
            ep_t3,
        }
    }

    #[inline]
    pub fn from_real(real: T) -> Dualt3<T>
    where
        T: Zero,
    {
        Dualt3 {
            real,
            ep_t: T::zero(),
            ep_t2: T::zero(),
            ep_t3: T::zero(),
        }
    }

    #[inline]
    pub fn map<F, U>(&self, f: F) -> Dualt3<U>
    where
        F: Fn(&T) -> U,
    {
        Dualt3 {
            real: f(&self.real),
            ep_t: f(&self.ep_t),
            ep_t2: f(&self.ep_t2),
            ep_t3: f(&self.ep_t3),
        }
    }

    #[inline]
    pub fn zip_map<F, U, K>(&self, rhs: &Dualt3<U>, f: F) -> Dualt3<K>
    where
        F: Fn(&T, &U) -> K,
    {
        Dualt3 {
            real: f(&self.real, &rhs.real),
            ep_t: f(&self.ep_t, &rhs.ep_t),
            ep_t2: f(&self.ep_t2, &rhs.ep_t2),
            ep_t3: f(&self.ep_t3, &rhs.ep_t3),
        }
    }
}

impl<T: Copy> Dualt3<T> {
    fn substitute_real(&self, new_real: T) -> Dualt3<T> {
        Dualt3 {
            real: new_real,
            ep_t: self.ep_t,
            ep_t2: self.ep_t2,
            ep_t3: self.ep_t3,
        }
    }
}

impl<T: Zero> From<T> for Dualt3<T> {
    #[inline]
    fn from(real: T) -> Dualt3<T> {
        Dualt3::from_real(real)
    }
}

impl<T: Display> Display for Dualt3<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let precision = f.precision().unwrap_or(4);
        write!(
            f,
            "({:.p$} + {:.p$}\u{03B5}_t + {:.p$}\u{03B5}_t^2 + {:.p$}\u{03B5}_t^3)",
            self.real,
            self.ep_t,
            self.ep_t2,
            self.ep_t3,
            p = precision
        )
    }
}

impl<T: Display + LowerExp> LowerExp for Dualt3<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let precision = f.precision().unwrap_or(4);
        write!(
            f,
            "({:.p$e} + {:.p$e}\u{03B5}_t + {:.p$e}\u{03B5}_t^2 + {:.p$e}\u{03B5}_t^3)",
            self.real,
            self.ep_t,
            self.ep_t2,
            self.ep_t3,
            p = precision
        )
    }
}

impl<T: PartialEq> PartialEq<Self> for Dualt3<T> {
    #[inline]
    fn eq(&self, rhs: &Self) -> bool {
        self.real == rhs.real
    }
}

impl<T: PartialOrd> PartialOrd<Self> for Dualt3<T> {
    #[inline]
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        PartialOrd::partial_cmp(&self.real, &rhs.real)
    }

    #[inline]
    fn lt(&self, rhs: &Self) -> bool {
        self.real < rhs.real
    }

    #[inline]
    fn gt(&self, rhs: &Self) -> bool {
        self.real > rhs.real
    }
}

impl<T: PartialEq> PartialEq<T> for Dualt3<T> {
    #[inline]
    fn eq(&self, rhs: &T) -> bool {
        self.real == *rhs
    }
}

impl<T: PartialOrd> PartialOrd<T> for Dualt3<T> {
    #[inline]
    fn partial_cmp(&self, rhs: &T) -> Option<Ordering> {
        PartialOrd::partial_cmp(&self.real, rhs)
    }

    #[inline]
    fn lt(&self, rhs: &T) -> bool {
        self.real < *rhs
    }

    #[inline]
    fn gt(&self, rhs: &T) -> bool {
        self.real > *rhs
    }
}

macro_rules! impl_to_primitive {
    ($($name:ident, $ty:ty),*) => {
        impl<T: ToPrimitive> ToPrimitive for Dualt3<T>
            {
            $(
                #[inline]
                fn $name(&self) -> Option<$ty> {
                    self.real.$name()
                }
            )*
        }
    }
}

macro_rules! impl_from_primitive {
    ($($name:ident, $ty:ty),*) => {
        impl<T: FromPrimitive + Zero> FromPrimitive for Dualt3<T>
           {
            $(
                #[inline]
                fn $name(n: $ty) -> Option<Dualt3<T>> {
                    T::$name(n).map(Dualt3::from_real)
                }
            )*
        }
    }
}

macro_rules! impl_primitive_cast {
    ($($to:ident, $from:ident - $ty:ty),*) => {
        impl_to_primitive!($($to, $ty),*);
        impl_from_primitive!($($from, $ty),*);
    }
}

impl_primitive_cast! {
    to_isize,   from_isize  - isize,
    to_i8,      from_i8     - i8,
    to_i16,     from_i16    - i16,
    to_i32,     from_i32    - i32,
    to_i64,     from_i64    - i64,
    to_usize,   from_usize  - usize,
    to_u8,      from_u8     - u8,
    to_u16,     from_u16    - u16,
    to_u32,     from_u32    - u32,
    to_u64,     from_u64    - u64,
    to_f32,     from_f32    - f32,
    to_f64,     from_f64    - f64
}

impl<T: Num + Copy> Add<T> for Dualt3<T> {
    type Output = Dualt3<T>;

    #[inline]
    fn add(self, rhs: T) -> Dualt3<T> {
        self.substitute_real(self.real + rhs)
    }
}

impl<T: Num + Copy> AddAssign<T> for Dualt3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        *self = (*self) + Dualt3::from_real(rhs)
    }
}

impl<T: Num + Copy> Sub<T> for Dualt3<T> {
    type Output = Dualt3<T>;

    #[inline]
    fn sub(self, rhs: T) -> Dualt3<T> {
        self.substitute_real(self.real - rhs)
    }
}

impl<T: Num + Copy> SubAssign<T> for Dualt3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: T) {
        *self = (*self) - Dualt3::from_real(rhs)
    }
}

impl<T: Num + Copy> Mul<T> for Dualt3<T> {
    type Output = Dualt3<T>;

    #[inline]
    fn mul(self, rhs: T) -> Dualt3<T> {
        self.map(|x| *x * rhs)
    }
}

impl<T: Num + Copy> MulAssign<T> for Dualt3<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = (*self) * Dualt3::from_real(rhs)
    }
}

impl<T: Num + Copy> Div<T> for Dualt3<T> {
    type Output = Dualt3<T>;

    #[inline]
    fn div(self, rhs: T) -> Dualt3<T> {
        self / Dualt3::from_real(rhs)
    }
}

impl<T: Num + Copy> DivAssign<T> for Dualt3<T> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = (*self) / Dualt3::from_real(rhs)
    }
}

impl<T: Signed + Copy> Neg for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        self.map(|x| x.neg())
    }
}

impl<T: Num + Copy> Add<Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        self.zip_map(&rhs, |x, y| *x + *y)
    }
}

impl<T: Num + Copy> Add<&Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: &Self) -> Self {
        self.zip_map(&rhs, |x, y| *x + *y)
    }
}

impl<T: Num + Copy> AddAssign<Self> for Dualt3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = (*self) + rhs
    }
}

impl<T: Num + Copy> Sub<Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        self.zip_map(&rhs, |x, y| *x - *y)
    }
}

impl<T: Num + Copy> Sub<&Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: &Self) -> Self {
        self.zip_map(&rhs, |x, y| *x - *y)
    }
}

impl<T: Num + Copy> SubAssign<Self> for Dualt3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = (*self) - rhs
    }
}

impl<T: Num + Copy> Mul<&Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Self) -> Self {
        Dualt3 {
            real: self.real * rhs.real,
            ep_t: self.real * rhs.ep_t + self.ep_t * rhs.real,
            ep_t2: self.ep_t * rhs.ep_t + self.real * rhs.ep_t2 + self.ep_t2 * rhs.real,
            ep_t3: self.ep_t * rhs.ep_t2
                + self.ep_t2 * rhs.ep_t
                + self.real * rhs.ep_t3
                + rhs.real * self.ep_t3,
        }
    }
}

impl<T: Num + Copy> Mul<Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        self * (&rhs)
    }
}

impl<T: Num + Copy> MulAssign<Self> for Dualt3<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = (*self) * rhs
    }
}

macro_rules! impl_mul_add {
    ($(<$a:ident, $b:ident>),*) => {
        $(
            impl<T: Num + Copy + Mul + Add> MulAdd<$a, $b> for Dualt3<T> {
                type Output = Dualt3<T>;

                #[inline]
                fn mul_add(self, a: $a, b: $b) -> Dualt3<T> {
                    (self * a) + b
                }
            }

            impl<T: Num + Copy + Mul + Add> MulAddAssign<$a, $b> for Dualt3<T> {
                #[inline]
                fn mul_add_assign(&mut self, a: $a, b: $b) {
                    *self = (*self * a) + b;
                }
            }
        )*
    }
}

impl_mul_add! {
    <Self, Self>,
    <T, Self>,
    <Self, T>,
    <T, T>
}

impl<T: Num + Copy> Div<&Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: &Self) -> Self {
        let r1 = T::one() / rhs.real;
        let r2 = r1 * r1;
        let r3 = r2 * r1;
        let r4 = r3 * r1;

        Dualt3 {
            real: self.real * r1,
            ep_t: self.ep_t * r1 - self.real * rhs.ep_t * r2,
            ep_t2: self.ep_t2 * r1 - (self.ep_t * rhs.ep_t + self.real * rhs.ep_t2) * r2
                + self.real * rhs.ep_t * rhs.ep_t * r3,
            ep_t3: self.ep_t3 * r1
                - (rhs.ep_t * self.ep_t2 + self.ep_t * rhs.ep_t2 + self.real * rhs.ep_t3) * r2
                + (self.ep_t * rhs.ep_t * rhs.ep_t
                    + (T::one() + T::one()) * self.real * rhs.ep_t * rhs.ep_t2)
                    * r3
                - self.real * rhs.ep_t * rhs.ep_t * rhs.ep_t * r4,
        }
    }
}

impl<T: Num + Copy> Div<Self> for Dualt3<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        self / (&rhs)
    }
}

impl<T: Num + Copy> DivAssign<Self> for Dualt3<T> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = (*self) / rhs
    }
}

impl<T: Num + Copy> Rem<Self> for Dualt3<T> {
    type Output = Self;

    fn rem(self, _: Self) -> Self {
        unimplemented!()
    }
}

impl<T: Num + Copy> Rem<&Self> for Dualt3<T> {
    type Output = Self;

    fn rem(self, _: &Self) -> Self {
        unimplemented!()
    }
}

impl<T: Num + Copy> RemAssign<Self> for Dualt3<T> {
    fn rem_assign(&mut self, _: Self) {
        unimplemented!()
    }
}

impl<T, P: Into<Dualt3<T>>> Pow<P> for Dualt3<T>
where
    Dualt3<T>: Float,
{
    type Output = Dualt3<T>;

    #[inline]
    fn pow(self, rhs: P) -> Dualt3<T> {
        self.powf(rhs.into())
    }
}

impl<T> Inv for Dualt3<T>
where
    Self: One + Div<Output = Self>,
{
    type Output = Dualt3<T>;

    #[inline]
    fn inv(self) -> Dualt3<T> {
        Dualt3::one() / self
    }
}

impl<T> Signed for Dualt3<T>
where
    T: Signed + PartialOrd + Copy,
{
    #[inline]
    fn abs(&self) -> Self {
        let s = self.real.signum();
        self.map(|x| s * *x)
    }

    #[inline]
    fn abs_sub(&self, rhs: &Self) -> Self {
        if self.real > rhs.real {
            self.sub(*rhs)
        } else {
            Self::zero()
        }
    }

    #[inline]
    fn signum(&self) -> Self {
        Dualt3::from_real(self.real.signum())
    }

    #[inline]
    fn is_positive(&self) -> bool {
        self.real.is_positive()
    }

    #[inline]
    fn is_negative(&self) -> bool {
        self.real.is_negative()
    }
}

impl<T: Unsigned + Copy> Unsigned for Dualt3<T> {}

impl<T: Num + Zero + Copy> Zero for Dualt3<T> {
    #[inline]
    fn zero() -> Dualt3<T> {
        Dualt3::from_real(T::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.real.is_zero()
    }
}

impl<T: Num + One + Copy> One for Dualt3<T> {
    #[inline]
    fn one() -> Dualt3<T> {
        Dualt3::from_real(T::one())
    }

    #[inline]
    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        self.real.is_one()
    }
}

impl<T: Num + Copy> Num for Dualt3<T> {
    type FromStrRadixErr = <T as Num>::FromStrRadixErr;

    #[inline]
    fn from_str_radix(str: &str, radix: u32) -> Result<Dualt3<T>, Self::FromStrRadixErr> {
        <T as Num>::from_str_radix(str, radix).map(Dualt3::from_real)
    }
}

impl<T: Float> NumCast for Dualt3<T> {
    #[inline]
    fn from<P: ToPrimitive>(n: P) -> Option<Dualt3<T>> {
        <T as NumCast>::from(n).map(Dualt3::from_real)
    }
}

macro_rules! impl_float_const {
    ($($c:ident),*) => {
        $(
            fn $c() -> Dualt3<T> { Dualt3::from_real(T::$c()) }
        )*
    }
}

impl<T: FloatConst + Zero> FloatConst for Dualt3<T> {
    impl_float_const!(
        E,
        FRAC_1_PI,
        FRAC_1_SQRT_2,
        FRAC_2_PI,
        FRAC_2_SQRT_PI,
        FRAC_PI_2,
        FRAC_PI_3,
        FRAC_PI_4,
        FRAC_PI_6,
        FRAC_PI_8,
        LN_10,
        LN_2,
        LOG10_E,
        LOG2_E,
        PI,
        SQRT_2
    );
}

macro_rules! impl_real_constant {
    ($($prop:ident),*) => {
        $(
            fn $prop() -> Self { Dualt3::from_real(<T as Float>::$prop()) }
        )*
    }
}

macro_rules! impl_single_boolean_op {
    ($op:ident REAL) => {
        fn $op(self) -> bool {
            self.real.$op()
        }
    };
    ($op:ident OR) => {
        fn $op(self) -> bool {
            unimplemented!();
            /*let mut b = self.real.$op();
            // TODO: propagate to other components
            for x in self.iter().skip(1) {
                b |= x.$op();
            }
            b*/
        }
    };
    ($op:ident AND) => {
        fn $op(self) -> bool {
            unimplemented!();
            /*let mut b = self.real.$op();
            for x in self.iter().skip(1) {
                b &= x.$op();
            }
            b*/
        }
    };
}

macro_rules! impl_boolean_op {
    ($($op:ident $t:tt),*) => {
        $(impl_single_boolean_op!($op $t);)*
    };
}

macro_rules! impl_real_op {
    ($($op:ident),*) => {
        $(
            fn $op(self) -> Self { Dualt3::from_real(self.real.$op()) }
        )*
    }
}

impl<T: Signed + Float + Debug> Float for Dualt3<T> {
    impl_real_constant!(
        nan,
        infinity,
        neg_infinity,
        neg_zero,
        min_positive_value,
        epsilon,
        min_value,
        max_value
    );

    impl_boolean_op!(
        is_nan              OR,
        is_infinite         OR,
        is_finite           AND,
        is_normal           AND,
        is_sign_positive    REAL,
        is_sign_negative    REAL
    );

    #[inline]
    fn classify(self) -> FpCategory {
        self.real.classify()
    }

    impl_real_op!(floor, ceil, round, trunc);

    #[inline]
    fn fract(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn signum(self) -> Self {
        Dualt3::from_real(self.real.signum())
    }

    #[inline]
    fn abs(self) -> Self {
        let s = self.real.signum();
        self.map(|x| *x * s)
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        if self.real > other.real {
            self
        } else {
            other
        }
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        if self.real < other.real {
            self
        } else {
            other
        }
    }

    #[inline]
    fn abs_sub(self, _rhs: Self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn mul_add(self, _a: Self, _b: Self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn recip(self) -> Self {
        Self::one() / self
    }

    #[inline]
    fn powi(self, n: i32) -> Self {
        let (r, dr, ddr, dddr) = if self.real.is_zero() && n >= 0 && n < 3 {
            (
                if n == 0 { T::one() } else { T::zero() },
                if n == 1 { T::one() } else { T::zero() },
                if n == 2 {
                    <T as NumCast>::from(n).unwrap()
                } else {
                    T::zero()
                },
                T::zero(),
            )
        } else {
            let rmmm = self.real.powi(n - 3);
            let rmm = rmmm * self.real;
            let rm = rmm * self.real;
            let r = rm * self.real;
            (
                r,
                <T as NumCast>::from(n).unwrap() * rm,
                <T as NumCast>::from(n).unwrap() * rmm,
                <T as NumCast>::from(n).unwrap() * rmmm,
            )
        };

        Dualt3 {
            real: r,
            ep_t: self.ep_t * dr,
            ep_t2: self.ep_t2 * dr
                + self.ep_t * self.ep_t / <T as NumCast>::from(2).unwrap()
                    * <T as NumCast>::from(n - 1).unwrap()
                    * ddr,
            ep_t3: self.ep_t3 * dr
                + self.ep_t * self.ep_t2 * <T as NumCast>::from(n - 1).unwrap() * ddr
                + self.ep_t * self.ep_t * self.ep_t / <T as NumCast>::from(6).unwrap()
                    * <T as NumCast>::from((n - 1) * (n - 2)).unwrap()
                    * dddr,
        }
    }

    #[inline]
    fn powf(self, _n: Self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn exp(self) -> Self {
        let r = self.real.exp();
        Dualt3 {
            real: r,
            ep_t: r * self.ep_t,
            ep_t2: r * (self.ep_t * self.ep_t / <T as NumCast>::from(2).unwrap() + self.ep_t2),
            ep_t3: r
                * (self.ep_t * self.ep_t * self.ep_t / <T as NumCast>::from(6).unwrap()
                    + self.ep_t3
                    + self.ep_t2 * self.ep_t),
        }
    }

    #[inline]
    fn exp2(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn ln(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn log(self, _base: Self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn log2(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn log10(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn sqrt(self) -> Self {
        let half = T::one() / (T::one() + T::one());
        let r = self.real.sqrt();
        let ir = T::one() / r * half;
        let irr = ir / self.real * half;
        let irrr = irr / self.real * half;
        Dualt3 {
            real: r,
            ep_t: ir * self.ep_t,
            ep_t2: ir * self.ep_t2 - irr * half * self.ep_t * self.ep_t,
            ep_t3: irrr * half * self.ep_t * self.ep_t * self.ep_t - irr * self.ep_t2 * self.ep_t
                + ir * self.ep_t3,
        }
    }

    #[inline]
    fn cbrt(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn hypot(self, _other: Self) -> Self {
        unimplemented!()
        //(self * self + other * other).sqrt()
    }

    #[inline]
    fn sin(self) -> Self {
        let (s, c) = self.real.sin_cos();
        Dualt3 {
            real: s,
            ep_t: c * self.ep_t,
            ep_t2: c * self.ep_t2 - s * self.ep_t * self.ep_t / (T::one() + T::one()),
            ep_t3: c
                * (self.ep_t3
                    - self.ep_t * self.ep_t * self.ep_t / <T as NumCast>::from(6).unwrap())
                - s * (self.ep_t * self.ep_t2),
        }
    }

    #[inline]
    fn cos(self) -> Self {
        let (s, c) = self.real.sin_cos();
        Dualt3 {
            real: c,
            ep_t: -s * self.ep_t,
            ep_t2: -s * self.ep_t2 - c * self.ep_t * self.ep_t / (T::one() + T::one()),
            ep_t3: -s
                * (self.ep_t3
                    - self.ep_t * self.ep_t * self.ep_t / <T as NumCast>::from(6).unwrap())
                - c * (self.ep_t * self.ep_t2),
        }
    }

    #[inline]
    fn tan(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn asin(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn acos(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn atan(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn atan2(self, _other: Self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn sin_cos(self) -> (Self, Self) {
        unimplemented!()
    }

    #[inline]
    fn exp_m1(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn ln_1p(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn sinh(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn cosh(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn tanh(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn asinh(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn acosh(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn atanh(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn integer_decode(self) -> (u64, i16, i8) {
        unimplemented!()
    }

    #[inline]
    fn to_degrees(self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn to_radians(self) -> Self {
        unimplemented!()
    }
}

impl<T> Sum for Dualt3<T> {
    fn sum<I: Iterator<Item = Self>>(_iter: I) -> Self {
        todo!()
    }
}

impl<T: Float + Display + Debug + Default + Signed> Field for Dualt3<T> {}

impl<T: Float + Display + Debug + Default + Signed> RealNumberLike for Dualt3<T> {}

#[test]
fn simple_tests() {
    let a = Dualt3::new(1., 2., 3., 4.);
    let s = a.sin();
    let c = a.cos();

    println!("{:e}", s.powi(4).sqrt() + (c * c * c) / c - Dualt3::one());
    println!("{:e}", (s * 2.).exp() / s.exp().powi(2));

    todo!();
}
