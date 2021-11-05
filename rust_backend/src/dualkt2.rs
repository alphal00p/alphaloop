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
pub struct Dualkt2<T> {
    pub real: T,
    pub ep_k0: T,
    pub ep_t: T,
    pub ep_k0_t: T,
    pub ep_t2: T,
}

impl<T> Dualkt2<T> {
    pub fn new(real: T, ep_k0: T, ep_t: T, ep_k0_t: T, ep_t2: T) -> Dualkt2<T> {
        Dualkt2 {
            real,
            ep_k0,
            ep_t,
            ep_k0_t,
            ep_t2,
        }
    }

    #[inline]
    pub fn from_real(real: T) -> Dualkt2<T>
    where
        T: Zero,
    {
        Dualkt2 {
            real,
            ep_k0: T::zero(),
            ep_t: T::zero(),
            ep_k0_t: T::zero(),
            ep_t2: T::zero(),
        }
    }

    #[inline]
    pub fn map<F, U>(&self, f: F) -> Dualkt2<U>
    where
        F: Fn(&T) -> U,
    {
        Dualkt2 {
            real: f(&self.real),
            ep_k0: f(&self.ep_k0),
            ep_t: f(&self.ep_t),
            ep_k0_t: f(&self.ep_k0_t),
            ep_t2: f(&self.ep_t2),
        }
    }

    #[inline]
    pub fn zip_map<F, U, K>(&self, rhs: &Dualkt2<U>, f: F) -> Dualkt2<K>
    where
        F: Fn(&T, &U) -> K,
    {
        Dualkt2 {
            real: f(&self.real, &rhs.real),
            ep_k0: f(&self.ep_k0, &rhs.ep_k0),
            ep_t: f(&self.ep_t, &rhs.ep_t),
            ep_k0_t: f(&self.ep_k0_t, &rhs.ep_k0_t),
            ep_t2: f(&self.ep_t2, &rhs.ep_t2),
        }
    }
}

impl<T: Copy> Dualkt2<T> {
    fn substitute_real(&self, new_real: T) -> Dualkt2<T> {
        Dualkt2 {
            real: new_real,
            ep_k0: self.ep_k0,
            ep_t: self.ep_t,
            ep_k0_t: self.ep_k0_t,
            ep_t2: self.ep_t2,
        }
    }
}

impl<T: Zero> From<T> for Dualkt2<T> {
    #[inline]
    fn from(real: T) -> Dualkt2<T> {
        Dualkt2::from_real(real)
    }
}

impl<T: Display> Display for Dualkt2<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let precision = f.precision().unwrap_or(4);
        write!(
            f,
            "({:.p$} + {:.p$}\u{03B5}_k + {:.p$}\u{03B5}_t  + {:.p$}\u{03B5}_k\u{03B5}_t + {:.p$}\u{03B5}_t^2)",
            self.real,
            self.ep_k0,
            self.ep_t,
            self.ep_k0_t,
            self.ep_t2,
            p = precision
        )
    }
}

impl<T: Display + LowerExp> LowerExp for Dualkt2<T> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        let precision = f.precision().unwrap_or(4);
        write!(
            f,
            "({:.p$e} + {:.p$e}\u{03B5}_k + {:.p$e}\u{03B5}_t + {:.p$e}\u{03B5}_k\u{03B5}_t + {:.p$e}\u{03B5}_t^2)",
            self.real,
            self.ep_k0,
            self.ep_t,
            self.ep_k0_t,
            self.ep_t2,
            p = precision
        )
    }
}

impl<T: PartialEq> PartialEq<Self> for Dualkt2<T> {
    #[inline]
    fn eq(&self, rhs: &Self) -> bool {
        self.real == rhs.real
    }
}

impl<T: PartialOrd> PartialOrd<Self> for Dualkt2<T> {
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

impl<T: PartialEq> PartialEq<T> for Dualkt2<T> {
    #[inline]
    fn eq(&self, rhs: &T) -> bool {
        self.real == *rhs
    }
}

impl<T: PartialOrd> PartialOrd<T> for Dualkt2<T> {
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
        impl<T: ToPrimitive> ToPrimitive for Dualkt2<T>
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
        impl<T: FromPrimitive + Zero> FromPrimitive for Dualkt2<T>
           {
            $(
                #[inline]
                fn $name(n: $ty) -> Option<Dualkt2<T>> {
                    T::$name(n).map(Dualkt2::from_real)
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

impl<T: Num + Copy> Add<T> for Dualkt2<T> {
    type Output = Dualkt2<T>;

    #[inline]
    fn add(self, rhs: T) -> Dualkt2<T> {
        self.substitute_real(self.real + rhs)
    }
}

impl<T: Num + Copy> AddAssign<T> for Dualkt2<T> {
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        *self = (*self) + Dualkt2::from_real(rhs)
    }
}

impl<T: Num + Copy> Sub<T> for Dualkt2<T> {
    type Output = Dualkt2<T>;

    #[inline]
    fn sub(self, rhs: T) -> Dualkt2<T> {
        self.substitute_real(self.real - rhs)
    }
}

impl<T: Num + Copy> SubAssign<T> for Dualkt2<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: T) {
        *self = (*self) - Dualkt2::from_real(rhs)
    }
}

impl<T: Num + Copy> Mul<T> for Dualkt2<T> {
    type Output = Dualkt2<T>;

    #[inline]
    fn mul(self, rhs: T) -> Dualkt2<T> {
        self.map(|x| *x * rhs)
    }
}

impl<T: Num + Copy> MulAssign<T> for Dualkt2<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = (*self) * Dualkt2::from_real(rhs)
    }
}

impl<T: Num + Copy> Div<T> for Dualkt2<T> {
    type Output = Dualkt2<T>;

    #[inline]
    fn div(self, rhs: T) -> Dualkt2<T> {
        self / Dualkt2::from_real(rhs)
    }
}

impl<T: Num + Copy> DivAssign<T> for Dualkt2<T> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = (*self) / Dualkt2::from_real(rhs)
    }
}

impl<T: Signed + Copy> Neg for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        self.map(|x| x.neg())
    }
}

impl<T: Num + Copy> Add<Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        self.zip_map(&rhs, |x, y| *x + *y)
    }
}

impl<T: Num + Copy> Add<&Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: &Self) -> Self {
        self.zip_map(&rhs, |x, y| *x + *y)
    }
}

impl<T: Num + Copy> AddAssign<Self> for Dualkt2<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = (*self) + rhs
    }
}

impl<T: Num + Copy> Sub<Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        self.zip_map(&rhs, |x, y| *x - *y)
    }
}

impl<T: Num + Copy> Sub<&Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: &Self) -> Self {
        self.zip_map(&rhs, |x, y| *x - *y)
    }
}

impl<T: Num + Copy> SubAssign<Self> for Dualkt2<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = (*self) - rhs
    }
}

impl<T: Num + Copy> Mul<&Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: &Self) -> Self {
        Dualkt2 {
            real: self.real * rhs.real,
            ep_k0: self.real * rhs.ep_k0 + self.ep_k0 * rhs.real,
            ep_t: self.real * rhs.ep_t + self.ep_t * rhs.real,
            ep_k0_t: self.real * rhs.ep_k0_t
                + self.ep_k0_t * rhs.real
                + self.ep_k0 * rhs.ep_t
                + self.ep_t * rhs.ep_k0,
            ep_t2: self.ep_t * rhs.ep_t + self.real * rhs.ep_t2 + self.ep_t2 * rhs.real,
        }
    }
}

impl<T: Num + Copy> Mul<Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        self * (&rhs)
    }
}

impl<T: Num + Copy> MulAssign<Self> for Dualkt2<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = (*self) * rhs
    }
}

macro_rules! impl_mul_add {
    ($(<$a:ident, $b:ident>),*) => {
        $(
            impl<T: Num + Copy + Mul + Add> MulAdd<$a, $b> for Dualkt2<T> {
                type Output = Dualkt2<T>;

                #[inline]
                fn mul_add(self, a: $a, b: $b) -> Dualkt2<T> {
                    (self * a) + b
                }
            }

            impl<T: Num + Copy + Mul + Add> MulAddAssign<$a, $b> for Dualkt2<T> {
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

impl<T: Num + Copy> Div<&Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: &Self) -> Self {
        let r1 = T::one() / rhs.real;
        let r2 = r1 * r1;
        let r3 = r2 * r1;

        Dualkt2 {
            real: self.real * r1,
            ep_k0: self.ep_k0 * r1 - self.real * rhs.ep_k0 * r2,
            ep_t: self.ep_t * r1 - self.real * rhs.ep_t * r2,
            ep_k0_t: self.ep_k0_t * r1
                - rhs.ep_k0_t * self.real * r2
                - self.ep_k0 * rhs.ep_t * r2
                - self.ep_t * rhs.ep_k0 * r2
                + (T::one() + T::one()) * self.real * rhs.ep_k0 * rhs.ep_t * r3,
            ep_t2: self.ep_t2 * r1 - self.ep_t * rhs.ep_t * r2 - self.real * rhs.ep_t2 * r2
                + self.real * rhs.ep_t * rhs.ep_t * r3,
        }
    }
}

impl<T: Num + Copy> Div<Self> for Dualkt2<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self {
        self / (&rhs)
    }
}

impl<T: Num + Copy> DivAssign<Self> for Dualkt2<T> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = (*self) / rhs
    }
}

impl<T: Num + Copy> Rem<Self> for Dualkt2<T> {
    type Output = Self;

    fn rem(self, _: Self) -> Self {
        unimplemented!()
    }
}

impl<T: Num + Copy> Rem<&Self> for Dualkt2<T> {
    type Output = Self;

    fn rem(self, _: &Self) -> Self {
        unimplemented!()
    }
}

impl<T: Num + Copy> RemAssign<Self> for Dualkt2<T> {
    fn rem_assign(&mut self, _: Self) {
        unimplemented!()
    }
}

impl<T, P: Into<Dualkt2<T>>> Pow<P> for Dualkt2<T>
where
    Dualkt2<T>: Float,
{
    type Output = Dualkt2<T>;

    #[inline]
    fn pow(self, rhs: P) -> Dualkt2<T> {
        self.powf(rhs.into())
    }
}

impl<T> Inv for Dualkt2<T>
where
    Self: One + Div<Output = Self>,
{
    type Output = Dualkt2<T>;

    #[inline]
    fn inv(self) -> Dualkt2<T> {
        Dualkt2::one() / self
    }
}

impl<T> Signed for Dualkt2<T>
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
        Dualkt2::from_real(self.real.signum())
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

impl<T: Unsigned + Copy> Unsigned for Dualkt2<T> {}

impl<T: Num + Zero + Copy> Zero for Dualkt2<T> {
    #[inline]
    fn zero() -> Dualkt2<T> {
        Dualkt2::from_real(T::zero())
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.real.is_zero()
    }
}

impl<T: Num + One + Copy> One for Dualkt2<T> {
    #[inline]
    fn one() -> Dualkt2<T> {
        Dualkt2::from_real(T::one())
    }

    #[inline]
    fn is_one(&self) -> bool
    where
        Self: PartialEq,
    {
        self.real.is_one()
    }
}

impl<T: Num + Copy> Num for Dualkt2<T> {
    type FromStrRadixErr = <T as Num>::FromStrRadixErr;

    #[inline]
    fn from_str_radix(str: &str, radix: u32) -> Result<Dualkt2<T>, Self::FromStrRadixErr> {
        <T as Num>::from_str_radix(str, radix).map(Dualkt2::from_real)
    }
}

impl<T: Float> NumCast for Dualkt2<T> {
    #[inline]
    fn from<P: ToPrimitive>(n: P) -> Option<Dualkt2<T>> {
        <T as NumCast>::from(n).map(Dualkt2::from_real)
    }
}

macro_rules! impl_float_const {
    ($($c:ident),*) => {
        $(
            fn $c() -> Dualkt2<T> { Dualkt2::from_real(T::$c()) }
        )*
    }
}

impl<T: FloatConst + Zero> FloatConst for Dualkt2<T> {
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
            fn $prop() -> Self { Dualkt2::from_real(<T as Float>::$prop()) }
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
            fn $op(self) -> Self { Dualkt2::from_real(self.real.$op()) }
        )*
    }
}

impl<T: Signed + Float + Debug> Float for Dualkt2<T> {
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
        Dualkt2::from_real(self.real.signum())
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
        let r = self.real.powi(n);
        let dr = <T as NumCast>::from(n).unwrap() * r / self.real;
        let ddr = <T as NumCast>::from(n).unwrap() * r / self.real / self.real;
        Dualkt2 {
            real: r,
            ep_k0: self.ep_k0 * dr,
            ep_t: self.ep_t * dr,
            ep_k0_t: self.ep_k0_t * dr
                + self.ep_k0 * self.ep_t * (<T as NumCast>::from(n).unwrap() * ddr - ddr),
            ep_t2: self.ep_t2 * dr
                + self.ep_t * self.ep_t / (T::one() + T::one())
                    * (<T as NumCast>::from(n).unwrap() * ddr - ddr),
        }
    }

    #[inline]
    fn powf(self, _n: Self) -> Self {
        unimplemented!()
    }

    #[inline]
    fn exp(self) -> Self {
        let r = self.real.exp();
        Dualkt2 {
            real: r,
            ep_k0: r * self.ep_k0,
            ep_t: r * self.ep_t,
            ep_k0_t: r * self.ep_k0_t + r * self.ep_k0 * self.ep_t,
            ep_t2: r * self.ep_t2 + r * self.ep_t * self.ep_t / (T::one() + T::one()),
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
        Dualkt2 {
            real: r,
            ep_k0: self.ep_k0 * ir,
            ep_t: self.ep_t * ir,
            ep_k0_t: self.ep_k0_t * ir - self.ep_k0 * self.ep_t * irr,
            ep_t2: self.ep_t2 * ir - self.ep_t * self.ep_t * irr * half,
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
        let c = self.real.cos();
        let s = self.real.sin();
        Dualkt2 {
            real: s,
            ep_k0: self.ep_k0 * c,
            ep_t: self.ep_t * c,
            ep_k0_t: self.ep_k0_t * c - self.ep_k0 * self.ep_t * s,
            ep_t2: self.ep_t2 * c - self.ep_t * self.ep_t * s / (T::one() + T::one()),
        }
    }

    #[inline]
    fn cos(self) -> Self {
        let c = self.real.cos();
        let s = self.real.sin();
        Dualkt2 {
            real: c,
            ep_k0: -self.ep_k0 * s,
            ep_t: -self.ep_t * s,
            ep_k0_t: -self.ep_k0_t * s - self.ep_k0 * self.ep_t * c,
            ep_t2: -self.ep_t2 * s - self.ep_t * self.ep_t * c / (T::one() + T::one()),
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

impl<T> Sum for Dualkt2<T> {
    fn sum<I: Iterator<Item = Self>>(_iter: I) -> Self {
        todo!()
    }
}

impl<T: Float + Display + Debug + Default + Signed> Field for Dualkt2<T> {}

impl<T: Float + Display + Debug + Default + Signed> RealNumberLike for Dualkt2<T> {}
