extern crate cpython;
extern crate dual_num;
extern crate num;
extern crate serde;

use dual_num::{Allocator, DefaultAllocator, Dim, DimName, DualN, Owned};
use num::traits::ops::mul_add::MulAdd;
use num::traits::Inv;
use num::traits::{NumAssign, NumOps, NumRef};
use num::Complex;
use num::Float;
use num::Num;
use num::Signed;
use num::{NumCast, ToPrimitive};
use std::fmt;
use std::fmt::{Debug, Display, LowerExp};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};

mod deserialize;

pub trait Field
where
    Self: Num,
    Self: Mul<Self, Output = Self>,
    Self: MulAssign<Self>,
    Self: AddAssign<Self>,
    Self: SubAssign<Self>,
    Self: Div<Self, Output = Self>,
    Self: Add<Self, Output = Self>,
    Self: Sub<Self, Output = Self>,
    Self: Neg<Output = Self>,
    Self: Inv<Output = Self>,
    Self: Sum<Self>,
    Self: PartialEq,
    Self: Copy,
    //Self: From<f64>,
    Self: Default,
    Self: Debug,
    Self: Display,
{
}

pub trait RealNumberLike
where
    Self: Field,
    Self: Num,
    Self: NumCast,
    Self: Float,
    Self: NumAssign,
    Self: NumCast,
    Self: NumOps,
    Self: NumRef,
{
}

impl Field for f32 {}
impl Field for f64 {}
impl Field for f128::f128 {}

impl RealNumberLike for f64 {}
impl RealNumberLike for f32 {}
impl RealNumberLike for f128::f128 {}
impl<U, T: RealNumberLike + dual_num::FloatConst + Signed + 'static> RealNumberLike for DualN<T, U>
where
    U: Dim + DimName,
    DefaultAllocator: Allocator<T, U>,
    Owned<T, U>: Copy,
{
}

impl<T: RealNumberLike> Field for num::Complex<T> {}
impl<U, T: RealNumberLike + Signed + 'static> Field for DualN<T, U>
where
    U: Dim + DimName,
    DefaultAllocator: Allocator<T, U>,
    Owned<T, U>: Copy,
{
}

/// A generalization of a field with the reals as a basic component, with partial ordering
/// An example of a `RealField` is a dual.
pub trait RealField
where
    Self: Field,
    Self: PartialOrd,
    Self: Mul<f64, Output = Self>,
    Self: MulAssign<f64>,
    Self: Add<f64, Output = Self>,
    Self: AddAssign<f64>,
    Self: From<f64>,
    Self: Sub<f64, Output = Self>,
    Self: SubAssign<f64>,
    Self: Div<f64, Output = Self>,
    Self: PartialOrd<f64>,
    Self: Inv<Output = Self>,
{
}

impl<U> RealField for DualN<f64, U>
where
    U: Dim + DimName,
    DefaultAllocator: Allocator<f64, U>,
    Owned<f64, U>: Copy,
{
}

#[derive(Debug, Copy, Clone)]
pub struct LorentzVector<T: Field> {
    pub t: T,
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Field> Default for LorentzVector<T> {
    fn default() -> LorentzVector<T> {
        LorentzVector {
            t: T::default(),
            x: T::default(),
            y: T::default(),
            z: T::default(),
        }
    }
}

impl<T: Field> Display for LorentzVector<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "(t:{}, x:{}, y:{}, z:{})",
            self.t, self.x, self.y, self.z
        )
    }
}

impl<T: Field + LowerExp> LowerExp for LorentzVector<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "(t:{:e}, x:{:e}, y:{:e}, z:{:e})",
            self.t, self.x, self.y, self.z
        )
    }
}

impl<T: Field> LorentzVector<T> {
    #[inline]
    pub fn new() -> LorentzVector<T> {
        LorentzVector {
            t: T::default(),
            x: T::default(),
            y: T::default(),
            z: T::default(),
        }
    }

    #[inline]
    pub fn from_args(t: T, x: T, y: T, z: T) -> LorentzVector<T> {
        LorentzVector { t, x, y, z }
    }

    #[inline]
    pub fn from_slice(v: &[T]) -> LorentzVector<T> {
        let (t, x, y, z) = (v[0], v[1], v[2], v[3]);
        LorentzVector { t, x, y, z }
    }

    #[inline]
    pub fn from_vec(v: Vec<T>) -> LorentzVector<T> {
        let (t, x, y, z) = (v[0], v[1], v[2], v[3]);
        LorentzVector { t, x, y, z }
    }

    #[inline]
    pub fn dual(&self) -> LorentzVector<T> {
        LorentzVector {
            t: self.t,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    #[inline]
    pub fn square(&self) -> T {
        self.t * self.t - self.x * self.x - self.y * self.y - self.z * self.z
    }

    #[inline]
    pub fn dot(&self, other: &LorentzVector<T>) -> T {
        self.t * other.t - self.x * other.x - self.y * other.y - self.z * other.z
    }

    #[inline]
    pub fn spatial_squared(&self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[inline]
    pub fn spatial_dot(&self, other: &LorentzVector<T>) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    #[inline]
    pub fn euclidean_square(&self) -> T {
        self.t * self.t + self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[inline]
    pub fn comp_mul(&self, other: &LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector::from_args(
            self.t * other.t,
            self.x * other.x,
            self.y * other.y,
            self.z * other.z,
        )
    }

    #[inline]
    pub fn euclidean_dot(&self, other: &LorentzVector<T>) -> T {
        self.t * other.t + self.x * other.x + self.y * other.y + self.z * other.z
    }

    #[inline]
    pub fn map<F, U: Field>(&self, map: F) -> LorentzVector<U>
    where
        F: Fn(T) -> U,
    {
        LorentzVector {
            t: map(self.t),
            x: map(self.x),
            y: map(self.y),
            z: map(self.z),
        }
    }

    #[inline]
    pub fn from<U: Field + Into<T>>(a: LorentzVector<U>) -> Self {
        LorentzVector {
            t: a.t.into(),
            x: a.x.into(),
            y: a.y.into(),
            z: a.z.into(),
        }
    }

    #[inline]
    pub fn convert<U: Field + From<T>>(&self) -> LorentzVector<U> {
        LorentzVector {
            t: self.t.into(),
            x: self.x.into(),
            y: self.y.into(),
            z: self.z.into(),
        }
    }
}

impl<T: Field + ToPrimitive> LorentzVector<T> {
    #[inline]
    pub fn cast<U: Field + NumCast>(&self) -> LorentzVector<U> {
        LorentzVector {
            t: <U as NumCast>::from(self.t).unwrap(),
            x: <U as NumCast>::from(self.x).unwrap(),
            y: <U as NumCast>::from(self.y).unwrap(),
            z: <U as NumCast>::from(self.z).unwrap(),
        }
    }
}

impl<'a, T: Field + MulAdd<Output = T>> LorentzVector<T> {
    #[inline]
    pub fn square_impr(&self) -> T {
        self.t.mul_add(self.t, -self.spatial_squared_impr())
    }

    #[inline]
    pub fn spatial_squared_impr(&self) -> T {
        self.x
            .mul_add(self.x, self.y.mul_add(self.y, self.z * self.z))
    }

    #[inline]
    pub fn dot_impr(&self, other: &LorentzVector<T>) -> T {
        self.t.mul_add(other.t, -self.spatial_dot_impr(other))
    }

    #[inline]
    pub fn spatial_dot_impr(&self, other: &LorentzVector<T>) -> T {
        self.x
            .mul_add(other.x, self.y.mul_add(other.y, self.z * other.z))
    }
}

impl<'a, T: Field> Neg for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn neg(self) -> LorentzVector<T> {
        LorentzVector {
            t: -self.t,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<'a, T: Field> Neg for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn neg(self) -> LorentzVector<T> {
        LorentzVector {
            t: -self.t,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<'a, T: Field> Add<&'a LorentzVector<T>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: &'a LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: Field> Add<LorentzVector<T>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: LorentzVector<T>) -> LorentzVector<T> {
        self.add(&other)
    }
}

impl<'a, T: Field> Add<&'a LorentzVector<T>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: &'a LorentzVector<T>) -> LorentzVector<T> {
        &self + other
    }
}

impl<'a, T: Field> Add<LorentzVector<T>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: LorentzVector<T>) -> LorentzVector<T> {
        &self + &other
    }
}

impl<'a, T: Field> AddAssign<LorentzVector<T>> for LorentzVector<T> {
    #[inline]
    fn add_assign(&mut self, other: LorentzVector<T>) {
        self.t += other.t;
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl<'a, T: Field> SubAssign<LorentzVector<T>> for LorentzVector<T> {
    #[inline]
    fn sub_assign(&mut self, other: LorentzVector<T>) {
        self.t -= other.t;
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl<'a, T: Field> Sub<&'a LorentzVector<T>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: Field> Sub<LorentzVector<T>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: LorentzVector<T>) -> LorentzVector<T> {
        self.sub(&other)
    }
}

impl<'a, T: Field> Sub<LorentzVector<T>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: LorentzVector<T>) -> LorentzVector<T> {
        &self - &other
    }
}

impl<'a, T: Field> Sub<&'a LorentzVector<T>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<T>) -> LorentzVector<T> {
        &self - other
    }
}

impl<'a, T: Field> Mul<T> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn mul(self, other: T) -> LorentzVector<T> {
        LorentzVector {
            t: self.t * other,
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<'a, T: RealField> Mul<f64> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn mul(self, other: f64) -> LorentzVector<T> {
        LorentzVector {
            t: self.t * other,
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<'a, T: RealField> Mul<f64> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn mul(self, other: f64) -> LorentzVector<T> {
        LorentzVector {
            t: self.t * other,
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<'a, T: Field> Mul<T> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn mul(self, other: T) -> LorentzVector<T> {
        LorentzVector {
            t: self.t * other,
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<'a, T: Field> Div<T> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn div(self, other: T) -> LorentzVector<T> {
        let o = other.inv();
        self * o
    }
}

impl<'a, T: Field> Div<T> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn div(self, other: T) -> LorentzVector<T> {
        let o = other.inv();
        self * o
    }
}

impl<'a, T: Field> Inv for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn inv(self) -> LorentzVector<T> {
        LorentzVector {
            t: self.t.inv(),
            x: self.x.inv(),
            y: self.y.inv(),
            z: self.z.inv(),
        }
    }
}

impl<'a, T: RealNumberLike> Sub<&'a LorentzVector<T>> for &'a LorentzVector<Complex<T>> {
    type Output = LorentzVector<Complex<T>>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<T>) -> LorentzVector<Complex<T>> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: RealNumberLike> Sub<&'a LorentzVector<T>> for LorentzVector<Complex<T>> {
    type Output = LorentzVector<Complex<T>>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<T>) -> LorentzVector<Complex<T>> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: RealNumberLike> Sub<LorentzVector<T>> for LorentzVector<Complex<T>> {
    type Output = LorentzVector<Complex<T>>;

    #[inline]
    fn sub(self, other: LorentzVector<T>) -> LorentzVector<Complex<T>> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: RealNumberLike> Add<&'a LorentzVector<T>> for &'a LorentzVector<Complex<T>> {
    type Output = LorentzVector<Complex<T>>;

    #[inline]
    fn add(self, other: &'a LorentzVector<T>) -> LorentzVector<Complex<T>> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: RealNumberLike> Add<&'a LorentzVector<T>> for LorentzVector<Complex<T>> {
    type Output = LorentzVector<Complex<T>>;

    #[inline]
    fn add(self, other: &'a LorentzVector<T>) -> LorentzVector<Complex<T>> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: RealNumberLike> Add<LorentzVector<T>> for LorentzVector<Complex<T>> {
    type Output = LorentzVector<Complex<T>>;

    #[inline]
    fn add(self, other: LorentzVector<T>) -> LorentzVector<Complex<T>> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: Field> MulAssign<T> for LorentzVector<T> {
    #[inline]
    fn mul_assign(&mut self, other: T) {
        self.t *= other;
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

impl<'a, T: RealField> MulAssign<f64> for LorentzVector<T> {
    #[inline]
    fn mul_assign(&mut self, other: f64) {
        self.t *= other;
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }
}

impl<'a, T: RealField> Sub<&'a LorentzVector<f64>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: RealField> Sub<&'a LorentzVector<f64>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: RealField> Sub<LorentzVector<f64>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<T: RealField> Sub<LorentzVector<f64>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn sub(self, other: LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: RealField> Add<LorentzVector<f64>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: RealField> Add<&'a LorentzVector<f64>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: &'a LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: RealField> Add<&'a LorentzVector<f64>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: &'a LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<T: RealField> Add<LorentzVector<f64>> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn add(self, other: LorentzVector<f64>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: RealField> AddAssign<LorentzVector<f64>> for LorentzVector<T> {
    #[inline]
    fn add_assign(&mut self, other: LorentzVector<f64>) {
        self.t += other.t;
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl<'a, T: RealField> SubAssign<LorentzVector<f64>> for LorentzVector<T> {
    #[inline]
    fn sub_assign(&mut self, other: LorentzVector<f64>) {
        self.t -= other.t;
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl<T: Field> Index<usize> for LorentzVector<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &T {
        match index {
            0 => &self.t,
            1 => &self.x,
            2 => &self.y,
            3 => &self.z,
            _ => panic!("Index is not between 0 and 3"),
        }
    }
}

impl<T: Field> IndexMut<usize> for LorentzVector<T> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut T {
        match index {
            0 => &mut self.t,
            1 => &mut self.x,
            2 => &mut self.y,
            3 => &mut self.z,
            _ => panic!("Index is not between 0 and 3"),
        }
    }
}

impl<T: Float + Field> LorentzVector<T> {
    #[inline]
    pub fn spatial_distance(&self) -> T {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    #[inline]
    pub fn euclidean_distance(&self) -> T {
        (self.t * self.t + self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn boost(&self, boost_vector: &LorentzVector<T>) -> LorentzVector<T> {
        let b2 = boost_vector.spatial_squared();
        let gamma = (T::one() - b2).sqrt().inv();

        let bp = self.spatial_dot(boost_vector);
        let gamma2 = if b2 > T::zero() {
            (gamma - T::one()) / b2
        } else {
            T::zero()
        };
        let factor = gamma2 * bp + gamma * self.t;
        LorentzVector::from_args(
            gamma * (self.t + bp),
            boost_vector.x.mul_add(factor, self.x),
            boost_vector.y.mul_add(factor, self.y),
            boost_vector.z.mul_add(factor, self.z),
        )
    }

    /// Compute transverse momentum.
    #[inline]
    pub fn pt(&self) -> T {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Compute pseudorapidity.
    #[inline]
    pub fn pseudo_rap(&self) -> T {
        let pt = self.pt();
        if pt < T::epsilon() && self.z.abs() < T::epsilon() {
            if self.z > T::zero() {
                return T::max_value();
            } else {
                return T::min_value();
            }
        }
        let th = pt.atan2(self.z);
        -(th / (T::one() + T::one())).tan().ln()
    }

    /// Compute the phi-angle separation with p2.
    pub fn getdelphi(&self, p2: &LorentzVector<T>) -> T {
        let pt1 = self.pt();
        let pt2 = p2.pt();
        if pt1 == T::zero() || pt2 == T::zero() {
            return T::max_value();
        }

        let mut tmp = self.x * p2.x + self.y * p2.y;
        tmp = tmp / (pt1 * pt2);
        if tmp.abs() > T::one() + T::epsilon() {
            panic!("Cosine larger than 1. in phase-space cuts.")
        }
        if tmp.abs() > T::one() {
            (tmp / tmp.abs()).acos()
        } else {
            tmp.acos()
        }
    }

    /// Compute the deltaR separation with momentum p2.
    #[inline]
    pub fn delta_r(&self, p2: &LorentzVector<T>) -> T {
        let delta_eta = self.pseudo_rap() - p2.pseudo_rap();
        let delta_phi = self.getdelphi(p2);
        (delta_eta * delta_eta + delta_phi * delta_phi).sqrt()
    }

    /// Apply a pure boost that sends p into q to this `LorentzVector`.
    /// For details, see appendix A.2.2 of Simone Lionetti's PhD thesis.
    ///
    /// * `p` - Starting Lorentz vector to define the boost.
    /// * `q` - Target Lorentz vector to define the boost.
    pub fn boost_from_to(&self, p: &LorentzVector<T>, q: &LorentzVector<T>) -> LorentzVector<T> {
        let eps = T::epsilon() + T::epsilon();
        let p_abs = p.euclidean_distance();
        let q_abs = q.euclidean_distance();

        if (p - q).spatial_distance() < eps * eps {
            return LorentzVector::new();
        }

        let mut n_vec = q - p;
        n_vec = n_vec / n_vec.spatial_distance();

        let na = LorentzVector::from_args(T::one(), n_vec.x, n_vec.y, n_vec.z);
        let nb = LorentzVector::from_args(T::one(), -n_vec.x, -n_vec.y, -n_vec.z);

        let p_plus = p.dot(&nb);
        let p_minus = p.dot(&na);
        let q_plus = q.dot(&nb);
        let q_minus = q.dot(&na);
        let ratioa;
        let ratiob;
        if p_minus / p_abs < eps && q_minus / q_abs < eps {
            if p_plus / p_abs < eps && q_plus / q_abs < eps {
                ratioa = T::one();
                ratiob = T::one();
            } else {
                ratiob = q_plus / p_plus;
                ratioa = T::one() / ratiob;
            }
        } else {
            if p_plus / p_abs < eps && q_plus / q_abs < eps {
                ratioa = q_minus / p_minus;
                ratiob = T::one() / ratioa;
            } else {
                ratioa = q_minus / p_minus;
                ratiob = q_plus / p_plus;
            }
        }

        let plus = self.dot(&nb);
        let minus = self.dot(&na);

        self + na * (ratiob - T::one()) / (T::one() + T::one()) * plus
            + nb * (ratioa - T::one()) / (T::one() + T::one()) * minus
    }
}

impl LorentzVector<f64> {
    /// Boost this kinematic configuration from the center of mass frame to the lab frame
    /// given specified Bjorken x's x1 and x2.
    /// This function needs to be cleaned up and built in a smarter way as the boost vector can be written
    /// down explicitly as a function of x1, x2 and the beam energies.
    pub fn boost_from_com_to_lab_frame(
        momenta: &mut [LorentzVector<f64>],
        x1: f64,
        x2: f64,
        ebeam1: f64,
        ebeam2: f64,
    ) {
        let target_summed =
            LorentzVector::from_args(x1 * ebeam1, 0., 0., (x1 * ebeam1).copysign(momenta[0].z))
                + LorentzVector::from_args(
                    x2 * ebeam2,
                    0.,
                    0.,
                    (x2 * ebeam2).copysign(momenta[1].z),
                );

        let source_summed =
            LorentzVector::from_args(2. * (x1 * x2 * ebeam1 * ebeam2).sqrt(), 0., 0., 0.);

        // We want to send the source to the target
        for vec in momenta {
            *vec = vec.boost_from_to(&source_summed, &target_summed);
        }
    }
}

impl<T: RealNumberLike> LorentzVector<T> {
    #[inline]
    pub fn to_complex(&self, real: bool) -> LorentzVector<Complex<T>> {
        if real {
            LorentzVector {
                t: Complex::new(self.t, T::zero()),
                x: Complex::new(self.x, T::zero()),
                y: Complex::new(self.y, T::zero()),
                z: Complex::new(self.z, T::zero()),
            }
        } else {
            LorentzVector {
                t: Complex::new(T::zero(), self.t),
                x: Complex::new(T::zero(), self.x),
                y: Complex::new(T::zero(), self.y),
                z: Complex::new(T::zero(), self.z),
            }
        }
    }
}

impl<T: RealNumberLike> LorentzVector<Complex<T>> {
    #[inline]
    pub fn real(&self) -> LorentzVector<T> {
        LorentzVector {
            t: self.t.re,
            x: self.x.re,
            y: self.y.re,
            z: self.z.re,
        }
    }

    #[inline]
    pub fn imag(&self) -> LorentzVector<T> {
        LorentzVector {
            t: self.t.im,
            x: self.x.im,
            y: self.y.im,
            z: self.z.im,
        }
    }
}

impl<U: dual_num::Dim + dual_num::DimName, T: RealNumberLike + Signed + 'static>
    LorentzVector<DualN<T, U>>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    #[inline]
    pub fn real(&self) -> LorentzVector<T> {
        self.map(|x| x.real())
    }
}
