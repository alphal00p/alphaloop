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
        LorentzVector {
            t: t,
            x: x,
            y: y,
            z: z,
        }
    }

    #[inline]
    pub fn from_vec(v: Vec<T>) -> LorentzVector<T> {
        let (t, x, y, z) = (v[0], v[1], v[2], v[3]);
        LorentzVector {
            t: t,
            x: x,
            y: y,
            z: z,
        }
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
