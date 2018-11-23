extern crate dual_num;
extern crate num;

use dual_num::Dual;
use num::traits::Inv;
use num::Float;
use std::fmt;
use std::fmt::{Debug, Display};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, Sub};

type Complex = num::Complex<f64>;

pub trait Field
where
    Self: Mul<Self, Output = Self>,
    Self: MulAssign<Self>,
    Self: AddAssign<Self>,
    Self: Div<Self, Output = Self>,
    Self: Add<Self, Output = Self>,
    Self: Sub<Self, Output = Self>,
    Self: Neg<Output = Self>,
    Self: Sum<Self>,
    Self: PartialEq,
    Self: Copy,
    Self: From<f64>,
    Self: Default,
    Self: Debug,
    Self: Display,
{
}

impl Field for f64 {}
impl Field for Complex {}
impl Field for Dual<f64> {}

/// A generalization of a field with the reals as a basic component, with partial ordering
/// An example of a `RealField` is a dual.
pub trait RealField
where
    Self: Field,
    Self: PartialOrd,
    Self: Mul<f64, Output = Self>,
    Self: Add<f64, Output = Self>,
    Self: Sub<f64, Output = Self>,
    Self: Div<f64, Output = Self>,
    Self: PartialOrd<f64>,
    Self: Inv<Output = Self>,
{
}

impl RealField for f64 {}
impl RealField for Dual<f64> {}

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
    pub fn euclidean_square(&self) -> T {
        self.t * self.t + self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[inline]
    pub fn euclidean_dot(&self, other: &LorentzVector<T>) -> T {
        self.t * other.t + self.x * other.x + self.y * other.y + self.z * other.z
    }

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
        LorentzVector {
            t: self.t / other,
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl<'a, T: Field> Div<T> for LorentzVector<T> {
    type Output = LorentzVector<T>;

    #[inline]
    fn div(self, other: T) -> LorentzVector<T> {
        LorentzVector {
            t: self.t / other,
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl<T: RealField> LorentzVector<T> {
    pub fn from_f64(a: LorentzVector<f64>) -> Self {
        LorentzVector {
            t: From::from(a.t),
            x: From::from(a.x),
            y: From::from(a.y),
            z: From::from(a.z),
        }
    }
}

impl<'a> Sub<&'a LorentzVector<f64>> for &'a LorentzVector<Complex> {
    type Output = LorentzVector<Complex>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<f64>) -> LorentzVector<Complex> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a> Sub<&'a LorentzVector<f64>> for LorentzVector<Complex> {
    type Output = LorentzVector<Complex>;

    #[inline]
    fn sub(self, other: &'a LorentzVector<f64>) -> LorentzVector<Complex> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a> Add<&'a LorentzVector<f64>> for &'a LorentzVector<Complex> {
    type Output = LorentzVector<Complex>;

    #[inline]
    fn add(self, other: &'a LorentzVector<f64>) -> LorentzVector<Complex> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a> Add<&'a LorentzVector<f64>> for LorentzVector<Complex> {
    type Output = LorentzVector<Complex>;

    #[inline]
    fn add(self, other: &'a LorentzVector<f64>) -> LorentzVector<Complex> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
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

impl LorentzVector<f64> {
    #[inline]
    pub fn to_complex(&self, real: bool) -> LorentzVector<Complex> {
        if real {
            LorentzVector {
                t: Complex::new(self.t, 0.0),
                x: Complex::new(self.x, 0.0),
                y: Complex::new(self.y, 0.0),
                z: Complex::new(self.z, 0.0),
            }
        } else {
            LorentzVector {
                t: Complex::new(0.0, self.t),
                x: Complex::new(0.0, self.x),
                y: Complex::new(0.0, self.y),
                z: Complex::new(0.0, self.z),
            }
        }
    }
}

impl LorentzVector<Complex> {
    #[inline]
    pub fn real(&self) -> LorentzVector<f64> {
        LorentzVector {
            t: self.t.re,
            x: self.x.re,
            y: self.y.re,
            z: self.z.re,
        }
    }

    #[inline]
    pub fn imag(&self) -> LorentzVector<f64> {
        LorentzVector {
            t: self.t.im,
            x: self.x.im,
            y: self.y.im,
            z: self.z.im,
        }
    }
}
