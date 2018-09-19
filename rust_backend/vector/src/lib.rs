extern crate num;

use std::ops::{Add, Mul, Sub};
use std::fmt::{Debug, Display};

type Complex = num::Complex<f64>;

pub trait Field where
    Self: Mul<Self, Output=Self>,
    Self: Add<Self, Output=Self>,
    Self: Sub<Self, Output=Self>,
    Self: Copy,
    Self: Default,
    Self: Debug,
    Self: Display
{
}

impl Field for f64 {}
impl Field for Complex {}

#[derive(Debug, Clone)]
pub struct LorentzVector<T: Field> {
    x0: T,
    x1: T,
    x2: T,
    x3: T,
}

impl<T: Field> LorentzVector<T> {
    pub fn new() -> LorentzVector<T> {
        LorentzVector {
            x0: T::default(),
            x1: T::default(),
            x2: T::default(),
            x3: T::default(),
        }
    }

    pub fn from(v: Vec<T>) -> LorentzVector<T> {
        let (t, x, y, z) = (v[0], v[1], v[2], v[3]);
        LorentzVector {
            x0: t,
            x1: x,
            x2: y,
            x3: z
        }
    }

    pub fn square(&self) -> T {
        self.x1 * self.x1 + self.x2 * self.x2 + self.x3 * self.x3 - self.x0 * self.x0
    }

    pub fn dot(&self, other: &LorentzVector<T>) -> T {
        self.x0 * other.x0 + self.x1 * other.x1 + self.x2 * other.x2 + self.x3 * other.x3
    }
}

impl<T: Field> Add<&LorentzVector<T>> for &LorentzVector<T> {
    type Output = LorentzVector<T>;

    fn add(self, other: &LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector {
            x0: self.x0 + other.x0,
            x1: self.x1 + other.x1,
            x2: self.x2 + other.x2,
            x3: self.x3 + other.x3,
        }
    }
}

impl<T: Field> Sub<&LorentzVector<T>> for &LorentzVector<T> {
    type Output = LorentzVector<T>;

    fn sub(self, other: &LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector {
            x0: self.x0 - other.x0,
            x1: self.x1 - other.x1,
            x2: self.x2 - other.x2,
            x3: self.x3 - other.x3,
        }
    }
}

impl<T: Field> Mul<T> for &LorentzVector<T> {
    type Output = LorentzVector<T>;

    fn mul(self, other: T) -> LorentzVector<T> {
        LorentzVector {
            x0: self.x0 * other,
            x1: self.x1 * other,
            x2: self.x2 * other,
            x3: self.x3 * other,
        }
    }
}

impl Sub<&LorentzVector<f64>> for &LorentzVector<Complex> {
    type Output = LorentzVector<Complex>;

    fn sub(self, other: &LorentzVector<f64>) -> LorentzVector<Complex> {
        LorentzVector {
            x0: self.x0 - other.x0,
            x1: self.x1 - other.x1,
            x2: self.x2 - other.x2,
            x3: self.x3 - other.x2,
        }
    }
}
