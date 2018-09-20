extern crate num;

use std::fmt::{Debug, Display};
use std::ops::{Add, Mul, Sub};

type Complex = num::Complex<f64>;

pub trait Field
where
    Self: Mul<Self, Output = Self>,
    Self: Add<Self, Output = Self>,
    Self: Sub<Self, Output = Self>,
    Self: Copy,
    Self: Default,
    Self: Debug,
    Self: Display,
{
}

impl Field for f64 {}
impl Field for Complex {}

#[derive(Debug, Clone)]
pub struct LorentzVector<T: Field> {
    pub t: T,
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Field> LorentzVector<T> {
    pub fn new() -> LorentzVector<T> {
        LorentzVector {
            t: T::default(),
            x: T::default(),
            y: T::default(),
            z: T::default(),
        }
    }

    pub fn from(t: T, x: T, y: T, z: T) -> LorentzVector<T> {
        LorentzVector { t, x, y, z }
    }

    pub fn from_vec(v: Vec<T>) -> LorentzVector<T> {
        let (t, x, y, z) = (v[0], v[1], v[2], v[3]);
        LorentzVector {
            t: t,
            x: x,
            y: y,
            z: z,
        }
    }

    pub fn square(&self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z - self.t * self.t
    }

    pub fn dot(&self, other: &LorentzVector<T>) -> T {
        self.x * other.x + self.y * other.y + self.z * other.z - self.t * other.t
    }
}

impl<'a, T: Field> Add<&'a LorentzVector<T>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    fn add(self, other: &'a LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t + other.t,
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<'a, T: Field> Sub<&'a LorentzVector<T>> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    fn sub(self, other: &'a LorentzVector<T>) -> LorentzVector<T> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<'a, T: Field> Mul<T> for &'a LorentzVector<T> {
    type Output = LorentzVector<T>;

    fn mul(self, other: T) -> LorentzVector<T> {
        LorentzVector {
            t: self.t * other,
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl<'a> Sub<&'a LorentzVector<f64>> for &'a LorentzVector<Complex> {
    type Output = LorentzVector<Complex>;

    fn sub(self, other: &'a LorentzVector<f64>) -> LorentzVector<Complex> {
        LorentzVector {
            t: self.t - other.t,
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.y,
        }
    }
}

impl LorentzVector<f64> {
    pub fn spatial_distance(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
}
