use dual_num::DualN;
use f128::f128;
use itertools::Itertools;
use num::Complex;
use num_traits::{Float, Num, NumAssign, NumCast};
use num_traits::{Inv, One, Zero};
use std::ops::Neg;
use vector::RealNumberLike;
use vector::{Field, LorentzVector};
use FloatLike;
use MAX_LOOP;

const MAX_DIMENSION: usize = MAX_LOOP * 3;

/// Format a mean ± sdev as mean(sdev) with the correct number of digits.
/// Based on the Python package gvar.
pub fn format_uncertainty(mean: f64, sdev: f64) -> String {
    fn ndec(x: f64, offset: usize) -> i32 {
        let mut ans = (offset as f64 - x.log10()) as i32;
        if ans > 0 && x * 10.0.powi(ans) >= [0.5, 9.5, 99.5][offset] {
            ans -= 1;
        }
        if ans < 0 {
            0
        } else {
            ans
        }
    }
    let v = mean;
    let dv = sdev.abs();

    // special cases
    if v.is_nan() || dv.is_nan() {
        format!("{:e} ± {:e}", v, dv)
    } else if dv.is_infinite() {
        format!("{:e} ± inf", v)
    } else if v == 0. && (dv >= 1e5 || dv < 1e-4) {
        if dv == 0. {
            "0(0)".to_owned()
        } else {
            let e = format!("{:.1e}", dv);
            let mut ans = e.split('e');
            let e1 = ans.next().unwrap();
            let e2 = ans.next().unwrap();
            "0.0(".to_owned() + e1 + ")e" + e2
        }
    } else if v == 0. {
        if dv >= 9.95 {
            format!("0({:.0})", dv)
        } else if dv >= 0.995 {
            format!("0.0({:.1})", dv)
        } else {
            let ndecimal = ndec(dv, 2);
            format!(
                "{:.*}({:.0})",
                ndecimal as usize,
                v,
                dv * 10.0.powi(ndecimal)
            )
        }
    } else if dv == 0. {
        let e = format!("{:e}", v);
        let mut ans = e.split('e');
        let e1 = ans.next().unwrap();
        let e2 = ans.next().unwrap();
        if e2 != "0" {
            e1.to_owned() + "(0)e" + e2
        } else {
            e1.to_owned() + "(0)"
        }
    } else if dv > 1e4 * v.abs() {
        format!("{:.1e} ± {:.2e}", v, dv)
    } else if v.abs() >= 1e6 || v.abs() < 1e-5 {
        // exponential notation for large |self.mean|
        let exponent = v.abs().log10().floor();
        let fac = 10.0.powf(exponent);
        let mantissa = format_uncertainty(v / fac, dv / fac);
        let e = format!("{:.0e}", fac);
        let mut ee = e.split('e');
        mantissa + "e" + ee.nth(1).unwrap()
    }
    // normal cases
    else if dv >= 9.95 {
        if v.abs() >= 9.5 {
            format!("{:.0}({:.0})", v, dv)
        } else {
            let ndecimal = ndec(v.abs(), 1);
            format!("{:.*}({:.*})", ndecimal as usize, v, ndecimal as usize, dv)
        }
    } else if dv >= 0.995 {
        if v.abs() >= 0.95 {
            format!("{:.1}({:.1})", v, dv)
        } else {
            let ndecimal = ndec(v.abs(), 1);
            format!("{:.*}({:.*})", ndecimal as usize, v, ndecimal as usize, dv)
        }
    } else {
        let ndecimal = ndec(v.abs(), 1).max(ndec(dv, 2));
        format!(
            "{:.*}({:.0})",
            ndecimal as usize,
            v,
            dv * 10.0.powi(ndecimal)
        )
    }
}

pub trait Signum {
    fn multiply_sign(&self, sign: i8) -> Self;
}

impl Signum for f128 {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> f128 {
        match sign {
            1 => self.clone(),
            0 => f128::zero(),
            -1 => self.neg(),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl Signum for f64 {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> f64 {
        match sign {
            1 => self.clone(),
            0 => f64::zero(),
            -1 => self.neg(),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: Num + Neg<Output = T> + Copy> Signum for Complex<T> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> Complex<T> {
        match sign {
            1 => *self,
            0 => Complex::zero(),
            -1 => Complex::new(-self.re, -self.im),
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<U: dual_num::Dim + dual_num::DimName, T: FloatLike> Signum for DualN<T, U>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    #[inline]
    fn multiply_sign(&self, sign: i8) -> DualN<T, U> {
        match sign {
            1 => *self,
            0 => DualN::zero(),
            -1 => -*self,
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

impl<T: Field> Signum for LorentzVector<T> {
    #[inline]
    fn multiply_sign(&self, sign: i8) -> LorentzVector<T> {
        match sign {
            1 => *self,
            0 => LorentzVector::default(),
            -1 => -self,
            _ => unreachable!("Sign should be -1,0,1"),
        }
    }
}

#[inline]
/// Invert with better precision
pub fn finv<T: Float>(c: Complex<T>) -> Complex<T> {
    let norm = c.norm();
    c.conj() / norm / norm
}

#[inline]
pub fn powi<T: Float + NumAssign>(c: Complex<T>, n: usize) -> Complex<T> {
    let mut c1 = Complex::<T>::one();
    for _ in 0..n {
        c1 *= c;
    }
    c1
}

pub fn evaluate_signature<T: RealNumberLike>(
    signature: &[i8],
    momenta: &[LorentzVector<T>],
) -> LorentzVector<T> {
    let mut momentum = LorentzVector::default();
    for (&sign, mom) in signature.iter().zip_eq(momenta) {
        if sign != 0 {
            momentum += mom.multiply_sign(sign);
        }
    }

    momentum
}

/// Calculate the determinant of any complex-valued input matrix using LU-decomposition.
/// Original C-code by W. Gong and D.E. Soper.
pub fn determinant<T: Float + RealNumberLike>(
    bb: &Vec<Complex<T>>,
    dimension: usize,
) -> Complex<T> {
    // Define matrix related variables.
    let mut determinant = Complex::new(T::one(), T::zero());
    let mut indx = [0; MAX_DIMENSION];
    let mut d = 1; // initialize parity parameter

    // Inintialize the matrix to be decomposed with the transferred matrix b.
    let mut aa = bb.clone();

    // Define parameters used in decomposition.
    let mut imax = 0;
    let mut flag = 1;
    let mut dumc;
    let mut sum;

    let mut aamax;
    let mut dumr;
    let mut vv = [T::zero(); MAX_DIMENSION];

    // Get the implicit scaling information.
    for i in 0..dimension {
        aamax = T::zero();
        for j in 0..dimension {
            let r = aa[i * dimension + j].norm_sqr();
            if r > aamax {
                aamax = r;
            }
        }
        // Set a flag to check if the determinant is zero.
        if aamax.is_zero() {
            flag = 0;
        }
        // Save the scaling.
        vv[i] = aamax.inv();
    }
    if flag == 1 {
        for j in 0..dimension {
            for i in 0..j {
                sum = aa[i * dimension + j];
                for k in 0..i {
                    sum = sum - aa[i * dimension + k] * aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
            }
            //Initialize for the search for largest pivot element.
            aamax = T::zero();
            for i in j..dimension {
                sum = aa[i * dimension + j];
                for k in 0..j {
                    sum = sum - aa[i * dimension + k] * aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
                // Figure of merit for the pivot.
                dumr = vv[i] * sum.norm_sqr();
                // Is it better than the best so far?
                if dumr >= aamax {
                    imax = i;
                    aamax = dumr;
                }
            }
            // See if we need to interchange rows.
            if j != imax {
                for k in 0..dimension {
                    dumc = aa[imax * dimension + k];
                    aa[imax * dimension + k] = aa[j * dimension + k];
                    aa[j * dimension + k] = dumc;
                }
                // Change the parity of d.
                d = -d;
                // Interchange the scale factor.
                vv[imax] = vv[j];
            }
            indx[j] = imax;
            if j + 1 != dimension {
                dumc = aa[j * dimension + j].inv();
                for i in j + 1..dimension {
                    aa[i * dimension + j] = aa[i * dimension + j] * dumc;
                }
            }
        }
    }
    // Calculate the determinant using the decomposed matrix.
    if flag == 0 {
        determinant = Complex::default();
    } else {
        // Multiply the diagonal elements.
        for diagonal in 0..dimension {
            determinant = determinant * aa[diagonal * dimension + diagonal];
        }
        determinant = determinant * <T as NumCast>::from(d).unwrap();
    }
    determinant
}

pub fn next_combination_with_replacement(state: &mut [usize], max_entry: usize) -> bool {
    for i in (0..state.len()).rev() {
        if state[i] < max_entry {
            state[i] += 1;
            for j in i + 1..state.len() {
                state[j] = state[i]
            }
            return true;
        }
    }
    false
}
