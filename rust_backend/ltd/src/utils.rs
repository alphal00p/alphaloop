use dual_num::DualN;
use f128::f128;
use itertools::Itertools;
use num::Complex;
use num_traits::{Float, Num, NumAssign, NumCast};
use num_traits::{Inv, One, Zero};
use std::cmp::{Ord, Ordering};
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

/// Compare two slices, selecting on length first
pub fn compare_slice<T: Ord>(slice1: &[T], slice2: &[T]) -> Ordering {
    match slice1.len().cmp(&slice2.len()) {
        Ordering::Equal => (),
        non_eq => return non_eq,
    }

    let l = slice1.len();
    // Slice to the loop iteration range to enable bound check
    // elimination in the compiler
    let lhs = &slice1[..l];
    let rhs = &slice2[..l];

    for i in 0..l {
        match lhs[i].cmp(&rhs[i]) {
            Ordering::Equal => (),
            non_eq => return non_eq,
        }
    }

    Ordering::Equal
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

pub fn evaluate_signature<T: FloatLike>(
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

#[derive(Clone, Default)]
pub struct PolynomialReconstruction {
    matrices: Vec<Vec<f64>>,
    var_sample: Vec<f64>,
    sample_buffer: Vec<Complex<f64>>,
    scale: f64,
}

impl PolynomialReconstruction {
    /// Construct new polynomial reconstruction Vandermonde matrices for points
    /// `(0, 1, ...)*scale`.
    pub fn new(max_rank: usize, max_vars: usize, scale: f64) -> PolynomialReconstruction {
        // construct the inverse of the Vandermonde matrices
        let mut matrices = vec![];
        for rank in 0..=max_rank {
            let mat: Vec<f64> = (0..rank)
                .map(|r| {
                    (0..rank)
                        .map(|p| (r as f64 * scale).powi(p as i32))
                        .collect::<Vec<_>>()
                })
                .flatten()
                .collect();
            let m = nalgebra::DMatrix::from_row_slice(rank as usize, rank as usize, &mat);
            let inv = m.try_inverse().unwrap();
            matrices.push(inv.transpose().iter().cloned().collect());
        }

        PolynomialReconstruction {
            matrices,
            var_sample: vec![0.; max_vars],
            sample_buffer: Vec::with_capacity(max_rank),
            scale,
        }
    }

    pub fn reconstruct(
        &mut self,
        f: &mut dyn FnMut(&[f64]) -> Complex<f64>,
        rank: usize,
        n_vars: usize,
        result: &mut [Complex<f64>],
    ) {
        self.sample_buffer.clear();
        self.reconstruct_rec(f, rank, n_vars, 0, result);
    }

    fn reconstruct_rec(
        &mut self,
        f: &mut dyn FnMut(&[f64]) -> Complex<f64>,
        rank: usize,
        n_vars: usize,
        v: usize,
        result: &mut [Complex<f64>],
    ) {
        if v == n_vars {
            let r = f(&self.var_sample[..n_vars]);
            // return a constant polynomial
            for r in result.iter_mut() {
                *r = Complex::zero();
            }
            result[0] = r;
            return;
        }

        let sample_start = self.sample_buffer.len();
        for r in 0..rank {
            self.var_sample[v] = r as f64 * self.scale;
            self.reconstruct_rec(f, rank, n_vars, v + 1, result);
            self.sample_buffer.extend_from_slice(result);
        }

        // now do the reconstruction
        for r in result.iter_mut() {
            *r = Complex::zero();
        }

        for (v_pow, row) in self.matrices[rank].chunks(rank).enumerate() {
            for (r, sample) in row
                .iter()
                .zip(self.sample_buffer[sample_start..].chunks(result.len()))
            {
                // add the sample to the result
                for (powi, monomial) in sample.iter().enumerate() {
                    if !monomial.is_zero() {
                        // find the new place of the monomial
                        let new_index = powi + rank.pow((n_vars - 1 - v) as u32) * v_pow;
                        result[new_index] += *monomial * r;
                    }
                }
            }
        }

        self.sample_buffer.truncate(sample_start);
    }
}

pub mod test_utils {
    extern crate eyre;

    use color_eyre::Help;
    use num::Complex;
    use num_traits::cast::FromPrimitive;
    use num_traits::Float;
    use std::fmt::{Debug, LowerExp};
    use topologies::Topology;
    use {DeformationStrategy, Settings};

    pub fn get_test_topology(topology_file: &str, topology_name: &str) -> Topology {
        // Import settings and disable the deformation
        let mut settings = Settings::from_file("../../LTD/hyperparameters.yaml").unwrap();
        settings.general.deformation_strategy = DeformationStrategy::None;
        settings.general.multi_channeling = false;
        settings.integrator.dashboard = false;

        // Import topology from folder
        let mut topologies = Topology::from_file(topology_file, &settings).unwrap();
        settings.general.topology = topology_name.to_string();
        topologies
            .remove(&settings.general.topology)
            .ok_or_else(|| eyre::eyre!("Could not find topology {}", settings.general.topology))
            .suggestion("Check if this topology is in the specified topology file.")
            .unwrap()
    }

    pub fn numeriacal_eq<T: Float + FromPrimitive + LowerExp + Debug>(
        x: Complex<T>,
        y: Complex<T>,
    ) -> bool {
        let cutoff = T::epsilon().powf(T::from_f64(0.75).unwrap());
        let real_check = if x.re == y.re {
            true
        } else {
            ((x.re - y.re) / x.re) < cutoff
        };
        let imag_check = if x.im == y.im {
            true
        } else {
            ((x.im - y.im) / x.im) < cutoff
        };

        println!(
            "REAL: {:e} < {:e} ? {}",
            ((x.re - y.re) / x.re).abs(),
            cutoff,
            real_check
        );
        println!(
            "IMAG: {:e} < {:e} ? {}",
            ((x.im - y.im) / x.im).abs(),
            cutoff,
            imag_check
        );
        real_check && imag_check
    }
}
