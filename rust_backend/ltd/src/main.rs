extern crate arrayvec;
extern crate cuba;
extern crate dual_num;
extern crate nalgebra as na;
extern crate num;
extern crate num_traits;
extern crate vector;
use arrayvec::ArrayVec;
use cuba::{CubaIntegrator, CubaVerbosity};
use dual_num::{DualN, U4};
use na::Vector3;
use num_traits::Float;
use std::f64::consts::PI;
use vector::{Field, LorentzVector};

mod topologies;
mod utils;

const MAX_DIM: usize = 4;
const MAX_PROP: usize = 6;
const MAX_LOOP: usize = 2;

type Dual4 = DualN<f64, U4>;
type Complex = num::Complex<f64>;

/// A loop line are all propagators with the same loop-momenta in it.
/// It is the basic component of LTD.
#[derive(Debug, Clone)]
pub struct LoopLine {
    loop_momenta: ArrayVec<[(usize, bool); MAX_LOOP]>,
    q_and_mass: Vec<(LorentzVector<f64>, f64)>,
}

impl LoopLine {
    pub fn new(
        loop_momenta: &[(usize, bool)],
        q_and_mass: Vec<(LorentzVector<f64>, f64)>,
    ) -> LoopLine {
        let lt: ArrayVec<[(usize, bool); MAX_LOOP]> = loop_momenta.iter().cloned().collect();
        LoopLine {
            loop_momenta: lt,
            q_and_mass,
        }
    }

    #[inline]
    fn spatial_vec<T: Field + 'static>(a: &LorentzVector<T>) -> Vector3<T> {
        Vector3::new(a.x, a.y, a.z)
    }

    #[inline]
    pub fn q_plus(qi: &Vector3<Complex>, mass: f64) -> Complex {
        (qi.dot(qi) + mass * mass).sqrt()
    }

    #[inline]
    pub fn inv_g_d(q_i0p: Complex, q_j0p: Complex, k_ji0: Complex) -> Complex {
        (q_i0p + k_ji0) * (q_i0p + k_ji0) - q_j0p * q_j0p
    }

    #[inline]
    pub fn q_plus_dual(qi: &Vector3<Dual4>, mass: f64) -> Dual4 {
        (qi.dot(qi) + mass * mass).sqrt()
    }

    #[inline]
    pub fn inv_g_d_dual(q_i0p: Dual4, q_j0p: Dual4, k_ji0: Dual4) -> Dual4 {
        (q_i0p + k_ji0) * (q_i0p + k_ji0) - q_j0p * q_j0p
    }

    /// Apply LTD to a loop line
    pub fn evaluate(&self, loop_momenta_values: &[Vector3<Complex>]) -> Complex {
        let mut loop_momenta_eval = Vector3::zeros();
        for &(index, sign) in self.loop_momenta.iter() {
            if sign {
                loop_momenta_eval += loop_momenta_values[index];
            } else {
                loop_momenta_eval -= loop_momenta_values[index];
            }
        }

        let ki_spatial: ArrayVec<[_; MAX_PROP]> = self
            .q_and_mass
            .iter()
            .map(|(x, _mass)| LoopLine::spatial_vec(&x.to_complex(true)) + loop_momenta_eval)
            .collect();
        let ki_plus: ArrayVec<[_; MAX_PROP]> = ki_spatial
            .iter()
            .zip(&self.q_and_mass)
            .map(|(x, (_q, m))| LoopLine::q_plus(x, *m))
            .collect();

        let mut res = Complex::new(0., 0.);
        for index in 0..self.q_and_mass.len() {
            let mut denom = Complex::new(1., 0.);
            for (i, ((q, _mass), qplus)) in self.q_and_mass.iter().zip(ki_plus.iter()).enumerate() {
                if i != index {
                    denom *= LoopLine::inv_g_d(
                        ki_plus[index],
                        *qplus,
                        Complex::new(q[0] - self.q_and_mass[index].0[0], 0.),
                    );
                } else {
                    denom *= 2. * ki_plus[index];
                }
            }
            res += Complex::new(0., -2. * PI) / denom;
        }

        let factor = Complex::new(0., -1. / (2. * PI).powi(loop_momenta_eval.len() as i32 + 1));
        factor * res
    }

    /// Compute the deformation vector and the deformation Jacobian wrt the loop momentum
    /// in the loop line
    pub fn deform(
        &self,
        e_cm: f64,
        loop_momenta_values: &[Vector3<f64>],
    ) -> (Vector3<f64>, Complex) {
        let mut loop_momenta_eval = Vector3::zeros();
        for &(index, sign) in self.loop_momenta.iter() {
            if sign {
                loop_momenta_eval += loop_momenta_values[index];
            } else {
                loop_momenta_eval -= loop_momenta_values[index];
            }
        }

        let mut dual_x = loop_momenta_eval.map(|y| Dual4::from_real(y));
        for (i, x) in dual_x.iter_mut().enumerate() {
            x[i + 1] = 1.0;
        }

        let ki_spatial: ArrayVec<[_; MAX_PROP]> = self
            .q_and_mass
            .iter()
            .map(|(x, _mass)| x.map(|y| Dual4::from_real(y)))
            .map(|x| LoopLine::spatial_vec(&x) + dual_x)
            .collect();
        let ki_plus: ArrayVec<[_; MAX_PROP]> = ki_spatial
            .iter()
            .zip(&self.q_and_mass)
            .map(|(x, (_q, m))| LoopLine::q_plus_dual(x, *m))
            .collect();

        let dirs: ArrayVec<[_; MAX_PROP]> =
            ki_spatial.iter().map(|x| x / x.dot(x).sqrt()).collect();
        let lambda = 0.06;
        let aij = 0.01;

        // TODO: we are deforming too much
        let mut kappa = Vector3::zeros();
        for (i, ((d, _m), dplus)) in self.q_and_mass.iter().zip(ki_plus.iter()).enumerate() {
            for (j, ((p, _m), pplus)) in self.q_and_mass.iter().zip(ki_plus.iter()).enumerate() {
                if i != j {
                    let inv = LoopLine::inv_g_d_dual(*dplus, *pplus, Dual4::from_real(p[0] - d[0]));
                    let e = (-inv * inv / (aij * e_cm.powi(4))).exp();
                    kappa += -(dirs[i] + dirs[j]) * e * Dual4::from_real(lambda);
                }
            }
        }

        let mut jac_mat = [[Complex::new(0., 0.); 3]; 3];
        for (i, k) in kappa.iter().enumerate() {
            jac_mat[i][i] = Complex::new(1., 0.);
            for (j, kd) in k.iter().skip(1).enumerate() {
                jac_mat[i][j] += Complex::new(0., *kd);
            }
        }

        let jac = utils::determinant(&jac_mat);
        (kappa.map(|x| x.real()), jac)
    }
}

#[derive(Debug, Clone)]
struct LTD {
    e_cm: f64,
    loop_lines: Vec<LoopLine>,
}

impl LTD {
    pub fn new(e_cm: f64, loop_lines: Vec<LoopLine>) -> LTD {
        LTD { e_cm, loop_lines }
    }

    /// Map from unit hypercube to infinite hypercube in N-d
    pub fn parameterize(&self, x: &[f64]) -> (ArrayVec<[f64; MAX_DIM]>, f64) {
        let mut jac = 1.;
        let radius = self.e_cm * x[0] / (1. - x[0]); // in [0,inf)
        jac *= (self.e_cm + radius).powi(2) / self.e_cm;
        assert!(x.len() > 1);
        let phi = 2. * PI * x[1];
        jac *= 2. * PI;

        match x.len() {
            2 => {
                let mut l_space = ArrayVec::new();
                l_space.push(radius * phi.cos());
                l_space.push(radius * phi.sin());
                jac *= radius;
                (l_space, jac)
            }
            3 => {
                let cos_theta = -1. + 2. * x[2];
                jac *= 2.;
                let sin_theta = (1. - cos_theta * cos_theta).sqrt();
                let mut l_space = ArrayVec::new();
                l_space.push(radius * sin_theta * phi.cos());
                l_space.push(radius * sin_theta * phi.sin());
                l_space.push(radius * cos_theta);
                jac *= radius * radius; // spherical coord
                (l_space, jac)
            }
            x => unimplemented!("Unknown dimension {}", x),
        }
    }

    #[inline]
    pub fn evaluate(&self, x: &[f64], deform: bool) -> Complex {
        match self.loop_lines.len() {
            1 => {
                let (l_space, jac) = self.parameterize(x);
                let l_vec = Vector3::from_row_slice(&l_space);

                let loop_line = self.loop_lines.first().unwrap();
                let (kappa, def_jac) = if deform {
                    loop_line.deform(self.e_cm, &[l_vec])
                } else {
                    (Vector3::zeros(), Complex::new(1., 0.))
                };

                let k = l_vec.map(|x| Complex::new(x, 0.)) + kappa.map(|x| Complex::new(0., x));
                let res = loop_line.evaluate(&[k]);
                res * jac * def_jac
            }
            2 => unimplemented!("Two-loop LTD will be added soon!"),
            x => unimplemented!("LTD for {} loops is not supported yet", x),
        }
    }
}

#[derive(Debug)]
struct UserData {
    ltd: Vec<LTD>,
}

#[inline(always)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    _nvec: usize,
    core: i32,
) -> Result<(), &'static str> {
    let res = user_data.ltd[(core + 1) as usize].evaluate(x, true);
    f[0] = res.re;
    Ok(())
}

fn main() {
    let mut ci = CubaIntegrator::new(integrand);
    let cores = 4;
    ci.set_mineval(10)
        .set_nstart(1000000)
        .set_maxeval(100000000)
        .set_epsabs(0.)
        .set_epsrel(1e-15)
        .set_seed(1)
        .set_cores(cores, 1000);

    let (e_cm, loop_lines) = topologies::create_topology("P1");

    let r = ci.vegas(
        3,
        1,
        CubaVerbosity::Progress,
        0,
        UserData {
            ltd: vec![LTD::new(e_cm, loop_lines.clone()); cores + 1],
        },
    );
    println!("{:#?}", r);
}
