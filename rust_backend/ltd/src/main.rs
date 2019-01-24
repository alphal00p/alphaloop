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
use vector::{Field, LorentzVector, RealField};

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
    pub loop_momenta: ArrayVec<[(usize, bool); MAX_LOOP]>,
    pub q_and_mass: Vec<(LorentzVector<f64>, f64)>,
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

    // Evaluate loop line with a Feynman prescription
    pub fn evaluate_feynman(&self, loop_momenta_values: &[Vector3<Complex>]) -> Complex {
        let mut loop_momenta_eval = Vector3::zeros();
        for &(index, sign) in self.loop_momenta.iter() {
            if sign {
                loop_momenta_eval += loop_momenta_values[index];
            } else {
                loop_momenta_eval -= loop_momenta_values[index];
            }
        }

        // TODO: what are we supposed to do?
        let mut denom = Complex::new(1., 0.);
        for (q, m) in self.q_and_mass.iter() {
            let p = loop_momenta_eval + LoopLine::spatial_vec(&q.to_complex(true));
            denom *= p.dot(&p) + m * m;
        }

        1. / denom
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

    #[inline]
    fn compute_lambda_factor<T: Field + RealField>(x: T, y: T) -> T {
        if x * 2. < y {
            y * 0.25
        } else if y < 0. {
            x - y * 0.5
        } else {
            x - y * 0.25
        }
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
        let aij = 0.01;

        // TODO: we are deforming too much
        // TODO: add lambda per kappa_i?
        let mut kappa = Vector3::zeros();
        for (i, ((d, _m), dplus)) in self.q_and_mass.iter().zip(ki_plus.iter()).enumerate() {
            for (j, ((p, _m), pplus)) in self.q_and_mass.iter().zip(ki_plus.iter()).enumerate() {
                if i != j {
                    let inv = LoopLine::inv_g_d_dual(*dplus, *pplus, Dual4::from_real(p[0] - d[0]));
                    let e = (-inv * inv / (aij * e_cm.powi(4))).exp();
                    kappa += -(dirs[i] + dirs[j]) * e;
                }
            }
        }

        // now determine the global lambda scaling
        let mut lambda_sq = Dual4::from_real(1.);
        for (i, (((q_cut, _mass_cut), ki_s_cut), &ki_plus_cut)) in self
            .q_and_mass
            .iter()
            .zip(ki_spatial.iter())
            .zip(ki_plus.iter())
            .enumerate()
        {
            for (j, (((q_dual, _mass_dual), ki_s_dual), &ki_plus_dual)) in self
                .q_and_mass
                .iter()
                .zip(ki_spatial.iter())
                .zip(ki_plus.iter())
                .enumerate()
            {
                let (a, b, c) = if i == j {
                    // "propagator" coming from delta function cut
                    let a = -kappa.dot(&kappa);
                    let b = ki_s_cut.dot(&kappa) * 2.;
                    let c = ki_plus_cut * ki_plus_cut;
                    (a, b, c)
                } else {
                    let e = (q_dual[0] - q_cut[0]).powi(2);
                    let c1 = ki_plus_cut * ki_plus_cut + e - ki_plus_dual * ki_plus_dual;
                    let b1 = (ki_s_cut - ki_s_dual).dot(&kappa) * 2.0;
                    let a = -b1 * b1 + kappa.dot(&kappa) * e * 4.;
                    let b = b1 * c1 * 2. - ki_s_cut.dot(&kappa) * e * 8.;
                    let c = c1 * c1 - ki_plus_cut * ki_plus_cut * e * 4.;
                    (a, b, c)
                };

                let x = b * b / (a * a * 4.);
                let y = -c / a;
                let ls = LoopLine::compute_lambda_factor(x, y);

                if ls < lambda_sq {
                    lambda_sq = ls;
                }
            }
        }

        kappa *= lambda_sq.sqrt();

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

    fn deform_two_loops(&self, x: &[f64]) -> (Vector3<f64>, Vector3<f64>, Complex) {
        (Vector3::zeros(), Vector3::zeros(), Complex::new(1., 0.))
    }

    #[inline]
    pub fn evaluate(&self, x: &[f64], deform: bool) -> Complex {
        match x.len() / 3 {
            1 => {
                let (l_space, jac) = self.parameterize(x);
                let l = Vector3::from_row_slice(&l_space);

                let loop_line = self.loop_lines.first().unwrap();
                let (kappa, def_jac) = if deform {
                    loop_line.deform(self.e_cm, &[l])
                } else {
                    (Vector3::zeros(), Complex::new(1., 0.))
                };

                let k = l.map(|x| Complex::new(x, 0.)) + kappa.map(|x| Complex::new(0., x));
                let res = loop_line.evaluate(&[k]);
                res * jac * def_jac
            }
            2 => {
                let (l_space, jac1) = self.parameterize(&x[..3]);
                let l1 = Vector3::from_row_slice(&l_space);
                let (l_space, jac2) = self.parameterize(&x[3..]);
                let l2 = Vector3::from_row_slice(&l_space);

                let (kappa1, kappa2, def_jac) = if deform {
                    self.deform_two_loops(x)
                } else {
                    (Vector3::zeros(), Vector3::zeros(), Complex::new(1., 0.))
                };

                let k1 = l1.map(|x| Complex::new(x, 0.)) + kappa1.map(|x| Complex::new(0., x));
                let k2 = l2.map(|x| Complex::new(x, 0.)) + kappa2.map(|x| Complex::new(0., x));

                // by convention, loop line 1 only has k1 and loop line 2 only has k2
                let (ll1, ll2, ll3) = (
                    &self.loop_lines[0],
                    &self.loop_lines[1],
                    &self.loop_lines[2],
                );
                assert!(
                    ll1.loop_momenta.len() == 1
                        && ll1.loop_momenta[0] == (0, true)
                        && ll2.loop_momenta.len() == 1
                        && ll2.loop_momenta[0] == (1, true)
                );

                let res = ll1.evaluate(&[k1, k2]) * ll2.evaluate(&[k1, k2]) * ll3.evaluate(&[k1, k2]) // triple cut
                + ll1.evaluate(&[-k1, k2]) * ll2.evaluate(&[k1, k2]) * ll3.evaluate(&[k1, k2]) // triple cut with sign flip
                + ll1.evaluate(&[k1, k2]) * ll2.evaluate(&[k1, k2]) * ll3.evaluate_feynman(&[k1, k2]) // dual cut
                + ll1.evaluate_feynman(&[-k1, k2]) * ll2.evaluate(&[k1, k2]) * ll3.evaluate(&[k1, k2]) // dual cut with (unnecessary) sign flip
                + ll1.evaluate(&[-k1, k2]) * ll2.evaluate_feynman(&[k1, k2]) * ll3.evaluate(&[k1, k2]); // dual cut with sign flip

                res * jac1 * jac2 * def_jac
            }
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

    let (loops, e_cm, loop_lines) = topologies::create_topology("double-triangle");

    let r = ci.vegas(
        3 * loops,
        1,
        CubaVerbosity::Progress,
        0,
        UserData {
            ltd: vec![LTD::new(e_cm, loop_lines.clone()); cores + 1],
        },
    );
    println!("{:#?}", r);
}
