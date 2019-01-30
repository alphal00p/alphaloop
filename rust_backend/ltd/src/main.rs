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

    // Evaluate loop line with Lorentz vectors and an optional cut index
    pub fn evaluate_4d(
        &self,
        cut_index: i32,
        loop_momenta_values: &[LorentzVector<Complex>],
    ) -> Complex {
        let mut loop_momenta_eval = LorentzVector::new();
        for &(index, sign) in self.loop_momenta.iter() {
            if sign {
                loop_momenta_eval = loop_momenta_eval + loop_momenta_values[index];
            } else {
                loop_momenta_eval = loop_momenta_eval - loop_momenta_values[index];
            }
        }

        let mut denom = Complex::new(1., 0.);
        for (i, (q, m)) in self.q_and_mass.iter().enumerate() {
            if i as i32 != cut_index {
                denom *= (loop_momenta_eval + q).square() - m * m;
            } else {
                // FIXME: what to do for triple cut?
                // It should be 1/(l_vec+ - l_vec-) where
                // l_vec+- is the solution of the vector part solving the triple cut
                denom *= 2. * (loop_momenta_eval.t + q.t);
            }
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

    fn deform_two_loops(
        &self,
        ll: &[LoopLine],
        loop_momenta: &[LorentzVector<f64>],
    ) -> (LorentzVector<f64>, LorentzVector<f64>, Complex) {
        (
            LorentzVector::new(),
            LorentzVector::new(),
            Complex::new(1., 0.),
        )
    }

    /// Parameterize, taking deltas into account
    pub fn parameterize_with_delta(
        &self,
        ll: &[LoopLine],
        cut_indices: &[i32],
        x: &[f64],
    ) -> (ArrayVec<[LorentzVector<f64>; MAX_LOOP]>, f64) {
        assert!(cut_indices.len() == 3); // only two loops for now

        // by convention, loop line 1 only has +-k1 and loop line 2 only has k2
        // and k3=k1+k2
        let (ll1, ll2, ll3) = (&ll[0], &ll[1], &ll[2]);
        assert!(
            ll1.loop_momenta.len() == 1
                && ll1.loop_momenta[0].0 == 0
                && ll2.loop_momenta.len() == 1
                && ll2.loop_momenta[0] == (1, true)
                && ll3.loop_momenta[0] == (0, true)
                && ll3.loop_momenta[1] == (1, true)
        );

        // generate the spatial parts of k and l
        // for the triple cut we will fix an angle later
        // x: r_k1, phi_k1, theta_k1, phi_k2, theta_k2, r_k2
        let radius = self.e_cm * x[0] / (1. - x[0]); // in [0,inf)
        let mut jac = (self.e_cm + radius).powi(2) / self.e_cm * 2. * PI * 2. * radius * radius;
        let phi = 2. * PI * x[1];
        let cos_theta = -1. + 2. * x[2];
        let sin_theta = (1. - cos_theta * cos_theta).sqrt();
        let mut k1 = LorentzVector::from_args(
            0.,
            radius * sin_theta * phi.cos(),
            radius * sin_theta * phi.sin(),
            radius * cos_theta,
        );

        let radius2 = if cut_indices.iter().all(|&x| x > -1) {
            // for the triple cut case we fix the radius
            1.0
        } else {
            let radius2 = self.e_cm * x[5] / (1. - x[5]); // in [0,inf)
            jac *= (self.e_cm + radius2).powi(2) / self.e_cm * radius2 * radius2;
            radius2
        };
        jac *= 2. * PI * 2.;
        let phi2 = 2. * PI * x[3];
        let cos_theta2 = -1. + 2. * x[4];
        let sin_theta2 = (1. - cos_theta2 * cos_theta2).sqrt();
        let mut k2 = LorentzVector::from_args(
            0.,
            radius2 * sin_theta2 * phi2.cos(),
            radius2 * sin_theta2 * phi2.sin(),
            radius2 * cos_theta2,
        );

        // FIXME: we assume k3=k1+k2 here!
        match (cut_indices[0], cut_indices[1], cut_indices[2]) {
            (-1, c2, c3) if c2 > -1 && c3 > -1 => {
                let (q, m) = ll2.q_and_mass[c2 as usize];
                k2[0] = -q.t + ((k2 + q).spatial_squared() + m * m).sqrt();
                let (q, m) = ll3.q_and_mass[c3 as usize];
                k1[0] = -k2.t - q.t + ((k1 + k2 + q).spatial_squared() + m * m).sqrt();
            }
            (c1, -1, c3) if c1 > -1 && c3 > -1 => {
                let (q, m) = ll1.q_and_mass[c1 as usize];
                k1[0] = -q.t + ((k1 + q).spatial_squared() + m * m).sqrt();
                let (q, m) = ll3.q_and_mass[c3 as usize];
                k2[0] = -k1.t - q.t + ((k1 + k2 + q).spatial_squared() + m * m).sqrt();
            }
            (c1, c2, -1) if c1 > -1 && c2 > -1 => {
                let (q, m) = ll1.q_and_mass[c1 as usize];
                k1[0] = -q.t + ((k1 + q).spatial_squared() + m * m).sqrt();
                let (q, m) = ll2.q_and_mass[c2 as usize];
                k2[0] = -q.t + ((k2 + q).spatial_squared() + m * m).sqrt();
            }
            (c1, c2, c3) if c1 > -1 && c2 > -1 && c3 > -1 => {
                // triple cut case
                let (q1, m1) = ll1.q_and_mass[c1 as usize];
                let (q2, m2) = ll2.q_and_mass[c2 as usize];
                let (q3, m3) = ll3.q_and_mass[c3 as usize];

                k1[0] = -q1.t + ((k1 + q1).spatial_squared() + m1 * m1).sqrt();

                let r = k1 + q3;
                let k = q2.square() - r.square() - m2 * m2 - m3 * m3;

                // determine new radius
                let k2r = k2.spatial_dot(&r);
                let a = (k2r / r.t).powi(2) - 1.;
                let b = 2. * k2r * q2.t / r.t - k * k2r / r.t / r.t - 2. * k2.spatial_dot(&q2);
                let c = k * k / (4. * r.t * r.t) - k * q2.t / r.t + q2.square() - m2 * m2;
                let d = b * b - 4. * a * c;
                if d < 0. {
                    return (ArrayVec::from([k1, k2]), 0.);
                }
                let radius2 = (-b + d.sqrt()) / (2. * a); // new radius
                jac *= radius2 * radius2;

                k2.x *= radius2;
                k2.y *= radius2;
                k2.z *= radius2;
                k2.t = (2. * k2.spatial_dot(&r) - k) / (2. * r.t);
            }
            x => panic!("Not enough cuts: {:?}", x),
        }

        (ArrayVec::from([k1, k2]), jac)
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
                // all LTD cuts at two-loops for k, l, and k+l
                // 1 means positive cut, 0 means no cut and -1 negative cut
                let cut_structures = [[1, 1, 0], [0, 1, 1], [-1, 0, 1], [1, 1, 1], [-1, 1, 1]];

                let mut result = Complex::default();
                for o in &cut_structures {
                    // set the loop line directions for this configuration
                    let mut ll = self.loop_lines.clone();
                    for (i, &loopline_cut) in o.iter().enumerate() {
                        if loopline_cut < 0 {
                            for (_, sign) in &mut ll[i].loop_momenta {
                                *sign = !*sign;
                            }
                        }
                    }

                    // TODO: cartesian product iterator
                    let mut cut_indices = [-1; 3];
                    for (i, _c1) in ll[0].q_and_mass.iter().enumerate() {
                        for (j, _c2) in ll[1].q_and_mass.iter().enumerate() {
                            for (k, _c3) in ll[2].q_and_mass.iter().enumerate() {
                                if i > 0 && o[0] == 0 || j > 0 && o[1] == 0 || k > 0 && o[2] == 0 {
                                    // skip configurations where we cut a loop line that we shouldn't cut
                                    continue;
                                }

                                cut_indices[0] = if o[0] == 0 { -1 } else { i as i32 };
                                cut_indices[1] = if o[1] == 0 { -1 } else { j as i32 };
                                cut_indices[2] = if o[2] == 0 { -1 } else { k as i32 };

                                // parameterize
                                let (momenta, par_jac) =
                                    self.parameterize_with_delta(&ll, &cut_indices, x);

                                if par_jac == 0. {
                                    continue;
                                }

                                // deform
                                let (kappa1, kappa2, def_jac) =
                                    self.deform_two_loops(&ll, &momenta);
                                let k1 = momenta[0].map(|x| Complex::new(x, 0.))
                                    + kappa1.map(|x| Complex::new(0., x));
                                let k2 = momenta[1].map(|x| Complex::new(x, 0.))
                                    + kappa2.map(|x| Complex::new(0., x));

                                // evaluate
                                let mut res = par_jac * def_jac;
                                for (l, &c) in ll.iter().zip(cut_indices.iter()) {
                                    res *= l.evaluate_4d(c, &[k1, k2]);
                                }

                                if !res.is_finite() {
                                    println!(
                                        "Bad result: k={}, l={}, par_jac={}, x={:?}, res={}",
                                        momenta[0], momenta[1], par_jac, x, res
                                    );
                                }

                                result += res;
                            }
                        }
                    }
                }

                // TODO: normalize
                result
            }
            x => unimplemented!("LTD for {} loops is not supported yet", x),
        }
    }
}

#[derive(Debug)]
struct UserData {
    ltd: Vec<LTD>,
    sample_count: usize,
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
    user_data.sample_count += 1;

    if res.re.is_finite() {
        f[0] = res.re;
    } else {
        println!("Bad point: {:?}", x);
        f[0] = 0.;
    }

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
            sample_count: 0,
        },
    );
    println!("{:#?}", r);
}
