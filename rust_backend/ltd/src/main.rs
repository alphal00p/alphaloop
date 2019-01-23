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

mod utils;

const MAX_DIM: usize = 4;

type Dual4 = DualN<f64, U4>;
type Complex = num::Complex<f64>;

#[derive(Debug, Clone)]
struct LTD {
    e_cm: f64,
    ext: Vec<LorentzVector<f64>>,
    qs: Vec<LorentzVector<f64>>,
    mass: Vec<f64>,
}

impl LTD {
    pub fn new(ext: &[LorentzVector<f64>], mass: &[f64]) -> LTD {
        let mut qs = vec![LorentzVector::default()];
        for e in &ext[1..] {
            let n = qs.last().unwrap() + e;
            qs.push(n);
        }

        let e_cm = (ext[0] + ext[1]).square().abs().sqrt();

        LTD {
            e_cm,
            ext: ext.to_vec(),
            qs: qs,
            mass: mass.to_vec(),
        }
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

    fn spatial_vec<T: Field + 'static>(a: &LorentzVector<T>) -> Vector3<T> {
        Vector3::new(a.x, a.y, a.z)
    }

    #[inline]
    pub fn q_plus(&self, qi: &Vector3<Complex>, mass: f64) -> Complex {
        (qi.dot(qi) + mass * mass).sqrt()
    }

    #[inline]
    pub fn q_plus_dual(&self, qi: &Vector3<Dual4>, mass: f64) -> Dual4 {
        (qi.dot(qi) + mass * mass).sqrt()
    }

    #[inline]
    pub fn inv_g_d(&self, q_i0p: Complex, q_j0p: Complex, k_ji0: Complex) -> Complex {
        (q_i0p + k_ji0) * (q_i0p + k_ji0) - q_j0p * q_j0p
    }

    #[inline]
    pub fn inv_g_d_dual(&self, q_i0p: Dual4, q_j0p: Dual4, k_ji0: Dual4) -> Dual4 {
        (q_i0p + k_ji0) * (q_i0p + k_ji0) - q_j0p * q_j0p
    }

    pub fn evaluate_dual(&self, index: usize, loop_momentum: &Vector3<Complex>) -> Complex {
        // FIXME: don't collect
        let ki_spatial: Vec<_> = self
            .qs
            .iter()
            .map(|x| LTD::spatial_vec(&x.to_complex(true)) + loop_momentum)
            .collect();
        let ki_plus: Vec<_> = ki_spatial
            .iter()
            .zip(&self.mass)
            .map(|(x, m)| self.q_plus(x, *m))
            .collect();

        let residue_factor = Complex::new(0., -2. * PI);
        let mut denom = Complex::new(1., 0.);
        for (i, (q, qplus)) in self.qs.iter().zip(ki_plus.iter()).enumerate() {
            if i != index {
                denom *= self.inv_g_d(
                    ki_plus[index],
                    *qplus,
                    Complex::new(q[0] - self.qs[index][0], 0.),
                );
            } else {
                denom *= 2. * ki_plus[index];
            }
        }

        residue_factor / denom
    }

    #[inline]
    pub fn evaluate_dual_integrand(&self, x: &[f64], deform: bool) -> Complex {
        let (l_space, jac) = self.parameterize(x);
        let l_vec = Vector3::from_row_slice(&l_space);

        let (kappa, def_jac) = if deform {
            self.deform(l_vec)
        } else {
            (Vector3::zeros(), Complex::new(1., 0.))
        };

        let k = l_vec.map(|x| Complex::new(x, 0.)) + kappa.map(|x| Complex::new(0., x));

        let mut res = Complex::new(0., 0.);
        for i in 0..self.qs.len() {
            res += self.evaluate_dual(i, &k);
        }
        let factor = Complex::new(0., -1. / (2. * PI).powi(x.len() as i32 + 1));
        factor * res * jac * def_jac
    }

    /// Compute the deformation vector and the deformation Jacobian
    pub fn deform(&self, x: Vector3<f64>) -> (Vector3<f64>, Complex) {
        let mut dual_x = x.map(|y| Dual4::from_real(y));
        for (i, x) in dual_x.iter_mut().enumerate() {
            x[i] = 1.0;
        }

        // TODO: store in struct
        let ki_spatial: Vec<_> = self
            .qs
            .iter()
            .map(|x| x.map(|y| Dual4::from_real(y)))
            .map(|x| LTD::spatial_vec(&x) + dual_x)
            .collect();

        let ki_plus: Vec<_> = ki_spatial
            .iter()
            .zip(&self.mass)
            .map(|(x, m)| self.q_plus_dual(&x, *m))
            .collect();

        let dirs: Vec<_> = ki_spatial.iter().map(|x| x / x.dot(x).sqrt()).collect();
        let lambda = 0.1;
        let aij = 0.1;

        // TODO: we are deforming too much
        let mut kappa = Vector3::zeros();
        for (i, (d, dplus)) in self.qs.iter().zip(ki_plus.iter()).enumerate() {
            for (j, (p, pplus)) in self.qs.iter().zip(ki_plus.iter()).enumerate() {
                if i != j {
                    let inv = self.inv_g_d_dual(*dplus, *pplus, Dual4::from_real(p[0] - d[0]));
                    let e = (-inv * inv / (aij * self.e_cm.powi(4))).exp();
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
    let res = user_data.ltd[(core + 1) as usize].evaluate_dual_integrand(x, true);
    f[0] = res.re;
    Ok(())
}

fn main() {
    let mut ci = CubaIntegrator::new(integrand);
    let cores = 1;
    ci.set_mineval(10)
        .set_nstart(1000000)
        .set_maxeval(10000000)
        .set_epsabs(0.)
        .set_epsrel(1e-15)
        .set_seed(1)
        .set_cores(cores, 1000);

    let topology = "P1";
    let (ext, mass) = match topology {
        "P1" => {
            // does not need deformation
            let p1 = LorentzVector::from_args(5.23923, -4.18858, 0.74966, -3.05669);
            let p2 = LorentzVector::from_args(6.99881, -2.93659, 5.03338, 3.87619);
            let m = 7.73358;
            (vec![p1, p2, -p1 - p2], vec![m, m, m])
        }
        "P5" => {
            // does not need deformation
            let p1 = LorentzVector::from_args(31.54872, -322.40325, 300.53015, -385.58013);
            let p2 = LorentzVector::from_args(103.90430, 202.00974, -451.27794, -435.12848);
            let p3 = LorentzVector::from_args(294.76653, 252.88958, 447.09194, 311.71630);
            let m = 4.68481;
            (vec![p1, p2, p3, -p1 - p2 - p3], vec![m; 4])
        }
        "P3" => {
            //P3 in https://arxiv.org/pdf/1510.00187.pdf
            let p1 = LorentzVector::from_args(10.51284, 6.89159, -7.40660, -2.85795);
            let p2 = LorentzVector::from_args(6.45709, 2.46635, 5.84093, 1.22257);
            let m = 0.52559;
            (vec![p1, p2, -p1 - p2], vec![m; 3])
        }
        x => unimplemented!("Unknown topology {}", x),
    };

    let r = ci.vegas(
        3,
        1,
        CubaVerbosity::Progress,
        0,
        UserData {
            ltd: vec![LTD::new(&ext, &mass); cores + 1],
        },
    );
    println!("{:#?}", r);
}
