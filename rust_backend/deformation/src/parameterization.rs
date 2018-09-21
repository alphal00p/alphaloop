use deformation::Deformer;
use num::traits::FloatConst;
use vector::LorentzVector;

pub enum ParameterizationMode {
    Log,
    Linear,
    Weinzierl,
}

pub struct Parameterizer {
    mode: ParameterizationMode,
    alpha: f64, // scaling parameter
    e_cm_sq: f64,
    channel: usize,
    qs: Vec<LorentzVector<f64>>,
}

impl Parameterizer {
    /// Map a vector in the unit hypercube to the infinite hypercube.
    /// Also compute the Jacobian.
    pub fn map(&self, mom: &LorentzVector<f64>) -> Result<(LorentzVector<f64>, f64), &str> {
        match self.mode {
            ParameterizationMode::Log => {
                let r = mom.map(|x| self.alpha * (x / (1.0 - x)).ln());

                let mut jacobian = 1.0;
                for i in 0..4 {
                    jacobian *= self.alpha / (mom[i] * (1.0 - mom[i]));
                }

                Ok((r, jacobian))
            }
            ParameterizationMode::Linear => {
                let r = mom.map(|x| self.alpha * (1.0 / (1.0 - x) - 1.0 / x));

                let mut jacobian = 1.0;
                for i in 0..4 {
                    jacobian *= self.alpha
                        * (1.0 / (mom[i] * mom[i]) + 1.0 / ((1.0 - mom[i]) * (1.0 - mom[i])));
                }

                Ok((r, jacobian))
            }
            ParameterizationMode::Weinzierl => {
                if self.channel == 0 {
                    return Err("A channel must be specified for Weinzierl's conformal map");
                }

                // TODO: cache x-independent block
                // TODO: check if this is the correct p
                let p = &self.qs[self.channel + 1] - &self.qs[self.channel];
                let p_abs = p.euclidean_distance();

                let cos_theta_1 = p.t / p_abs;
                let sin_theta_1 = (1.0 - cos_theta_1 * cos_theta_1).sqrt();

                let cos_theta_2 = if p.x == 0.0 {
                    0.0
                } else {
                    1.0 / (1.0 + (p.y * p.y + p.z * p.z).sqrt() / (p.x * p.x))
                };
                let sin_theta_2 = if p.x == 0.0 {
                    1.0
                } else {
                    (p.y * p.y + p.z * p.z).sqrt() / p.x * cos_theta_2
                };

                let cos_phi_3 = if p.y == 0.0 {
                    0.0
                } else {
                    1.0 / (1.0 + p.z * p.z / (p.y * p.y)).sqrt()
                };
                let sin_phi_3 = if p.y == 0.0 {
                    1.0
                } else {
                    p.z / p.y * cos_phi_3
                };

                // construct the x-dependent part
                let rho =
                    (1.0 + self.e_cm_sq.sqrt() / p_abs * (f64::FRAC_PI_2() * mom.t).tan()).ln();
                let xi = f64::PI() * mom.x;
                let ep = rho.sinh() * xi.sin();
                let theta = if mom.y < 0.5 {
                    ((1.0 + ep) * ((1.0 + ep) / ep).powf(-2.0 * mom.y) - ep).acos()
                } else {
                    (ep - (1.0 + ep) * ((1.0 + ep) / ep).powf(-2.0 * (1.0 - mom.y))).acos()
                };
                let phi = f64::PI() * 2.0 * mom.z;

                let k0 = rho.cosh() * xi.cos();
                let k1 = ep * theta.cos();
                let k2 = ep * theta.sin() * phi.cos();
                let k3 = ep * theta.sin() * phi.sin();

                let mut k = LorentzVector::from(
                    k0 * cos_theta_1 - k1 * sin_theta_1,
                    k1 * cos_theta_1 * cos_theta_2 + k0 * cos_theta_2 * sin_theta_1
                        - k2 * sin_theta_2,
                    k2 * cos_theta_2 * cos_phi_3
                        + k1 * cos_theta_1 * cos_phi_3 * sin_theta_2
                        + k0 * cos_phi_3 * sin_theta_1 * sin_theta_2
                        - k3 * sin_phi_3,
                    k3 * cos_phi_3
                        + k2 * cos_theta_2 * sin_phi_3
                        + k1 * cos_theta_1 * sin_theta_2 * sin_phi_3
                        + k0 * sin_theta_1 * sin_theta_2 * sin_phi_3,
                );

                // TODO: add more operators to the vector so that we can remove &s
                k = &self.qs[self.channel] + &(&(&p * 0.5) + &(&k * (0.5 * p_abs)));

                // Compute the Jacobian
                let mut jac = 1.0 / 16.0
                    * p_abs.powi(4)
                    * ep
                    * ep
                    * theta.sin()
                    * (rho.sinh().powi(2) + xi.sin().powi(2));
                jac *= f64::FRAC_PI_2()
                    * (self.e_cm_sq / p_abs / p_abs + (f64::E().powf(rho) - 1.0).powi(2))
                    / (self.e_cm_sq.sqrt() / p_abs * f64::E().powf(rho));
                jac *= f64::PI();
                jac *= 2.0 * (ep + theta.cos().abs()) / theta.sin() * ((1.0 + ep) / ep).ln();
                jac *= 2.0 * f64::PI();

                Ok((k, jac))
            }
        }
    }
}
