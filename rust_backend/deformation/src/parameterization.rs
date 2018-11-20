use num::traits::FloatConst;
use vector::LorentzVector;
use REGION_EXT;

#[derive(Clone)]
pub enum ParameterizationMode {
    Log,
    Linear,
    Weinzierl,
    WeinzierlExternal,
}

#[derive(Clone)]
pub struct Parameterizer {
    mode: ParameterizationMode,
    alpha: f64, // scaling parameter
    e_cm_sq: f64,
    channel: usize,
    qs: Vec<LorentzVector<f64>>,
    region: usize,
}

impl Parameterizer {
    pub fn new(
        e_cm_sq: f64,
        alpha: Option<f64>,
        region: usize,
        channel: usize,
    ) -> Result<Parameterizer, &'static str> {
        Ok(Parameterizer {
            e_cm_sq,
            channel,
            qs: Vec::with_capacity(4),
            alpha: alpha.unwrap_or(e_cm_sq.sqrt()),
            mode: ParameterizationMode::Weinzierl,
            region,
        })
    }

    pub fn set_mode(&mut self, mode: &str) -> Result<(), &'static str> {
        match mode {
            "log" => self.mode = ParameterizationMode::Log,
            "linear" => self.mode = ParameterizationMode::Linear,
            "weinzierl" => self.mode = ParameterizationMode::Weinzierl,
            "weinzierl_ext" => self.mode = ParameterizationMode::WeinzierlExternal,
            _ => return Err("Unknown parameterization mode"),
        }
        Ok(())
    }

    /// Set new qs and update all parameters accordingly.
    pub fn set_qs(&mut self, qs: Vec<LorentzVector<f64>>) {
        self.qs = qs;
    }

    pub fn set_region(&mut self, region: usize) {
        self.region = region;
    }

    pub fn set_channel(&mut self, channel: usize) {
        self.channel = channel;
    }

    fn map_weinzierl_ext(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<f64>, f64), &str> {
        // TODO: mu1 can be anything of the same order as self.e_cm_sq
        let mu1 = self.e_cm_sq.sqrt();
        let k_e = mu1 * (f64::FRAC_PI_2() * mom.t).tan().sqrt();
        let cos_xi = 1.0 - 2.0 * mom.x;
        let sin_xi = (1.0 - cos_xi * cos_xi).sqrt();
        let cos_theta = 1.0 - 2.0 * mom.y;
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let phi = 2.0 * f64::PI() * mom.z;

        let k = LorentzVector::from_args(
            k_e * cos_xi,
            k_e * sin_xi * sin_theta * phi.sin(),
            k_e * sin_xi * sin_theta * phi.cos(),
            k_e * sin_xi * cos_theta,
        );

        let jac = 2.0 * f64::PI() * f64::PI() * k_e * k_e / (mu1 * mu1)
            * (k_e.powi(4) + mu1.powi(4))
            * sin_xi;
        Ok((k, jac))
    }

    fn map_weinzierl_int(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<f64>, f64), &str> {
        // channel starts at 1
        if self.channel == 0 {
            return Err("A channel must be specified for Weinzierl's conformal map");
        }

        // TODO: cache x-independent block
        let p = &self.qs[self.channel] - &self.qs[self.channel - 1];
        let p_abs = p.euclidean_distance();

        let cos_theta_1 = p.t / p_abs;
        let sin_theta_1 = (1.0 - cos_theta_1 * cos_theta_1).sqrt();

        let cos_theta_2 = if p.x == 0.0 {
            0.0
        } else {
            1.0 / (1.0 + (p.y * p.y + p.z * p.z) / (p.x * p.x)).sqrt()
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
        let rho = (1.0 + self.e_cm_sq.sqrt() / p_abs * (f64::FRAC_PI_2() * mom.t).tan()).ln();
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

        let mut k = LorentzVector::from_args(
            k0 * cos_theta_1 - k1 * sin_theta_1,
            k1 * cos_theta_1 * cos_theta_2 + k0 * cos_theta_2 * sin_theta_1 - k2 * sin_theta_2,
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
        k = &self.qs[self.channel - 1] + &(&(&p * 0.5) + &(&k * (0.5 * p_abs)));

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

    /// Map a vector in the unit hypercube to the infinite hypercube.
    /// Also compute the Jacobian.
    pub fn map(&self, mom: &LorentzVector<f64>) -> Result<(LorentzVector<f64>, f64), &str> {
        match self.mode {
            ParameterizationMode::Log => {
                let r = mom.map(|x| self.alpha * (x / (1.0 - x)).ln());

                let mut jacobian = 1.0;
                for i in 0..4 {
                    jacobian *= self.alpha / (mom[i] - mom[i] * mom[i]);
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
                if self.region == REGION_EXT {
                    self.map_weinzierl_ext(mom)
                } else {
                    self.map_weinzierl_int(mom)
                }
            }
            ParameterizationMode::WeinzierlExternal => self.map_weinzierl_ext(mom),
        }
    }
}
