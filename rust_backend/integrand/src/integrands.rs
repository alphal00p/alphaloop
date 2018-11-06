use num;
use num::traits::FloatConst;
use vector::LorentzVector;
type Complex = num::Complex<f64>;
use std::f64;

const REGION_ALL: usize = 0;
const REGION_EXT: usize = 1;
const REGION_INT: usize = 2;

const OFF_SHELL_BOX: usize = 1;
const ON_SHELL_BOX_SUBTRACTED: usize = 2;
const ONE_OFF_SHELL_BOX_SUBTRACTED: usize = 3;
const ON_SHELL_BOX_SUBTRACTED_UV_INT: usize = 4;
const OFF_SHELL_PENTAGON: usize = 5;
const OFF_SHELL_HEXAGON: usize = 6;
const OFF_SHELL_DOUBLE_BOX: usize = 7;
const OFF_SHELL_DOUBLE_TRIANGLE: usize = 8;
const OFF_SHELL_TRIANGLE_BOX: usize = 9;
const OFF_SHELL_TRIANGLE_BOX_ALTERNATIVE: usize = 10;

#[derive(Clone)]
pub struct Integrand {
    integrand_id: usize,
    channel: usize,
    region: usize,
    mu_sq: Complex,
    qs: Vec<LorentzVector<f64>>,
    ext: Vec<LorentzVector<f64>>,
    shift: LorentzVector<f64>,
}

impl Integrand {
    pub fn new(
        integrand_name: &str,
        channel: usize,
        region: usize,
        mu_sq: f64,
    ) -> Result<Integrand, &str> {
        let integrand_id = match integrand_name {
            "box1L_direct_integration" => OFF_SHELL_BOX,
            "pentagon1L_direct_integration" => OFF_SHELL_PENTAGON,
            "hexagon1L_direct_integration" => OFF_SHELL_HEXAGON,
            "box1L_direct_integration_subtracted" => ON_SHELL_BOX_SUBTRACTED,
            "box1L_direct_integration_one_offshell_subtracted" => ONE_OFF_SHELL_BOX_SUBTRACTED,
            "box1L_direct_integration_subtracted_uv_int" => ON_SHELL_BOX_SUBTRACTED_UV_INT,
            "box2L_direct_integration" => OFF_SHELL_DOUBLE_BOX,
            "triangle2L_direct_integration" => OFF_SHELL_DOUBLE_TRIANGLE,
            "trianglebox_direct_integration" => OFF_SHELL_TRIANGLE_BOX,
            "trianglebox_alternative_direct_integration" => OFF_SHELL_TRIANGLE_BOX_ALTERNATIVE,
            _ => return Err("Unknown integrand"),
        };

        Ok(Integrand {
            integrand_id,
            channel,
            region,
            mu_sq: Complex::new(0.0, mu_sq),
            qs: Vec::with_capacity(4),
            ext: Vec::with_capacity(4),
            shift: LorentzVector::new(),
        })
    }

    pub fn set_channel(&mut self, channel: usize) {
        self.channel = channel;
    }

    pub fn set_region(&mut self, region: usize) {
        self.region = region;
    }

    pub fn set_externals(&mut self, ext: Vec<LorentzVector<f64>>) {
        self.ext = ext;

        // compute the qs for one-loop computations
        self.qs.clear();
        self.qs.push(LorentzVector::new());
        for (i, e) in self.ext[1..].iter().enumerate() {
            let r = &self.qs[i] + e;
            self.qs.push(r);
        }

        // compute the shift
        self.shift = LorentzVector::new();
        for q in &self.qs {
            self.shift = &self.shift + q;
        }
        self.shift = &self.shift * 0.25;
    }

    pub fn evaluate(&self, mom: &[LorentzVector<Complex>]) -> Result<Complex, &str> {
        match self.integrand_id {
            OFF_SHELL_BOX | ON_SHELL_BOX_SUBTRACTED | OFF_SHELL_PENTAGON | OFF_SHELL_HEXAGON => {
                let mut factor = Complex::new(0.0, -f64::FRAC_1_PI() * f64::FRAC_1_PI());

                let mut denominator = Complex::new(1.0, 0.0);

                if self.region == REGION_ALL || self.region == REGION_INT {
                    for q in &self.qs {
                        denominator *= (&mom[0] - q).square();
                    }
                }

                if self.region != REGION_ALL {
                    let mut d = (&mom[0] - &self.shift).square() - self.mu_sq;
                    d = d.powf(self.qs.len() as f64); // TODO: use powi when it is implemented

                    if self.region == REGION_EXT {
                        denominator = d;
                    } else {
                        if d.is_finite() {
                            // the d could overflow
                            factor *= 1.0 - denominator / d;
                        }
                    }
                }

                // add the subtraction terms
                if self.integrand_id == ON_SHELL_BOX_SUBTRACTED {
                    let mut ct = Complex::new(1.0, 0.0);
                    let mut invariant;
                    for i in 0..4 {
                        // TODO: cache
                        invariant = (&self.qs[i % 2] - &self.qs[i % 2 + 2]).square();
                        ct -= (&mom[0] - &self.qs[i]).square() / invariant;
                    }
                    factor *= ct;
                }

                if self.channel > 0 && self.channel <= self.qs.len() - 1 {
                    let mut mc_factor = Complex::default();
                    let mut tmp;

                    for (i, q) in self.qs.windows(2).enumerate() {
                        tmp = 1.0
                            / ((&mom[0] - &q[0]).square().norm()
                                * (&mom[0] - &q[1]).square().norm());
                        mc_factor += tmp * tmp;

                        if i + 1 == self.channel as usize {
                            factor *= tmp * tmp;
                        }
                    }

                    Ok(factor / denominator / mc_factor)
                } else {
                    Ok(factor / denominator)
                }
            }
            OFF_SHELL_DOUBLE_BOX => {
                let mut factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let (k, l) = (&mom[0], &mom[1]);

                let denominator = k.square()
                    * l.square()
                    * (k - l - &self.ext[0]).square()
                    * (l - &self.ext[1]).square()
                    * (l - &self.ext[1] - &self.ext[2]).square()
                    * (k - &self.ext[0]).square()
                    * (k + &self.ext[3]).square();

                Ok(factor / denominator)
            }
            OFF_SHELL_DOUBLE_TRIANGLE => {
                let mut factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let (k, l) = (&mom[0], &mom[1]);

                let denominator = k.square()
                    * l.square()
                    * (k - l).square()
                    * (&(k - l) + &self.ext[0]).square()
                    * (k + &self.ext[0]).square();

                Ok(factor / denominator)
            }
            OFF_SHELL_TRIANGLE_BOX => {
                let mut factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let (k, l) = (&mom[0], &mom[1]);

                let denominator = k.square()
                    * l.square()
                    * (k - l).square()
                    * (&(k - l) - &self.ext[1]).square()
                    * (&(k - l) + &self.ext[0]).square()
                    * (k + &self.ext[0]).square();

                Ok(factor / denominator)
            }
            OFF_SHELL_TRIANGLE_BOX_ALTERNATIVE => {
                let mut factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let (k, l) = (&mom[0], &mom[1]);

                let denominator = k.square()
                    * l.square()
                    * (k - l - &self.ext[0]).square()
                    * (l - &self.ext[1]).square()
                    * (l + &self.ext[0]).square()
                    * (k - &self.ext[0]).square();

                Ok(factor / denominator)
            }
            _ => Err("Integrand is not implemented yet"),
        }
    }
}
