use num;
use num::traits::FloatConst;
use vector::LorentzVector;
type Complex = num::Complex<f64>;

const REGION_ALL: usize = 0;
const REGION_EXT: usize = 1;
const REGION_INT: usize = 2;

const OFF_SHELL_BOX: usize = 1;
const ON_SHELL_BOX_SUBTRACTED: usize = 2;
const ONE_OFF_SHELL_BOX_SUBTRACTED: usize = 3;
const ON_SHELL_BOX_SUBTRACTED_UV_INT: usize = 4;

pub struct Integrand {
    pub integrand_id: usize,
    pub channel: usize,
    pub region: usize,
    pub mu_sq: Complex,
    pub qs: Vec<LorentzVector<f64>>,
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
            "box1L_direct_integration_subtracted" => ON_SHELL_BOX_SUBTRACTED,
            "box1L_direct_integration_one_offshell_subtracted" => ONE_OFF_SHELL_BOX_SUBTRACTED,
            "box1L_direct_integration_subtracted_uv_int" => ON_SHELL_BOX_SUBTRACTED_UV_INT,
            _ => return Err("Unknown integrand"),
        };

        Ok(Integrand {
            integrand_id,
            channel,
            region,
            mu_sq: Complex::new(0.0, mu_sq),
            qs: Vec::with_capacity(4),
        })
    }

    pub fn evaluate(&self, mom: &LorentzVector<Complex>) -> Result<Complex, &str> {
        match self.integrand_id {
            OFF_SHELL_BOX | ON_SHELL_BOX_SUBTRACTED => {
                if self.qs.len() != 4 {
                    return Err("Four Qs are required for the box");
                }

                let mut factor = Complex::new(0.0, -f64::FRAC_1_PI() * f64::FRAC_1_PI());

                let mut denominator = Complex::new(1.0, 0.0);

                if self.region == REGION_ALL || self.region == REGION_INT {
                    for q in &self.qs {
                        denominator *= (mom - q).square();
                    }
                }

                if self.region != REGION_ALL {
                    // TODO: cache
                    let mut shift = LorentzVector::new();
                    for q in &self.qs {
                        shift = &shift + q;
                    }
                    shift = &shift * 0.25;
                    let mut d = (mom - &shift).square() - self.mu_sq;
                    d = d.powf(4.0);

                    if self.region == REGION_EXT {
                        denominator = d;
                    } else {
                        factor *= 1.0 - denominator / d;
                    }
                }

                // add the subtraction terms
                if self.integrand_id == ON_SHELL_BOX_SUBTRACTED {
                    let mut ct = Complex::new(1.0, 0.0);
                    let mut invariant;
                    for i in 0..4 {
                        // TODO: cache
                        invariant = (&self.qs[i % 2] - &self.qs[i % 2 + 2]).square();
                        ct -= (mom - &self.qs[i]).square() / invariant;
                    }
                    factor *= ct;
                }

                if self.channel > 0 && self.channel <= 3 {
                    let mut mc_factor = Complex::default();
                    let mut tmp;

                    for (i, q) in self.qs.windows(2).enumerate() {
                        tmp = 1.0 / ((mom - &q[0]).square().norm() * (mom - &q[1]).square().norm());
                        mc_factor += tmp * tmp;

                        if i == self.channel as usize {
                            factor *= tmp * tmp;
                        }
                    }

                    Ok(factor / denominator / mc_factor)
                } else {
                    Ok(factor / denominator)
                }
            }
            _ => Err("Integrand is not implemented yet"),
        }
    }
}
