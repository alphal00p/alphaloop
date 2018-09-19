use num;
use num::traits::FloatConst;
use vector::LorentzVector;
type Complex = num::Complex<f64>;

const OFF_SHELL_BOX: usize = 1;
const ON_SHELL_BOX_SUBTRACTED: usize = 2;
const ONE_OFF_SHELL_BOX_SUBTRACTED: usize = 3;
const ON_SHELL_BOX_SUBTRACTED_UV_INT: usize = 4;

pub struct Integrand {
    pub integrand_id: usize,
    pub channel: i32,
    pub region: i32,
    pub mu_sq: f64,
    pub qs: Vec<LorentzVector<f64>>,
}

impl Integrand {
    pub fn new(
        integrand_name: &str,
        channel: i32,
        region: i32,
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
            mu_sq,
            qs: Vec::with_capacity(4),
        })
    }

    pub fn evaluate(&self, mom: &LorentzVector<Complex>) -> Result<Complex, &str> {
        match self.integrand_id {
            OFF_SHELL_BOX => {
                if self.qs.len() != 4 {
                    return Err("Four Qs are required for the box");
                }

                let mut factor = -f64::FRAC_1_PI() * f64::FRAC_1_PI() * Complex::new(0.0, 1.0);

                let mut denominator = Complex::new(0.0, 0.0);
                for q in &self.qs {
                    denominator *= (mom - q).square()
                }

                if self.channel >= 0 && self.channel < 3 {
                    let mut mc_factor = Complex::default();
                    let mut tmp;

                    for (i, q) in self.qs.windows(2).enumerate() {
                        tmp = 1.0 / ((mom - &q[0]).square().norm() * (mom - &q[1]).square().norm());
                        mc_factor *= tmp * tmp;

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
