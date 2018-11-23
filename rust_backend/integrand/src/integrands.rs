use num;
use num::traits::FloatConst;
use vector::LorentzVector;
type Complex = num::Complex<f64>;
use num::traits::Inv;
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
const DIAGONAL_BOX: usize = 11;
const OFF_SHELL_DOUBLE_BOX_SB: usize = 12;

#[inline]
/// Invert with better precision
fn finv(c: Complex) -> Complex {
    let norm = c.norm();
    c.conj() / norm / norm
}

#[derive(Clone)]
pub struct Integrand {
    integrand_id: usize,
    channel: usize,
    region: usize,
    mu_sq: Complex,
    qs: Vec<LorentzVector<f64>>,
    ext: Vec<LorentzVector<f64>>,
    invs: Vec<f64>, // invariants
    shift: LorentzVector<f64>,
    on_shell_flag: usize,
}

impl Integrand {
    pub fn new(
        integrand_name: &str,
        channel: usize,
        region: usize,
        mu_sq: f64,
    ) -> Result<Integrand, &'static str> {
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
            "diagonalbox_direct_integration" => DIAGONAL_BOX,
            "box2L_direct_integration_SB" => OFF_SHELL_DOUBLE_BOX_SB,
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
            on_shell_flag: 0,
            invs: Vec::with_capacity(16),
        })
    }

    pub fn set_channel(&mut self, channel: usize) {
        self.channel = channel;
    }

    pub fn set_region(&mut self, region: usize) {
        self.region = region;
    }

    /// Get an invariant.
    #[inline(always)]
    fn inv(&self, i: usize, j: usize) -> f64 {
        self.invs[i * self.ext.len() + j]
    }
    // build bit-flags for the external momenta
    // this will indicate which counterterms should be used
    pub fn set_on_shell_flag(&mut self, flag: usize) {
        //TODO: print proper warning
        for (i, e) in self.ext.iter().enumerate() {
            if e.square() == 0. && flag & (1 << i) == 0 {
                println!("WARNING: External {} is an unflaged on-shell momentum!", i);
            }
        }
        self.on_shell_flag = flag;
    }

    pub fn set_externals(&mut self, ext: Vec<LorentzVector<f64>>) {
        self.ext = ext;

        if self.on_shell_flag != 0 {
            println!(
                "Using subtraction terms for configuration {}",
                self.on_shell_flag
            );
        }

        // compute the invariants
        self.invs.clear();
        for (i, x) in self.ext.iter().enumerate() {
            for (j, y) in self.ext.iter().enumerate() {
                if i == j {
                    self.invs.push(x.dot(y));
                } else {
                    self.invs.push((x + y).square());
                }
            }
        }

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

    #[inline]
    //Compute the xi's for the collinear limits, here mom should always by on-shell
    pub fn collinear_x(loopmom: &LorentzVector<Complex>, mom: &LorentzVector<f64>) -> f64 {
        //TODO: Check if eta.dot(mom) is not zero
        let eta = mom.dual();
        (eta[0] * loopmom[0].re
            - eta[1] * loopmom[1].re
            - eta[2] * loopmom[2].re
            - eta[3] * loopmom[3].re)
            / eta.dot(mom)
    }

    pub fn evaluate(&self, mom: &[LorentzVector<Complex>]) -> Result<Complex, &'static str> {
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
                if self.integrand_id == ON_SHELL_BOX_SUBTRACTED
                    || self.integrand_id == OFF_SHELL_BOX && self.on_shell_flag == 15
                {
                    let mut ct = Complex::new(1.0, 0.0);
                    for i in 0..4 {
                        ct -= (&mom[0] - &self.qs[i]).square() / self.inv(2, 1 + 2 * (i % 2));
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

                    Ok(factor * finv(denominator) / mc_factor)
                } else {
                    Ok(factor * finv(denominator))
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

                Ok(factor * finv(denominator))
            }
            OFF_SHELL_DOUBLE_TRIANGLE => {
                let mut factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let (k, l) = (&mom[0], &mom[1]);

                let denominator = k.square()
                    * l.square()
                    * (k - l).square()
                    * (&(k - l) + &self.ext[0]).square()
                    * (k + &self.ext[0]).square();

                Ok(factor * finv(denominator))
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

                Ok(factor * finv(denominator))
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

                if self.channel > 0 {
                    let alpha = 3;
                    let mc_factor = l.square().norm().powi(alpha).inv()
                        + (l + &self.ext[0]).square().norm().powi(alpha).inv();

                    if self.channel == 1 {
                        factor *= l.square().norm().powi(alpha).inv();
                    } else {
                        factor *= (l + &self.ext[0]).square().norm().powi(alpha).inv();
                    }

                    Ok(factor * finv(denominator) / mc_factor)
                } else {
                    Ok(factor * finv(denominator))
                }
            }

            DIAGONAL_BOX => {
                let mut factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let (k, l) = (&mom[0], &mom[1]);

                //Propagators
                let d1 = k.square(); //A1
                let d2 = (k - &self.ext[0]).square(); //A2
                let d3 = l.square(); //A3
                let d4 = (l - &self.ext[2]).square(); //A4
                let d5 = ((k - l) - &self.ext[0] - &self.ext[1]).square(); //A5

                let denominator = d1 * d2 * d3 * d4 * d5;

                //println!("Shell Flag: {:#b}", self.on_shell_flag);
                match (
                    self.on_shell_flag & (1 << 0) == 0,
                    self.on_shell_flag & (1 << 1) == 0,
                    self.on_shell_flag & (1 << 2) == 0,
                    self.on_shell_flag & (1 << 3) == 0,
                ) {
                    (false, false, false, false) => {
                        // p1,p2,p3,p4 : on-shell and p1,p3 : off-shell
                        Ok(factor * finv(denominator))
                    }
                    (true, _, true, _) => {
                        // p2,p4 : on-shell and p1,p3 : off-shell
                        //Collinear Limits
                        let x1 = Integrand::collinear_x(k, &self.ext[0]);
                        let x3 = Integrand::collinear_x(l, &self.ext[2]);

                        let k_c1 = self.ext[0] * x1;
                        let l_c3 = self.ext[0] * x3;

                        let c1d5 = (l - &k_c1 + &self.ext[0] + &self.ext[1]).square(); //A5
                        let c3d5 = (k - &l_c3 - &self.ext[0] - &self.ext[1]).square(); //A5
                        let c13d5 = (k_c1 - l_c3 - &self.ext[0] - &self.ext[1]).square(); //A5

                        //Compute collinear counterterms
                        //k -> x1 p1
                        let c1 = (finv(d1 * d2) - finv((d1 - self.mu_sq) * (d2 - self.mu_sq)))
                            * finv(d3 * d4 * c1d5);

                        //l -> x3 p3
                        let c3 = (finv(d3 * d4) - finv((d3 - self.mu_sq) * (d4 - self.mu_sq)))
                            * finv(d1 * d2 * c3d5);

                        //k -> x1 p1 and l -> x3 p3
                        let c13 = (finv(d1 * d2 * d3 * d4)
                            - finv(
                                (d1 - self.mu_sq)
                                    * (d2 - self.mu_sq)
                                    * (d3 - self.mu_sq)
                                    * (d4 - self.mu_sq),
                            ))
                            / c13d5;

                        Ok(factor * (finv(denominator) - c1 - c3 + c13))
                    }
                    _ => Err("Unknown Subtracted Term for SUBTRACTED_DIAGONAL_BOX!"),
                }
            }

            OFF_SHELL_DOUBLE_BOX_SB => {
                // let factor = Complex::new(-f64::FRAC_1_PI().powi(4), 0.0);
                let factor = Complex::new(1.0, 0.0);
                let (k, l) = (&mom[0], &mom[1]);
                let p = (&self.ext[0], &self.ext[1], &self.ext[2], &self.ext[3]);
                let denominator = (k + p.2).square()
                    * (k + &(p.1 + p.2)).square()
                    * (k - p.3).square()
                    * (l - p.3).square()
                    * (l).square()
                    * (l + p.2).square()
                    * (k - l).square();
                Ok(factor * finv(denominator))
            }

            _ => Err("Integrand is not implemented yet"),
        }
    }
}
