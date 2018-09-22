use num;
use std::mem;
use vector::LorentzVector;
use REGION_EXT;

type Complex = num::Complex<f64>;

pub struct Deformer {
    qs: Vec<LorentzVector<f64>>,
    masses: Vec<f64>,
    p_plus: LorentzVector<f64>,
    p_min: LorentzVector<f64>,
    e_cm_sq: f64,
    mu_p: f64,
    m1_fac: f64,
    m2_fac: f64,
    m3_fac: f64,
    m4_fac: f64,
    m1_sq: f64,
    m2_sq: f64,
    m3_sq: f64,
    m4_sq: f64,
    gamma1: f64,
    gamma2: f64,
    e_soft: f64,
    region: usize,
    mu_sq: Complex,
    uv_shift: LorentzVector<f64>,
}

impl Deformer {
    pub fn new(e_cm_sq: f64, mu_sq: f64, region: usize) -> Result<Deformer, &'static str> {
        Ok(Deformer {
            qs: Vec::with_capacity(4),
            masses: Vec::with_capacity(4),
            p_plus: LorentzVector::new(),
            p_min: LorentzVector::new(),
            e_cm_sq,
            mu_p: 0.0,
            m1_fac: 0.035,
            m2_fac: 0.7,
            m3_fac: 0.035,
            m4_fac: 0.035,
            m1_sq: 0.0,
            m2_sq: 0.0,
            m3_sq: 0.0,
            m4_sq: 0.0,
            gamma1: 0.7,
            gamma2: 0.008,
            e_soft: 0.03,
            mu_sq: Complex::new(0.0, mu_sq),
            uv_shift: LorentzVector::new(),
            region,
        })
    }

    /// Set new qs and update all parameters accordingly.
    pub fn set_qs(&mut self, qs: Vec<LorentzVector<f64>>) {
        self.qs = qs;

        if self.masses.len() != self.qs.len() {
            self.masses = vec![0.0; self.qs.len()];
        }

        self.uv_shift = LorentzVector::new();
        for q in &self.qs {
            self.uv_shift = &self.uv_shift + &q;
        }
        self.uv_shift = &self.uv_shift * (1.0 / (self.qs.len() as f64));

        self.p_plus = self.compute_p(true);
        self.p_min = self.compute_p(false);
        self.mu_p = (&self.p_min - &self.p_plus).square();

        self.m1_sq = self.m1_fac * self.m1_fac * self.e_cm_sq;
        self.m2_sq = self.m2_fac
            * self.m2_fac
            * (if self.mu_p > self.e_cm_sq {
                self.mu_p
            } else {
                self.e_cm_sq
            });
        self.m3_sq = self.m3_fac
            * self.m3_fac
            * (if self.mu_p > self.e_cm_sq {
                self.mu_p
            } else {
                self.e_cm_sq
            });
        self.m4_sq = self.m4_fac * self.m4_fac * self.e_cm_sq;
    }

    ///Find a vector P such that all external momenta `qs` are in
    ///the forward lightcone of P if `plus`, or are in the backwards lightcone
    ///if not `plus`.
    fn compute_p(&self, plus: bool) -> LorentzVector<f64> {
        /// Interpolate between vectors
        fn z(a: LorentzVector<f64>, b: LorentzVector<f64>, plus: bool) -> LorentzVector<f64> {
            // WARNING: watch out for configurations with input momenta that are back-to-back as then the spatial distance is 0
            let n = 1.0 / b.spatial_distance();

            if plus {
                &LorentzVector::from(
                    a.t + n * b.square() - n * b.t * b.t,
                    a.x - n * b.t * b.x,
                    a.y - n * b.t * b.y,
                    a.z - n * b.t * b.z,
                ) * 0.5
            } else {
                &LorentzVector::from(
                    a.t - n * b.square() + n * b.t * b.t,
                    a.x + n * b.t * b.x,
                    a.y + n * b.t * b.y,
                    a.z + n * b.t * b.z,
                ) * 0.5
            }
        }

        let mut orig_vec = self.qs.clone();
        let mut filtered_vec = vec![];
        loop {
            // filter all vectors that are in the forward/backward light-cone of another vector
            for (i, v) in orig_vec.iter().enumerate() {
                let mut spacelike = true;
                for (j, v1) in orig_vec.iter().enumerate() {
                    if i == j {
                        spacelike = false;
                        continue;
                    }
                    if (v - v1).square() > 0.0 && ((plus && v.t > v1.t) || (!plus && v.t < v1.t)) {
                        spacelike = false;
                        break;
                    }
                }
                if spacelike {
                    // this vector is space-like relative to all others
                    filtered_vec.push(v.clone())
                }
            }

            if filtered_vec.len() == 0 {
                panic!("Failed to determine P_{+/-}.");
            }
            if filtered_vec.len() == 1 {
                break;
            }

            // find the pair with the smallest space-like seperation
            let (mut first, mut sec, mut min_sep) = (0, 0, 0.0);
            for (i, v1) in filtered_vec.iter().enumerate() {
                for (j, v2) in filtered_vec.iter().enumerate() {
                    let sep = -(v1 - v2).square();
                    if min_sep == 0.0 || sep <= min_sep {
                        min_sep = sep;
                        first = i;
                        sec = j;
                    }
                }
            }

            // replace first vector and drop the second
            filtered_vec[first] = z(
                &filtered_vec[first] + &filtered_vec[sec],
                &filtered_vec[first] - &filtered_vec[sec],
                plus,
            );
            filtered_vec.swap_remove(sec);
            mem::swap(&mut orig_vec, &mut filtered_vec);
            if orig_vec.len() == 1 {
                break;
            }
        }

        assert!(orig_vec.len() == 1);
        orig_vec.pop().unwrap()
    }

    /// The three helper functions h_delta-, h_delta+, and h_delta0, indicated by `sign`.
    fn h_delta(sign: i32, k: &LorentzVector<f64>, mass: f64, m: f64) -> f64 {
        let v = if sign == 0 {
            (k.t.abs() - (k.spatial_squared() + mass).sqrt()).powi(2)
        } else {
            (sign as f64 * k.t.abs() - (k.spatial_squared() + mass).sqrt()).powi(2)
        };
        v / (v + m)
    }

    fn h_theta(t: f64, m: f64) -> f64 {
        if t <= 0.0 {
            0.0
        } else {
            t / (t + m)
        }
    }

    fn g(k: &LorentzVector<f64>, gamma: f64, m: f64) -> f64 {
        gamma * m / (k.t * k.t + k.spatial_squared() + m)
    }

    fn d(&self, mom: &LorentzVector<f64>, i: usize, l: usize) -> f64 {
        if self.masses[l] == 0.0 {
            if l == i {
                return 1.0;
            }

            if (&self.qs[i] - &self.qs[l]).square() == 0.0 {
                if self.qs[i].t < self.qs[l].t {
                    return Deformer::h_delta(1, &(mom - &self.qs[l]), 0.0, self.m1_sq);
                } else {
                    return Deformer::h_delta(-1, &(mom - &self.qs[l]), 0.0, self.m1_sq);
                }
            }
        }

        let hd = Deformer::h_delta(
            0,
            &(mom - &self.qs[l]),
            self.masses[l] * self.masses[l],
            self.m1_sq,
        );
        let ht = Deformer::h_theta(
            -2.0 * (mom - &self.qs[l]).dot(&(mom - &self.qs[i])),
            self.m1_sq,
        );
        if hd > ht {
            hd
        } else {
            ht
        }
    }

    fn deform_int(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, Complex), &str> {
        let (mut c_plus, mut c_minus) = (1.0, 1.0);

        for (qi, mi) in self.qs.iter().zip(&self.masses) {
            // note the sign reversal
            c_plus *= Deformer::h_delta(-1, &(mom - qi), mi * mi, self.m3_sq);
            c_minus *= Deformer::h_delta(1, &(mom - qi), mi * mi, self.m3_sq);
        }

        let k_plus = mom - &self.p_plus;
        let k_minus = mom - &self.p_min;
        let k_ext = LorentzVector::from(
            c_plus * k_plus.t + c_minus * k_minus.t,
            -c_plus * k_plus.x - c_minus * k_minus.x,
            -c_plus * k_plus.y - c_minus * k_minus.y,
            -c_plus * k_plus.z - c_minus * k_minus.z,
        );
        let k_centre = &(&k_plus + &k_minus) * 0.5;

        // TODO: cache qi= mom-qs[i]?

        let n = self.qs.len();
        let f = Deformer::g(&k_centre, self.gamma1, self.m2_sq);
        let mut c = vec![f; n];

        let mut k_int = LorentzVector::new();
        for i in 0..n {
            for l in 0..n {
                // note: in paper l starts at 1
                c[i] *= Deformer::d(self, mom, i, l);
            }

            k_int = &k_int + &(&(mom - &self.qs[i]) * -c[i]); // massless part
        }

        let k0 = &k_int + &k_ext;
        let mut lambda_cycle = 1.0; // maximum lambda value

        for j in 0..n {
            let xj = (k0.dot(&(mom - &self.qs[j])) / k0.square()).powi(2);
            let yj = ((mom - &self.qs[j]).square() - self.masses[j].powi(2)) / k0.square();

            if 2.0 * xj < yj {
                if (yj / 4.0).sqrt() < lambda_cycle {
                    lambda_cycle = (yj / 4.0).sqrt();
                }
            } else if yj < 0.0 {
                if (xj - yj / 2.0).sqrt() < lambda_cycle {
                    lambda_cycle = (xj - yj / 2.0).sqrt();
                }
            } else {
                if (xj - yj / 4.0).sqrt() < lambda_cycle {
                    lambda_cycle = (xj - yj / 4.0).sqrt();
                }
            };
        }

        let n1: f64 = c.iter().sum();
        lambda_cycle = 1.0 / (4.0 * n1);
        if 1.0 / (4.0 * n1) < lambda_cycle {}

        let uv_fac = 4.0 * (mom - &self.uv_shift).dot(&k0);
        if uv_fac > self.mu_sq.im {
            if lambda_cycle > self.mu_sq.im / uv_fac {
                lambda_cycle = self.mu_sq.im / uv_fac;
            }
        }

        let k = &mom.to_complex(true) + &(&(&k_int + &k_ext) * lambda_cycle).to_complex(false);

        // TODO: implement jacobian analytically or numerically with multi-loop support
        let jac = Complex::new(1.0, 0.0);
        Ok((k, jac))
    }

    fn deform_ext(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, Complex), &str> {
        let k = LorentzVector::from(
            Complex::new(mom.t, mom.t - self.uv_shift.t),
            Complex::new(mom.x, -mom.x - self.uv_shift.x),
            Complex::new(mom.y, -mom.y - self.uv_shift.y),
            Complex::new(mom.z, -mom.z - self.uv_shift.z),
        );

        Ok((k, -Complex::new(0.0, -4.0)))
    }

    pub fn deform(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, Complex), &str> {
        if self.region == REGION_EXT {
            self.deform_ext(mom)
        } else {
            self.deform_int(mom)
        }
    }
}
