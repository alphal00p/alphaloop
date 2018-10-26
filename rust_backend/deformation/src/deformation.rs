use std::f64::EPSILON;
use utils::{determinant, determinant8};
use vector::LorentzVector;
use Complex;
use REGION_EXT;

pub const DOUBLE_BOX_ID: usize = 0;
pub const DOUBLE_TRIANGLE_ID: usize = 1;

const DIRECTIONS: [LorentzVector<f64>; 4] = [
    LorentzVector {
        t: 1.0,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
    LorentzVector {
        t: 0.0,
        x: 1.0,
        y: 0.0,
        z: 0.0,
    },
    LorentzVector {
        t: 0.0,
        x: 0.0,
        y: 1.0,
        z: 0.0,
    },
    LorentzVector {
        t: 0.0,
        x: 0.0,
        y: 0.0,
        z: 1.0,
    },
];

pub struct Deformer {
    qs: Vec<LorentzVector<f64>>,
    ext: Vec<LorentzVector<f64>>, // external momenta
    masses: Vec<f64>,
    p_buf: Vec<LorentzVector<f64>>, // buffer for p computation
    p_plus: LorentzVector<f64>,
    p_min: LorentzVector<f64>,
    e_cm_sq: f64,
    mu_p_sq: f64,
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
    pub fn new(
        e_cm_sq: f64,
        mu_sq: f64,
        region: usize,
        masses: Vec<f64>,
    ) -> Result<Deformer, &'static str> {
        Ok(Deformer {
            qs: Vec::with_capacity(6),
            ext: Vec::with_capacity(4),
            masses,
            p_buf: Vec::with_capacity(6),
            p_plus: LorentzVector::new(),
            p_min: LorentzVector::new(),
            e_cm_sq,
            mu_p_sq: 0.0,
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

    /// Set external momenta. Only used for the double box at the moment.
    pub fn set_external_momenta(&mut self, ext: Vec<LorentzVector<f64>>) {
        self.ext = ext;
    }

    pub fn set_region(&mut self, region: usize) {
        self.region = region;
    }

    /// Set new qs and update all parameters accordingly.
    pub fn set_qs(&mut self, qs: &[LorentzVector<f64>]) {
        self.qs.clear();
        self.qs.extend_from_slice(qs);

        if self.masses.len() != self.qs.len() {
            self.masses.resize(self.qs.len(), 0.);
        }

        self.uv_shift = LorentzVector::new();
        for q in &self.qs {
            self.uv_shift = &self.uv_shift + q;
        }
        self.uv_shift = &self.uv_shift * (1.0 / (self.qs.len() as f64));

        self.p_plus = self.compute_p(true);
        self.p_min = self.compute_p(false);
        self.shift_and_check_p();

        self.mu_p_sq = (&self.p_min - &self.p_plus).square();

        self.m1_sq = self.m1_fac * self.m1_fac * self.e_cm_sq;
        self.m2_sq = self.m2_fac
            * self.m2_fac
            * (if self.mu_p_sq > self.e_cm_sq {
                self.mu_p_sq
            } else {
                self.e_cm_sq
            });
        self.m3_sq = self.m3_fac
            * self.m3_fac
            * (if self.mu_p_sq > self.e_cm_sq {
                self.mu_p_sq
            } else {
                self.e_cm_sq
            });
        self.m4_sq = self.m4_fac * self.m4_fac * self.e_cm_sq;
    }

    ///Find a vector P such that all external momenta `qs` are in
    ///the forward lightcone of P if `plus`, or are in the backwards lightcone
    ///if not `plus`.
    fn compute_p(&mut self, plus: bool) -> LorentzVector<f64> {
        /// Interpolate between vectors
        #[inline]
        fn z(a: LorentzVector<f64>, b: LorentzVector<f64>, plus: bool) -> LorentzVector<f64> {
            if b.euclidean_square() < 1e-9 {
                return a;
            }

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

        self.p_buf.clear();
        self.p_buf.extend_from_slice(&self.qs);
        let mut keep_map = [true; 10]; // we will not have more than 10 propagators
        assert!(self.p_buf.len() < 10);

        loop {
            unsafe {
                // filter all vectors that are in the forward/backward light-cone of another vector
                for (i, v) in self.p_buf.iter().enumerate() {
                    if *keep_map.get_unchecked(i) {
                        for (j, v1) in self.p_buf[i + 1..].iter().enumerate() {
                            if *keep_map.get_unchecked(i + 1 + j) && (v - v1).square() >= 0.0 {
                                *keep_map.get_unchecked_mut(i) &=
                                    (plus && v.t <= v1.t) || (!plus && v.t >= v1.t);
                                *keep_map.get_unchecked_mut(i + 1 + j) &=
                                    (plus && v.t >= v1.t) || (!plus && v.t <= v1.t);

                                if !*keep_map.get_unchecked(i) {
                                    break;
                                }
                            }
                        }
                    }
                }

                let mut i = 0;
                while i < self.p_buf.len() {
                    if !keep_map.get_unchecked(i) {
                        self.p_buf.swap_remove(i);
                        *keep_map.get_unchecked_mut(i) = *keep_map.get_unchecked(self.p_buf.len());
                    } else {
                        i += 1;
                    }
                }
            }

            if self.p_buf.len() == 0 {
                panic!("Failed to determine P_{+/-}.");
            }
            if self.p_buf.len() == 1 {
                break;
            }

            // find the pair with the smallest space-like separation
            let (mut first, mut sec, mut min_sep) = (0, 1, 0.0);
            for (i, v1) in self.p_buf.iter().enumerate() {
                for (j, v2) in self.p_buf[i + 1..].iter().enumerate() {
                    let sep = -(v1 - v2).square();
                    if min_sep == 0.0 || sep <= min_sep {
                        min_sep = sep;
                        first = i;
                        sec = i + j + 1;
                    }
                }
            }

            // replace first vector and drop the second
            self.p_buf[first] = z(
                &self.p_buf[first] + &self.p_buf[sec],
                &self.p_buf[first] - &self.p_buf[sec],
                plus,
            );
            self.p_buf.swap_remove(sec);
            if self.p_buf.len() == 1 {
                break;
            }
        }

        assert!(self.p_buf.len() == 1);
        self.p_buf.pop().unwrap()
    }

    /// Shift P+ and P- outward and check if their new values fulfills P^2 >= 0.0 and (+/-) * P+/-.t <= 0.0  
    /// This is necessary to counterbalance numerical error arising in the exact algorithm
    fn shift_and_check_p(&mut self) {
        //Perform the shift
        let r = (&self.p_min - &self.p_plus) * 0.5;
        let center = (&self.p_min + &self.p_plus) * 0.5;
        let shift_size = 1.0e-8;

        self.p_min = &center + &r * (1.0 + shift_size);
        self.p_plus = &center - &r * (1.0 + shift_size);

        //Check if P+/- has all the qs in its backward/forward light-cone
        let mut pq_diff;
        for q in &self.qs {
            pq_diff = &self.p_plus - q;
            if pq_diff.t > 0.0 || pq_diff.square() < 0.0 {
                panic!(
                    "P_plus is not correctly defined! P.t = {0}, P^2 = {1}, qs = {2:?}",
                    pq_diff.t,
                    pq_diff.square(),
                    self.qs
                );
            }

            pq_diff = &self.p_min - q;
            if pq_diff.t < 0.0 || pq_diff.square() < 0.0 {
                panic!(
                    "P_minus is not correctly defined! P.t = {0}, P^2 = {1}, qs = {2:?}",
                    pq_diff.t,
                    pq_diff.square(),
                    self.qs
                );
            }
        }
    }

    /// The three helper functions h_delta-, h_delta+, and h_delta0, indicated by `sign`.
    fn h_delta(sign: i32, k: &LorentzVector<f64>, mass: f64, m: f64) -> f64 {
        let v = if sign == 0 {
            (k.t.abs() - (k.spatial_squared() + mass).sqrt()).powi(2)
        } else {
            (sign as f64 * k.t - (k.spatial_squared() + mass).sqrt()).powi(2)
        };
        v / (v + m)
    }

    /// Derivative of h_delta wrt ln k
    fn h_delta_ln_grad(sign: i32, k: &LorentzVector<f64>, mass: f64, m: f64) -> LorentzVector<f64> {
        let w = (k.spatial_squared() + mass).sqrt();
        let mut w_grad = k * (-1.0 / w);

        let w1 = if sign == 0 {
            if k.t > 0.0 {
                w_grad.t = 1.0;
            } else {
                w_grad.t = -1.0;
            }

            k.t.abs() - w
        } else {
            w_grad.t = sign as f64;
            sign as f64 * k.t - w
        };

        &w_grad * (2.0 * m / w1 / (w1 * w1 + m))
    }

    #[inline]
    fn h_theta(t: f64, m: f64) -> f64 {
        if t <= 0.0 {
            0.0
        } else {
            t / (t + m)
        }
    }

    /// Derivative of h_theta wrt ln t
    #[inline]
    fn h_theta_ln_grad(t: f64, m: f64) -> f64 {
        if t <= 0.0 {
            0.0
        } else {
            m / t / (t + m)
        }
    }

    #[inline]
    fn g(k: &LorentzVector<f64>, gamma: f64, m: f64) -> f64 {
        gamma * m / (k.t * k.t + k.spatial_squared() + m)
    }

    /// Derivative of h_theta wrt ln k
    #[inline]
    fn g_ln_grad(k: &LorentzVector<f64>, m: f64) -> LorentzVector<f64> {
        k * (-2.0 / (k.euclidean_square() + m))
    }

    /// Return d and the log partial derivative
    fn d_and_grad(
        &self,
        mom: &LorentzVector<f64>,
        i: usize,
        l: usize,
    ) -> (f64, LorentzVector<f64>) {
        if self.masses[l] == 0.0 {
            if l == i {
                return (1.0, LorentzVector::new());
            }

            if (&self.qs[i] - &self.qs[l]).square() == 0.0 {
                if self.qs[i].t < self.qs[l].t {
                    return (
                        Deformer::h_delta(1, &(mom - &self.qs[l]), 0.0, self.m1_sq),
                        Deformer::h_delta_ln_grad(1, &(mom - &self.qs[l]), 0.0, self.m1_sq),
                    );
                } else {
                    return (
                        Deformer::h_delta(-1, &(mom - &self.qs[l]), 0.0, self.m1_sq),
                        Deformer::h_delta_ln_grad(-1, &(mom - &self.qs[l]), 0.0, self.m1_sq),
                    );
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
            (
                hd,
                Deformer::h_delta_ln_grad(
                    0,
                    &(mom - &self.qs[l]),
                    self.masses[l] * self.masses[l],
                    self.m1_sq,
                ),
            )
        } else {
            (
                ht,
                &(&(mom - &self.qs[l]).dual() + &(mom - &self.qs[i]).dual())
                    * (-2.0 * Deformer::h_theta_ln_grad(
                        -2.0 * (mom - &self.qs[l]).dot(&(mom - &self.qs[i])),
                        self.m1_sq,
                    )),
            )
        }
    }

    /// Return d
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
        let (mut c_plus_grad, mut c_minus_grad) = (LorentzVector::new(), LorentzVector::new());

        for (qi, mi) in self.qs.iter().zip(&self.masses) {
            // note the sign reversal
            c_plus *= Deformer::h_delta(-1, &(mom - qi), mi * mi, self.m3_sq);
            c_minus *= Deformer::h_delta(1, &(mom - qi), mi * mi, self.m3_sq);

            c_plus_grad =
                &c_plus_grad + &Deformer::h_delta_ln_grad(-1, &(mom - qi), mi * mi, self.m3_sq);
            c_minus_grad =
                &c_minus_grad + &Deformer::h_delta_ln_grad(1, &(mom - qi), mi * mi, self.m3_sq);
        }

        c_plus_grad = &c_plus_grad * c_plus;
        c_minus_grad = &c_minus_grad * c_minus;

        let k_plus = mom - &self.p_plus;
        let k_minus = mom - &self.p_min;
        let k_ext = LorentzVector::from(
            c_plus * k_plus.t + c_minus * k_minus.t,
            -c_plus * k_plus.x - c_minus * k_minus.x,
            -c_plus * k_plus.y - c_minus * k_minus.y,
            -c_plus * k_plus.z - c_minus * k_minus.z,
        );

        let mut k_grad = [[0.0; 4]; 4];
        let metric_diag = [1.0, -1.0, -1.0, -1.0];
        for i in 0..4 {
            k_grad[i][i] = metric_diag[i] * (c_plus + c_minus);
            for j in 0..4 {
                k_grad[i][j] += metric_diag[i] * c_plus_grad[j] * k_plus[i]
                    + metric_diag[i] * c_minus_grad[j] * k_minus[i];
            }
        }

        // internal part
        let k_centre = &(&k_plus + &k_minus) * 0.5;

        // TODO: cache qi= mom-qs[i]?

        let n = self.qs.len();
        assert!(n <= 10);
        let f = Deformer::g(&k_centre, self.gamma1, self.m2_sq);
        let mut c = [f; 10];
        let f_grad = Deformer::g_ln_grad(&k_centre, self.m2_sq);
        let mut c_grad : [LorentzVector<f64>; 10] = Default::default();

        let mut k_int = LorentzVector::new();
        for i in 0..n {
            for l in 0..n {
                // note: in paper l starts at 1
                let (c_c, c_grad_c) = Deformer::d_and_grad(self, mom, i, l);
                c[i] *= c_c;
                c_grad[i] = &c_grad[i] + &c_grad_c;
            }

            // add the internal contribution to the jacobian
            c_grad[i] = &c_grad[i] * c[i];
            let l = mom - &self.qs[i];

            for mu in 0..4 {
                k_grad[mu][mu] += -c[i];

                for nu in 0..4 {
                    k_grad[mu][nu] += l[mu] * -c_grad[i][nu];
                }
            }

            k_int = k_int + &l * -c[i]; // massless part
        }

        let k0 = k_int + k_ext;
        let mut lambda_cycle_sq = 1.0; // maximum lambda value
        let mut lambda_grad = LorentzVector::new();

        for j in 0..n {
            let k0_sq = k0.square();
            let lj = mom - &self.qs[j];
            let k0_lqj = k0.dot(&lj);
            let w = lj.square() - self.masses[j].powi(2);

            let xj = (k0_lqj / k0_sq).powi(2);
            let yj = w / k0_sq;

            let mut xj_grad = LorentzVector::new();
            let mut yj_grad = LorentzVector::new();

            for mu in 0..4 {
                xj_grad[mu] += metric_diag[mu] * 2.0 / k0_lqj * k0[mu];
                yj_grad[mu] += metric_diag[mu] * 2.0 / w * lj[mu];
                for nu in 0..4 {
                    xj_grad[mu] += metric_diag[nu]
                        * (-4.0 / k0_sq * k0[nu] + 2.0 / k0_lqj * lj[nu])
                        * k_grad[nu][mu];
                    yj_grad[mu] += metric_diag[nu] * -2.0 / k0_sq * k0[nu] * k_grad[nu][mu];
                }
            }
            xj_grad = &xj_grad * xj;
            yj_grad = &yj_grad * yj;

            if 2.0 * xj < yj {
                if yj / 4.0 < lambda_cycle_sq {
                    lambda_cycle_sq = yj * 0.25;
                    lambda_grad = &yj_grad * 0.25;
                }
            } else if yj < 0.0 {
                if xj - yj / 2.0 < lambda_cycle_sq {
                    lambda_cycle_sq = xj - yj * 0.5;
                    lambda_grad = &xj_grad - &(&yj_grad * 0.5);
                }
            } else {
                if xj - yj / 4.0 < lambda_cycle_sq {
                    lambda_cycle_sq = xj - yj * 0.25;
                    lambda_grad = &xj_grad - &(&yj_grad * 0.25);
                }
            }
        }

        let mut lambda_cycle = lambda_cycle_sq.sqrt();
        lambda_grad = &lambda_grad * (0.5 / lambda_cycle);

        // collinear contribution
        let n1: f64 = c.iter().sum();
        if 1.0 / (4.0 * n1) < lambda_cycle {
            lambda_cycle = 1.0 / (4.0 * n1);

            lambda_grad = LorentzVector::new();
            for x in &c_grad {
                lambda_grad = &lambda_grad + x;
            }

            lambda_grad = &lambda_grad * (-4.0 * lambda_cycle * lambda_cycle);
        }

        let l = mom - &self.uv_shift;
        let uv_fac = 4.0 * l.dot(&k0);
        if uv_fac <= self.mu_sq.im {
            if lambda_cycle > self.mu_sq.im / uv_fac {
                lambda_cycle = self.mu_sq.im / uv_fac;

                lambda_grad = k0.dual();
                for mu in 0..4 {
                    for nu in 0..4 {
                        // TODO: check if the indices are correct!
                        lambda_grad[mu] += metric_diag[nu] * k_grad[nu][mu] * l[nu];
                    }
                }

                lambda_grad = &lambda_grad * (-4.0 * self.mu_sq.im / uv_fac / uv_fac);
            }
        }

        let k = &mom.to_complex(true) + &(&k0 * lambda_cycle).to_complex(false);

        let mut grad = [[Complex::new(0.0, 0.0); 4]; 4];
        for mu in 0..4 {
            grad[mu][mu] += Complex::new(1.0, 0.0);
            for nu in 0..4 {
                grad[mu][nu] += Complex::new(
                    0.0,
                    lambda_cycle * k_grad[mu][nu] + lambda_grad[nu] * k0[mu],
                );
            }
        }

        let jac = determinant(&grad);
        Ok((k, jac))
    }

    fn deform_int_no_jacobian(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<LorentzVector<Complex>, &str> {
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

        // internal part
        let k_centre = &(&k_plus + &k_minus) * 0.5;

        // TODO: cache qi= mom-qs[i]?

        let n = self.qs.len();
        assert!(n <= 10);
        let f = Deformer::g(&k_centre, self.gamma1, self.m2_sq);
        let mut c = [f; 10]; // use static array, for performance

        let mut k_int = LorentzVector::new();
        for i in 0..n {
            for l in 0..n {
                // note: in paper l starts at 1
                c[i] *= Deformer::d(self, mom, i, l);
            }

            let l = mom - &self.qs[i];
            k_int = k_int + &l * -c[i]; // massless part
        }

        let k0 = k_int + k_ext;
        let mut lambda_cycle_sq = 1.0; // maximum lambda value

        for j in 0..n {
            let k0_sq = k0.square();
            let lj = mom - &self.qs[j];
            let k0_lqj = k0.dot(&lj);
            let w = lj.square() - self.masses[j].powi(2);

            let xj = (k0_lqj / k0_sq).powi(2);
            let yj = w / k0_sq;

            if 2.0 * xj < yj {
                if yj / 4.0 < lambda_cycle_sq {
                    lambda_cycle_sq = yj * 0.25;
                }
            } else if yj < 0.0 {
                if xj - yj / 2.0 < lambda_cycle_sq {
                    lambda_cycle_sq = xj - yj * 0.5;
                }
            } else {
                if xj - yj / 4.0 < lambda_cycle_sq {
                    lambda_cycle_sq = xj - yj * 0.25;
                }
            }
        }

        let mut lambda_cycle = lambda_cycle_sq.sqrt();

        // collinear contribution
        let n1: f64 = c.iter().sum();
        if 1.0 / (4.0 * n1) < lambda_cycle {
            lambda_cycle = 1.0 / (4.0 * n1);
        }

        let l = mom - &self.uv_shift;
        let uv_fac = 4.0 * l.dot(&k0);
        if uv_fac <= self.mu_sq.im {
            if lambda_cycle > self.mu_sq.im / uv_fac {
                lambda_cycle = self.mu_sq.im / uv_fac;
            }
        }

        let k = &mom.to_complex(true) + &(&k0 * lambda_cycle).to_complex(false);
        Ok(k)
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

        Ok((k, Complex::new(0.0, -4.0)))
    }

    #[inline]
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

    #[inline]
    pub fn deform_two_loops(
        &mut self,
        id: usize,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, LorentzVector<Complex>), &str> {
        match id {
            DOUBLE_BOX_ID => self.deform_doublebox(k, l),
            DOUBLE_TRIANGLE_ID => self.deform_doubletriangle(k, l),
            _ => Err("Unknown id"),
        }
    }

    fn deform_doubletriangle(
        &mut self,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, LorentzVector<Complex>), &str> {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - l, k - l + p1, k + p1
        let mut c12_qs = [LorentzVector::new(), k - l, -&self.ext[0]];

        self.set_qs(&c12_qs);
        let c12_k = self.deform_int_no_jacobian(k).unwrap().imag(); // get the direction

        c12_qs = [LorentzVector::new(), -k + l - &self.ext[0], -k + l];

        self.set_qs(&c12_qs);
        let c12_l = self.deform_int_no_jacobian(l).unwrap().imag();

        let c23_qs = [LorentzVector::new(), k + &self.ext[0], k.clone()];
        self.set_qs(&c23_qs);
        let c23_l = self.deform_int_no_jacobian(l).unwrap().imag();

        let c13_qs = [
            LorentzVector::new(),
            l.clone(),
            l - &self.ext[0],
            -&self.ext[0],
        ];
        self.set_qs(&c13_qs);
        let c13_k = self.deform_int_no_jacobian(k).unwrap().imag();

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        let mut lambda_sq = 1.;

        // propagators with substituted momenta split in real and imag
        // k, l, k - l, k - l + p1, k + p1
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l, &k_1 - &k_2),
            (k - l + &self.ext[0], &k_1 - &k_2),
            (k + &self.ext[0], k_1.clone()),
        ];

        for (qt, kt) in &props {
            let xj = (kt.dot(qt) / kt.square()).powi(2);
            let yj = (qt.square()) / kt.square(); // for the massless case

            if 2.0 * xj < yj {
                if yj / 4.0 < lambda_sq {
                    lambda_sq = yj * 0.25;
                }
            } else if yj < 0.0 {
                if xj - yj / 2.0 < lambda_sq {
                    lambda_sq = xj - yj * 0.5;
                }
            } else {
                if xj - yj / 4.0 < lambda_sq {
                    lambda_sq = xj - yj * 0.25;
                }
            }
        }
        let lambda = lambda_sq.sqrt();

        let k_1_full = &((&k_1 * lambda).to_complex(false)) + k;
        let k_2_full = &((&k_2 * lambda).to_complex(false)) + l;

        Ok((k_1_full, k_2_full))
    }

    fn deform_doublebox(
        &mut self,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, LorentzVector<Complex>), &str> {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - l, k - l - p2, k - l - p2 - p3, k - p2 - p3, k + p1
        let mut c12_qs = [
            LorentzVector::new(),
            k - l,
            &self.ext[1] + &self.ext[2],
            -&self.ext[0],
        ];

        self.set_qs(&c12_qs);
        let c12_k = self.deform_int_no_jacobian(k).unwrap().imag(); // get the direction

        c12_qs = [
            LorentzVector::new(),
            -k + l + &self.ext[1] + &self.ext[2],
            -k + l - &self.ext[0],
            -k + l,
        ];

        self.set_qs(&c12_qs);
        let c12_l = self.deform_int_no_jacobian(l).unwrap().imag();

        let c23_qs = [
            LorentzVector::new(),
            k - &self.ext[1] - &self.ext[2],
            k - &self.ext[1],
            k.clone(),
        ];
        self.set_qs(&c23_qs);
        let c23_l = self.deform_int_no_jacobian(l).unwrap().imag();

        let c13_qs = [
            LorentzVector::new(),
            l.clone(),
            l + &self.ext[1],
            l + &self.ext[1] + &self.ext[2],
            &self.ext[1] + &self.ext[2],
            -&self.ext[0],
        ];
        self.set_qs(&c13_qs);
        let c13_k = self.deform_int_no_jacobian(k).unwrap().imag();

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        let mut lambda_sq = 1.;

        // propagators with substituted momenta split in real and imag
        // k, l, k - l, k - l - p2, k - l - p2 - p3, k - p2 - p3, k + p1
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l, &k_1 - &k_2),
            (k - l - &self.ext[1], &k_1 - &k_2),
            (k - l - &self.ext[1] - &self.ext[2], &k_1 - &k_2),
            (k - &self.ext[1] - &self.ext[2], k_1.clone()),
            (k + &self.ext[0], k_1.clone()),
        ];

        for (qt, kt) in &props {
            let xj = (kt.dot(qt) / kt.square()).powi(2);
            let yj = (qt.square()) / kt.square(); // for the massless case

            if 2.0 * xj < yj {
                if yj / 4.0 < lambda_sq {
                    lambda_sq = yj * 0.25;
                }
            } else if yj < 0.0 {
                if xj - yj / 2.0 < lambda_sq {
                    lambda_sq = xj - yj * 0.5;
                }
            } else {
                if xj - yj / 4.0 < lambda_sq {
                    lambda_sq = xj - yj * 0.25;
                }
            }
        }
        let lambda = lambda_sq.sqrt();

        let k_1_full = &((&k_1 * lambda).to_complex(false)) + k;
        let k_2_full = &((&k_2 * lambda).to_complex(false)) + l;

        Ok((k_1_full, k_2_full))
    }

    pub fn numerical_jacobian(
        &self,
        k: &LorentzVector<f64>,
        center: &LorentzVector<Complex>,
    ) -> Complex {
        let eps = EPSILON.sqrt();
        let mut grad = [[Complex::new(0., 0.); 4]; 4];

        for i in 0..4 {
            let mut ep = eps;
            if k[i] > 1.0 {
                ep = k[i] * eps;
            }
            let k_p = k + &DIRECTIONS[i] * ep;

            let (res_k, _) = self.deform(&k_p).unwrap();
            let delta_k = (res_k - center) * Complex::new(1. / ep, 0.);
            for j in 0..4 {
                grad[i][j] = delta_k[j];
            }
        }

        determinant(&grad)
    }

    pub fn numerical_jacobian_two_loops(
        &mut self,
        id: usize,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
        center: (&LorentzVector<Complex>, &LorentzVector<Complex>),
    ) -> Complex {
        let eps = EPSILON.sqrt();
        let mut grad = [[Complex::new(0., 0.); 8]; 8];

        for i in 0..8 {
            let mut ep = eps;
            let mut k_n = k.clone();
            let mut l_n = l.clone();
            if i < 4 {
                if k_n[i] > 1.0 {
                    ep = k_n[i] * eps;
                }
                k_n = k_n + &DIRECTIONS[i] * ep;
            } else {
                if l_n[i - 4] > 1.0 {
                    ep = l_n[i - 4] * eps;
                }
                l_n = l_n + &DIRECTIONS[i - 4] * ep;
            }

            let (res_k, res_l) = self.deform_two_loops(id, &k_n, &l_n).unwrap();
            let delta_k = (res_k - center.0) * Complex::new(1. / ep, 0.);
            let delta_l = (res_l - center.1) * Complex::new(1. / ep, 0.);
            for j in 0..4 {
                grad[i][j] = delta_k[j];
                grad[i][j + 4] = delta_l[j];
            }
        }

        determinant8(&grad)
    }

    pub fn numerical_jacobian_center_two_loops(
        &mut self,
        id: usize,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
        center: (&LorentzVector<Complex>, &LorentzVector<Complex>),
    ) -> (Complex, Complex) {
        let eps = EPSILON.sqrt();
        let mut grad = [[Complex::new(0., 0.); 8]; 8];

        for i in 0..8 {
            let mut ep = eps;

            let (res_n_k, res_n_l, res_p_k, res_p_l) = if i < 4 {
                if k[i].abs() > 1.0 {
                    ep = k[i].abs() * eps;
                }
                let k_n = k - &DIRECTIONS[i] * ep;
                let k_p = k + &DIRECTIONS[i] * ep;

                let (res_n_k, res_n_l) = self.deform_two_loops(id, &k_n, l).unwrap();
                let (res_p_k, res_p_l) = self.deform_two_loops(id, &k_p, l).unwrap();
                (res_n_k, res_n_l, res_p_k, res_p_l)
            } else {
                if l[i - 4].abs() > 1.0 {
                    ep = l[i - 4].abs() * eps;
                }
                let l_n = l - &DIRECTIONS[i - 4] * ep;
                let l_p = l + &DIRECTIONS[i - 4] * ep;

                let (res_n_k, res_n_l) = self.deform_two_loops(id, k, &l_n).unwrap();
                let (res_p_k, res_p_l) = self.deform_two_loops(id, k, &l_p).unwrap();
                (res_n_k, res_n_l, res_p_k, res_p_l)
            };

            let delta_k = (res_p_k - res_n_k) * Complex::new(1. / ep / 2., 0.);
            let delta_l = (res_p_l - res_n_l) * Complex::new(1. / ep / 2., 0.);

            // TODO: component-wise error
            //let rel_error_k = max((res_p_k - center[0]).abs(), (res_n_k - center[0]).abs()) * ep / delta_k;
            //let rel_error_l = max((res_p_l - center[1]).abs(), (res_n_l - center[1]).abs()) * ep / delta_l;
            for j in 0..4 {
                grad[i][j] = delta_k[j];
                grad[i][j + 4] = delta_l[j];
            }
        }

        (determinant8(&grad), Complex::new(0., 0.))
    }
}
