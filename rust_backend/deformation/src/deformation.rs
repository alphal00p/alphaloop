use arrayvec::ArrayVec;
use dual_num::dual9::Dual9;
use dual_num::Dual;
use num::Float;
use std::convert::From;
use std::f64::EPSILON;
use utils::{determinant, determinant8};
use vector::LorentzVector;
use vector::RealField;
use Complex;
use REGION_EXT;

pub const Q_UB: usize = 8; // upper bound for number of propagators

pub const DOUBLE_BOX_ID: usize = 0;
pub const DOUBLE_TRIANGLE_ID: usize = 1;
pub const TRIANGLE_BOX_ID: usize = 2;
pub const TRIANGLE_BOX_ALTERNATIVE_ID: usize = 3;
pub const DIAGONAL_BOX_ID: usize = 4;
pub const DOUBLE_BOX_SB_ID: usize = 5;

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

#[derive(Clone)]
pub struct Deformer<F: Float + RealField> {
    qs: ArrayVec<[LorentzVector<F>; Q_UB]>,
    ext: ArrayVec<[LorentzVector<F>; Q_UB]>, // external momenta
    masses: ArrayVec<[f64; Q_UB]>,
    p_plus: LorentzVector<F>,
    p_min: LorentzVector<F>,
    e_cm_sq: f64,
    mu_p_sq: F,
    m1_fac: f64,
    m2_fac: f64,
    m3_fac: f64,
    m4_fac: f64,
    m1_sq: F,
    m2_sq: F,
    m3_sq: F,
    m4_sq: F,
    gamma1: f64,
    gamma2: f64,
    soft_fac: f64,
    e_soft: F,
    region: usize,
    mu_sq: Complex,
    uv_shift: LorentzVector<F>,
}

impl<F: Float + RealField> Deformer<F> {
    pub fn new(
        e_cm_sq: f64,
        mu_sq: f64,
        m1_fac: f64,
        m2_fac: f64,
        m3_fac: f64,
        m4_fac: f64,
        gamma1: f64,
        gamma2: f64,
        soft_fac: f64,
        region: usize,
        masses: &[f64],
    ) -> Result<Deformer<F>, &'static str> {
        Ok(Deformer {
            qs: ArrayVec::new(),
            ext: ArrayVec::new(),
            masses: {
                let mut m = ArrayVec::new();
                m.extend(masses.iter().cloned());
                m
            },
            p_plus: LorentzVector::new(),
            p_min: LorentzVector::new(),
            e_cm_sq,
            mu_p_sq: From::from(0.0),
            m1_fac,
            m2_fac,
            m3_fac,
            m4_fac,
            m1_sq: From::from(0.0),
            m2_sq: From::from(0.0),
            m3_sq: From::from(0.0),
            m4_sq: From::from(0.0),
            gamma1,
            gamma2,
            soft_fac,
            e_soft: From::from(soft_fac * e_cm_sq.sqrt()),
            mu_sq: Complex::new(0.0, mu_sq),
            uv_shift: LorentzVector::new(),
            region,
        })
    }

    /// Set external momenta. Only used for the double box at the moment.
    pub fn set_external_momenta_iter<T>(&mut self, ext: T)
    where
        T: Iterator<Item = LorentzVector<F>>,
    {
        self.ext.clear();
        self.ext.extend(ext);
    }

    pub fn set_external_momenta(&mut self, ext: &[LorentzVector<F>]) {
        self.set_external_momenta_iter(ext.iter().cloned());
    }

    pub fn set_region(&mut self, region: usize) {
        self.region = region;
    }

    fn set_qs_iter<T>(&mut self, qs: T)
    where
        T: Iterator<Item = LorentzVector<F>>,
    {
        self.qs.clear();
        self.qs.extend(qs);

        if self.masses.len() < self.qs.len() {
            let ld = self.qs.len() - self.masses.len();
            self.masses.extend((0..ld).map(|_| 0.));
        }

        self.uv_shift = LorentzVector::new();
        for q in &self.qs {
            self.uv_shift = &self.uv_shift + q;
        }
        self.uv_shift = &self.uv_shift * From::from(1.0 / (self.qs.len() as f64));

        self.p_plus = self.compute_p(true);
        self.p_min = self.compute_p(false);
        self.shift_and_check_p();

        self.mu_p_sq = (&self.p_min - &self.p_plus).square();

        //println!("mu_p={}, e_cm_sq={}", self.mu_p_sq, self.e_cm_sq);
        //self.mu_p_sq = From::from(self.e_cm_sq);

        self.m1_sq = From::from(self.m1_fac * self.m1_fac * self.e_cm_sq);
        self.m2_sq = (if self.mu_p_sq > self.e_cm_sq {
            self.mu_p_sq
        } else {
            From::from(self.e_cm_sq)
        }) * self.m2_fac
            * self.m2_fac;
        self.m3_sq = (if self.mu_p_sq > self.e_cm_sq {
            self.mu_p_sq
        } else {
            From::from(self.e_cm_sq)
        }) * self.m3_fac
            * self.m3_fac;
        self.m4_sq = From::from(self.m4_fac * self.m4_fac * self.e_cm_sq);
    }

    /// Set new qs and update all parameters accordingly.
    pub fn set_qs(&mut self, qs: &[LorentzVector<F>]) {
        self.set_qs_iter(qs.iter().cloned());
    }

    /// Interpolate between vectors
    #[inline]
    fn z(a: LorentzVector<F>, b: LorentzVector<F>, plus: bool) -> LorentzVector<F> {
        if b.euclidean_square() < 1e-9 {
            return a;
        }

        let n = b.spatial_distance().inv();

        if plus {
            &LorentzVector::from_args(
                a.t + n * b.square() - n * b.t * b.t,
                a.x - n * b.t * b.x,
                a.y - n * b.t * b.y,
                a.z - n * b.t * b.z,
            ) * From::from(0.5)
        } else {
            &LorentzVector::from_args(
                a.t - n * b.square() + n * b.t * b.t,
                a.x + n * b.t * b.x,
                a.y + n * b.t * b.y,
                a.z + n * b.t * b.z,
            ) * From::from(0.5)
        }
    }

    ///Find a vector P such that all external momenta `qs` are in
    ///the forward lightcone of P if `plus`, or are in the backwards lightcone
    ///if not `plus`.
    ///This function uses the algorithm from Becker's thesis
    fn compute_p(&mut self, plus: bool) -> LorentzVector<F> {
        let mut p = self.qs[0].clone();

        for q in &self.qs[1..] {
            let s = (p - q).square();
            if s <= 0. {
                p = Deformer::z(p + q, p - q, plus);
            } else {
                if (plus && p.t > q.t) || (!plus && p.t < q.t) {
                    p = q.clone();
                }
            }
        }
        p
    }

    /// Shift P+ and P- outward and check if their new values fulfills P^2 >= 0.0 and (+/-) * P+/-.t <= 0.0  
    /// This is necessary to counterbalance numerical error arising in the exact algorithm
    fn shift_and_check_p(&mut self) {
        //Perform the shift
        let shift_size = 1e-2;
        let r: F = (self.p_min.t - self.p_plus.t) * 0.5;
        let center: F = (self.p_min.t + self.p_plus.t) * 0.5;
        self.p_min.t = center + r * (1.0 + shift_size);
        self.p_plus.t = center - r * (1.0 + shift_size);

        //Check if P+/- has all the qs in its backward/forward light-cone
        let mut pq_diff;
        for q in &self.qs {
            pq_diff = &self.p_plus - q;
            if pq_diff.t > 0.0 || pq_diff.square() < 0.0 {
                panic!(
                    "P_plus is not correctly defined! (P-q).t = {}, (P-q)^2 = {}, P = {}, qs = {:#?}",
                    pq_diff.t,
                    pq_diff.square(),
                    self.p_plus,
                    self.qs
                );
            }

            pq_diff = &self.p_min - q;
            if pq_diff.t < 0.0 || pq_diff.square() < 0.0 {
                panic!(
                    "P_minus is not correctly defined! (P-q).t = {}, (P-q)^2 = {}, P = {}, qs = {:#?}",
                    pq_diff.t,
                    pq_diff.square(),
                    self.p_min,
                    self.qs
                );
            }
        }
    }

    /// The three helper functions h_delta-, h_delta+, and h_delta0, indicated by `sign`.
    fn h_delta(sign: i32, k: &LorentzVector<F>, mass: f64, m: F) -> F {
        let a = k.spatial_squared() + <F as num::NumCast>::from(mass).unwrap();
        let v = if sign == 0 {
            (k.t.abs() - a.sqrt()).powi(2)
        } else {
            (k.t * <F as num::NumCast>::from(sign).unwrap() - a.sqrt()).powi(2)
        };
        v / (v + m)
    }

    #[inline]
    fn h_theta(t: F, m: F) -> F {
        if t <= 0. {
            num::NumCast::from(0.).unwrap()
        } else {
            t / (t + m)
        }
    }

    #[inline]
    fn g(k: &LorentzVector<F>, gamma: F, m: F) -> F {
        (k.euclidean_square() + m).inv() * gamma * m
    }

    /// Return d
    fn d(&self, mom: &LorentzVector<F>, i: usize, l: usize) -> F {
        if self.masses[l] == 0.0 {
            if l == i {
                return num::NumCast::from(1.0).unwrap();
            }

            if (&self.qs[i] - &self.qs[l]).square() == 0.0 {
                let a: LorentzVector<F> = mom - self.qs[l];
                if self.qs[i].t < self.qs[l].t {
                    return Deformer::h_delta(1, &a, 0.0, self.m1_sq);
                } else {
                    return Deformer::h_delta(-1, &a, 0.0, self.m1_sq);
                }
            }
        }

        let a: LorentzVector<F> = mom - self.qs[l];
        let hd = Deformer::h_delta(0, &a, self.masses[l] * self.masses[l], self.m1_sq);
        let ht = Deformer::h_theta(a.dot(&(mom - self.qs[i])) * -2.0, self.m1_sq);
        hd.max(ht)
    }

    fn deform_int_no_jacobian(
        &self,
        mom: &LorentzVector<F>,
        include_cij: bool,
        include_ca: bool,
        include_uv: bool,
    ) -> LorentzVector<F> {
        let (mut c_plus, mut c_minus): (F, F) = (From::from(1.), From::from(1.));

        for (qi, mi) in self.qs.iter().zip(&self.masses) {
            // note the sign reversal
            c_plus *= Deformer::h_delta(-1, &(mom - qi), mi * mi, self.m3_sq);
            c_minus *= Deformer::h_delta(1, &(mom - qi), mi * mi, self.m3_sq);
        }

        let k_plus = mom - self.p_plus;
        let k_minus = mom - self.p_min;
        let k_ext = LorentzVector::from_args(
            c_plus * k_plus.t + c_minus * k_minus.t,
            -c_plus * k_plus.x - c_minus * k_minus.x,
            -c_plus * k_plus.y - c_minus * k_minus.y,
            -c_plus * k_plus.z - c_minus * k_minus.z,
        );

        // internal part
        let k_centre = (k_plus + k_minus) * From::from(0.5);

        let n = self.qs.len();
        assert!(n <= Q_UB);
        let f = Deformer::g(&k_centre, From::from(self.gamma1), self.m2_sq);

        // the sum of all the cs, used for lambda scaling
        let mut c_sum: F = From::from(0.);
        let mut k_int = LorentzVector::new();

        // compute c
        for i in 0..n {
            let mut c = f;
            for l in 0..n {
                // note: in paper l starts at 1
                c *= Deformer::d(self, mom, i, l);
            }

            let l = mom - self.qs[i];
            k_int = k_int + l * -c;
            c_sum += c;
        }

        // compute cij
        if include_cij {
            for i in 0..n {
                for j in i + 1..n {
                    let mut c = f;

                    let zij = (self.qs[i] - self.qs[j]).square();
                    if zij > 0. {
                        let kij = if true {
                            // multi-loop paper approach
                            // TODO: check if we should drop the ^n
                            c *= Deformer::h_theta(zij, self.m1_sq).powi(n as i32);
                            mom - self.qs[i] * From::from(0.5) - self.qs[j] * From::from(0.5)
                        } else {
                            // thesis approach
                            let ki = mom - self.qs[i];
                            let kj = mom - self.qs[j];
                            (mom - (self.qs[i] * kj.t - self.qs[j] * ki.t)
                                * (self.qs[i].t - self.qs[j].t).inv())
                                * (zij / (zij + self.m2_sq))
                        };

                        for l in 0..n {
                            let l1 = Deformer::h_delta(
                                0,
                                &(mom - self.qs[l]),
                                self.masses[l],
                                self.m1_sq,
                            );

                            let a = (mom - self.qs[l]).dot(&kij) * -2.;
                            let l2 = Deformer::h_theta(a, self.m1_sq);

                            c *= l1.max(l2);
                        }

                        k_int = k_int - kij * c;
                        c_sum += c;
                    }
                }
            }
        }

        // compute ca
        if include_ca {
            let nul = num::NumCast::from(0.).unwrap();
            let ka = [
                LorentzVector::from_args(self.e_soft, nul, nul, nul),
                LorentzVector::from_args(nul, self.e_soft, nul, nul),
                LorentzVector::from_args(nul, nul, self.e_soft, nul),
                LorentzVector::from_args(nul, nul, nul, self.e_soft),
            ];

            for a in 0..4 {
                let mut c = f;
                let mut da_plus: F = num::NumCast::from(1.).unwrap();
                let mut da_min: F = num::NumCast::from(1.).unwrap();
                for l in 0..n {
                    let l1 =
                        Deformer::h_delta(0, &(mom - self.qs[l]), 0., self.m1_sq * self.gamma2);
                    let lp2 = Deformer::h_theta(
                        (mom - self.qs[l]).dot(&ka[a]) * 2.,
                        self.m1_sq * self.gamma2,
                    );
                    let lm2 = Deformer::h_theta(
                        (mom - self.qs[l]).dot(&ka[a]) * -2.,
                        self.m1_sq * self.gamma2,
                    );

                    da_plus *= l1.max(lp2);
                    da_min *= l1.max(lm2);
                }

                c *= da_plus - da_min;
                k_int = k_int + ka[a] * c;
                c_sum += c.abs();
            }
        }

        let k0 = k_int + k_ext;
        let mut lambda_cycle_sq: F = From::from(1.); // maximum lambda value

        let mut restriction = false;
        let k0_sq = k0.square();
        for j in 0..n {
            let lj = mom - self.qs[j];
            let k0_lqj = k0.dot(&lj);
            let w = lj.square() - self.masses[j].powi(2);

            let xj = (k0_lqj / k0_sq).powi(2);
            let yj = w / k0_sq;

            if xj * 2. < yj {
                if yj * 0.25 < lambda_cycle_sq {
                    lambda_cycle_sq = yj * 0.25;
                }
            } else if yj < 0. {
                if xj - yj * 0.5 < lambda_cycle_sq {
                    lambda_cycle_sq = xj - yj * 0.5;
                    restriction = true;
                }
            } else {
                if xj - yj * 0.25 < lambda_cycle_sq {
                    lambda_cycle_sq = xj - yj * 0.25;
                    restriction = true;
                }
            }
        }

        let mut lambda_cycle = lambda_cycle_sq.sqrt();

        // collinear contribution
        if c_sum.inv() * 0.25 < lambda_cycle {
            if restriction {
                // TODO: we are not near any collinear limit. We could consider not scaling the lambda down!
                //println!("We are not near the collinear limit, yet we rescale lambda from {} to {}", lambda_cycle, c_sum.inv() * 0.25);
            }

            //println!("Coll scaling from {} to {}", lambda_cycle, c_sum.inv() * 0.25);
            lambda_cycle = c_sum.inv() * 0.25;
        }

        if include_uv {
            let l = mom - self.uv_shift;
            let uv_fac: F = l.dot(&k0) * 4.;
            if uv_fac <= self.mu_sq.im {
                if lambda_cycle > uv_fac.inv() * self.mu_sq.im {
                    lambda_cycle = uv_fac.inv() * self.mu_sq.im;
                }
            }

            // each propagator could be UV since it regularises the collinear CT
            for q in &self.qs {
                let uv_fac: F = (mom - q).dot(&k0) * 4.;
                if uv_fac <= self.mu_sq.im {
                    if lambda_cycle > uv_fac.inv() * self.mu_sq.im {
                        lambda_cycle = uv_fac.inv() * self.mu_sq.im;
                    }
                }
            }
        }

        k0 * lambda_cycle
    }

    fn deform_two_loops_no_jac(
        &mut self,
        id: usize,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> Result<(LorentzVector<F>, LorentzVector<F>), &'static str> {
        match id {
            DOUBLE_BOX_ID => Ok(self.deform_doublebox(k, l)),
            DOUBLE_TRIANGLE_ID => Ok(self.deform_doubletriangle(k, l)),
            TRIANGLE_BOX_ID => Ok(self.deform_trianglebox(k, l)),
            TRIANGLE_BOX_ALTERNATIVE_ID => Ok(self.deform_trianglebox_alternative(k, l)),
            DIAGONAL_BOX_ID => Ok(self.deform_diagonalbox(k, l)),
            DOUBLE_BOX_SB_ID => Ok(self.deform_doublebox_sb(k, l)),
            _ => Err("Unknown id"),
        }
    }

    #[inline]
    fn compute_overall_lambda(
        &self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
        k_1: &LorentzVector<F>,
        k_2: &LorentzVector<F>,
        props: &[(LorentzVector<F>, LorentzVector<F>)],
    ) -> F {
        let mut lambda_sq = From::from(1.);

        for (qt, kt) in props {
            let xj = (kt.dot(qt) / kt.square()).powi(2);
            let yj = (qt.square()) / kt.square(); // for the massless case

            if xj * 2.0 < yj {
                if yj * 0.25 < lambda_sq {
                    lambda_sq = yj * 0.25;
                }
            } else if yj < 0. {
                if xj - yj * 0.5 < lambda_sq {
                    lambda_sq = xj - yj * 0.5;
                }
            } else {
                if xj - yj * 0.25 < lambda_sq {
                    lambda_sq = xj - yj * 0.25;
                }
            }
        }

        // we assume the UV propagators to be 1/(k^2+mu^2) and 1/(l^2+mu^2)
        let mut lambda = lambda_sq.sqrt();

        let uv_fac_k = k.dot(&k_1) * 4.0;
        let uv_fac_l = l.dot(&k_2) * 4.0;
        if uv_fac_k <= self.mu_sq.im && lambda > uv_fac_k.inv() * self.mu_sq.im {
            lambda = uv_fac_k.inv() * self.mu_sq.im;
        }
        if uv_fac_l <= self.mu_sq.im && lambda > uv_fac_l.inv() * self.mu_sq.im {
            lambda = uv_fac_l.inv() * self.mu_sq.im;
        }

        lambda
    }

    fn deform_doubletriangle(
        &mut self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> (LorentzVector<F>, LorentzVector<F>) {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - l, k - l + p1, k + p1
        let mut c12_qs = [LorentzVector::new(), k - l, -&self.ext[0]];

        self.set_qs(&c12_qs);
        let c12_k = self.deform_int_no_jacobian(k, true, false, false); // get the direction

        c12_qs = [LorentzVector::new(), -k + l - &self.ext[0], -k + l];

        self.set_qs(&c12_qs);
        let c12_l = self.deform_int_no_jacobian(l, true, false, false);

        let c23_qs = [LorentzVector::new(), k + &self.ext[0], k.clone()];
        self.set_qs(&c23_qs);
        let c23_l = self.deform_int_no_jacobian(l, true, false, false);

        let c13_qs = [
            LorentzVector::new(),
            l.clone(),
            l - &self.ext[0],
            -&self.ext[0],
        ];
        self.set_qs(&c13_qs);
        let c13_k = self.deform_int_no_jacobian(k, true, false, false);

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        // propagators with substituted momenta split in real and imag
        // k, l, k - l, k - l + p1, k + p1
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l, &k_1 - &k_2),
            (k - l + &self.ext[0], &k_1 - &k_2),
            (k + &self.ext[0], k_1.clone()),
        ];

        let lambda = self.compute_overall_lambda(k, l, &k_1, &k_2, &props);
        (k_1 * lambda, k_2 * lambda)
    }

    fn deform_trianglebox_alternative(
        &mut self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> (LorentzVector<F>, LorentzVector<F>) {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - l, l - p2, l + p1, k + p1
        let momenta_12 = [
            k.clone(),
            k - self.ext[0],
            l.clone(),
            l - self.ext[1],
            l + self.ext[0],
        ];
        let momenta_13 = [k - self.ext[0], k.clone(), k - l - self.ext[0]];
        let momenta_23 = [
            l.clone(),
            l - self.ext[1],
            l + self.ext[0],
            l - k + self.ext[0],
        ];

        //Compute qs
        let c12_qs_k: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_12.iter().map(|mom| k - mom).collect();
        let c12_qs_l: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_12.iter().map(|mom| l - mom).collect();
        let c13_qs_k: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_13.iter().map(|mom| k - mom).collect();
        let c23_qs_l: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_23.iter().map(|mom| l - mom).collect();

        //Compute kappas
        self.set_qs(&c12_qs_k);
        let c12_k = self.deform_int_no_jacobian(k, true, false, false);

        self.set_qs(&c12_qs_l);
        let c12_l = self.deform_int_no_jacobian(l, true, false, false);

        self.set_qs(&c13_qs_k);
        let c13_k = self.deform_int_no_jacobian(k, true, false, false);

        self.set_qs(&c23_qs_l);
        let c23_l = self.deform_int_no_jacobian(l, true, false, false);

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        // propagators with substituted momenta split in real and imag
        // k, l, k - l, k - l - p2, k - l + p1, k + p1
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l - self.ext[0], &k_1 - &k_2),
            (l - self.ext[1], k_2.clone()),
            (l + self.ext[0], k_2.clone()),
            (k - self.ext[0], k_1.clone()),
        ];

        let lambda = self.compute_overall_lambda(k, l, &k_1, &k_2, &props);
        (k_1 * lambda, k_2 * lambda)
    }

    fn deform_trianglebox(
        &mut self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> (LorentzVector<F>, LorentzVector<F>) {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - l, k - l - p2, k - l + p1, k + p1
        let mut c12_qs = [LorentzVector::new(), k - l, -&self.ext[0]];

        self.set_qs(&c12_qs);
        let c12_k = self.deform_int_no_jacobian(k, true, false, false); // get the direction

        c12_qs = [LorentzVector::new(), -k + l - &self.ext[0], -k + l];

        self.set_qs(&c12_qs);
        let c12_l = self.deform_int_no_jacobian(l, true, false, false);

        let c23_qs = [
            LorentzVector::new(),
            k + &self.ext[0],
            k - &self.ext[1],
            k.clone(),
        ];
        self.set_qs(&c23_qs);
        let c23_l = self.deform_int_no_jacobian(l, true, false, false);

        let c13_qs = [
            LorentzVector::new(),
            l.clone(),
            l + &self.ext[1],
            l - &self.ext[0],
            -&self.ext[0],
        ];
        self.set_qs(&c13_qs);
        let c13_k = self.deform_int_no_jacobian(k, true, false, false);

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        // propagators with substituted momenta split in real and imag
        // k, l, k - l, k - l - p2, k - l + p1, k + p1
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l, &k_1 - &k_2),
            (k - l - &self.ext[1], &k_1 - &k_2),
            (k - l + &self.ext[0], &k_1 - &k_2),
            (k + &self.ext[0], k_1.clone()),
        ];

        let lambda = self.compute_overall_lambda(k, l, &k_1, &k_2, &props);
        (k_1 * lambda, k_2 * lambda)
    }

    fn deform_diagonalbox(
        &mut self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> (LorentzVector<F>, LorentzVector<F>) {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - p1 - p2 - l, k - p1, l - p3,
        let momenta_12 = [k.clone(), k - &self.ext[0], l.clone(), l - &self.ext[2]];
        let momenta_13 = [
            k.clone(),
            k - &self.ext[0],
            k - l - &self.ext[0] - &self.ext[1],
        ];
        let momenta_23 = [
            l.clone(),
            l - &self.ext[2],
            l - k + &self.ext[0] + &self.ext[1],
        ];

        //Compute qs
        let c12_qs_k: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_12.iter().map(|mom| k - mom).collect();
        let c12_qs_l: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_12.iter().map(|mom| l - mom).collect();
        let c13_qs_k: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_13.iter().map(|mom| k - mom).collect();
        let c23_qs_l: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_23.iter().map(|mom| l - mom).collect();

        //Compute kappas
        self.set_qs(&c12_qs_k);
        let c12_k = self.deform_int_no_jacobian(k, true, false, false);

        self.set_qs(&c12_qs_l);
        let c12_l = self.deform_int_no_jacobian(l, true, false, false);

        self.set_qs(&c13_qs_k);
        let c13_k = self.deform_int_no_jacobian(k, true, false, false);

        self.set_qs(&c23_qs_l);
        let c23_l = self.deform_int_no_jacobian(l, true, false, false);

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        // propagators with substituted momenta split in real and imag
        // k, l, k - p1 - p2 - l, k - p1, l - p3,
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l - &self.ext[0] - &self.ext[1], &k_1 - &k_2),
            (k - &self.ext[0], k_1.clone()),
            (l - &self.ext[2], k_2.clone()),
        ];

        let mut lambda = self.compute_overall_lambda(k, l, &k_1, &k_2, &props);

        // for the counterterms, we have two additional UV propagators
        let uv_fac_k = (k - &self.ext[0]).dot(&k_1) * 4.0;
        let uv_fac_l = (l - &self.ext[2]).dot(&k_2) * 4.0;
        if uv_fac_k <= self.mu_sq.im && lambda > uv_fac_k.inv() * self.mu_sq.im {
            lambda = uv_fac_k.inv() * self.mu_sq.im;
        }
        if uv_fac_l <= self.mu_sq.im && lambda > uv_fac_l.inv() * self.mu_sq.im {
            lambda = uv_fac_l.inv() * self.mu_sq.im;
        }

        (k_1 * lambda, k_2 * lambda)
    }

    fn deform_doublebox_sb(
        &mut self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> (LorentzVector<F>, LorentzVector<F>) {
        // First let us set up the momenta q's for all four lambda_cycles

        let (cycle_12_qs, cycle_21_qs, cycle_13_qs, cycle_23_qs) = {
            let p = [&self.ext[0], &self.ext[1], &self.ext[2], &self.ext[3]];
            let zero_q = LorentzVector::<F>::default();
            let cycle_12_qs: ArrayVec<[LorentzVector<F>; 6]> = ArrayVec::from([
                -p[2],
                -p[1] - p[2],
                p[3].clone(),
                k - l + p[3],
                k - l,
                k - l - p[2],
            ]);
            let cycle_21_qs: ArrayVec<[LorentzVector<F>; 6]> = ArrayVec::from([
                p[3].clone(),
                zero_q.clone(),
                -p[2],
                l - k - p[2],
                l - k - p[1] - p[2],
                l - k + p[3],
            ]);
            let mut cycle_13_qs: ArrayVec<[LorentzVector<F>; 6]> = ArrayVec::new();
            cycle_13_qs.extend(
                [-p[2], -p[1] - p[2], p[3].clone(), l.clone()]
                    .iter()
                    .cloned(),
            );
            let mut cycle_23_qs: ArrayVec<[LorentzVector<F>; 6]> = ArrayVec::new();
            cycle_23_qs.extend(
                [p[3].clone(), zero_q.clone(), -p[2], k.clone()]
                    .iter()
                    .cloned(),
            );
            (cycle_12_qs, cycle_21_qs, cycle_13_qs, cycle_23_qs)
        };

        // Perform sanity check
        /*
        for (i, cycle_qs) in [&cycle_12_qs, &cycle_21_qs, &cycle_13_qs, &cycle_23_qs].iter().enumerate() {
            let test_left = cycle_qs
                .iter()
                .fold(LorentzVector::default(), |a, b| a + b)
                .euclidean_square()
                .sqrt();
            let test_right = 1.0e-Q_UB * cycle_qs
                .iter()
                .map(|mom| mom.euclidean_square().sqrt())
                .sum::<F>();
            assert!(
                test_left < test_right,
                "Incorrect definition of qs for cycle #{} ({:#?}). sum(qs) = {:?} > 1e-Q_UB*sum(|qs|) = {:?}",
                i+1, cycle_qs,
                test_left,
                test_right,
            );
        }
        */

        // Now compute the corresponding Kappas
        self.set_qs(&cycle_12_qs);
        let cycle_12_kappa = self.deform_int_no_jacobian(k, true, false, false);
        self.set_qs(&cycle_21_qs);
        let cycle_21_kappa = self.deform_int_no_jacobian(l, true, false, false);
        self.set_qs(&cycle_13_qs);
        let cycle_13_kappa = self.deform_int_no_jacobian(k, true, false, false);
        self.set_qs(&cycle_23_qs);
        let cycle_23_kappa = self.deform_int_no_jacobian(l, true, false, false);

        // Add the kappas of the cycles contributing to each loop momentum
        let kappa_k = cycle_12_kappa + cycle_13_kappa;
        let kappa_l = cycle_21_kappa + cycle_23_kappa;

        // Now address the scaling of the deformation, a.k.a lambda determination
        let mut lambda_sq = From::from(1.0);

        // Listing the propagators of the topology, separating the real and imaginary parts
        let propagators = {
            let p = (&self.ext[0], &self.ext[1], &self.ext[2], &self.ext[3]);
            [
                ((k + p.2), kappa_k.clone()),
                ((k + &(p.1 + p.2)), kappa_k.clone()),
                ((k - p.3), kappa_k.clone()),
                ((l - p.3), kappa_l.clone()),
                (l.clone(), kappa_l.clone()),
                ((l + p.2), kappa_l.clone()),
                ((k - l), (kappa_k - kappa_l)),
            ]
        };

        for (prop_real, prop_imag) in &propagators {
            let xj = (prop_imag.dot(prop_real) / prop_imag.square()).powi(2);
            // Warning, the line below is only correct for the massless case
            let yj = (prop_real.square()) / prop_imag.square();

            if xj * 2.0 < yj {
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

        // UV rescaling is ignored for now
        let lambda = lambda_sq.sqrt();

        (kappa_k * lambda, kappa_l * lambda)
    }

    fn deform_doublebox(
        &mut self,
        k: &LorentzVector<F>,
        l: &LorentzVector<F>,
    ) -> (LorentzVector<F>, LorentzVector<F>) {
        // compute the external momenta for the C12(k) cycle and express the qs in terms of that.
        // This is not a shift. We also set q[0] = 0 at the 1/k^2 line
        // k, l, k - l -p1, l - p2, l - p2 - p3, k + p4, k - p1
        let momenta_12 = [
            k.clone(),
            k - &self.ext[0],
            l.clone(),
            l - &self.ext[1],
            l - &self.ext[1] - &self.ext[2],
            k + &self.ext[3],
        ];
        let momenta_13 = [
            k.clone(),
            k - &self.ext[0],
            k - l - &self.ext[0],
            k + &self.ext[3],
        ];
        let momenta_23 = [
            l.clone(),
            l - &self.ext[1],
            l - &self.ext[1] - &self.ext[2],
            l - k + &self.ext[0],
        ];

        //Compute qs
        let c12_qs_k: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_12.iter().map(|mom| k - mom).collect();
        let c12_qs_l: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_12.iter().map(|mom| l - mom).collect();
        let c13_qs_k: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_13.iter().map(|mom| k - mom).collect();
        let c23_qs_l: ArrayVec<[LorentzVector<F>; Q_UB]> =
            momenta_23.iter().map(|mom| l - mom).collect();

        //Compute kappas
        self.set_qs(&c12_qs_k);
        let c12_k = self.deform_int_no_jacobian(k, true, false, false);

        self.set_qs(&c12_qs_l);
        let c12_l = self.deform_int_no_jacobian(l, true, false, false);

        self.set_qs(&c13_qs_k);
        let c13_k = self.deform_int_no_jacobian(k, true, false, false);

        self.set_qs(&c23_qs_l);
        let c23_l = self.deform_int_no_jacobian(l, true, false, false);

        let k_1 = c12_k + c13_k;
        let k_2 = c12_l + c23_l;

        // propagators with substituted momenta split in real and imag
        // k, l, k - l, k - l - p2, k - l - p2 - p3, k - p2 - p3, k + p1
        let props = [
            (k.clone(), k_1.clone()),
            (l.clone(), k_2.clone()),
            (k - l - &self.ext[0], &k_1 - &k_2),
            (l - &self.ext[1], k_2.clone()),
            (l - &self.ext[1] - &self.ext[2], k_2.clone()),
            (k + &self.ext[3], k_1.clone()),
            (k - &self.ext[0], k_1.clone()),
        ];

        let lambda = self.compute_overall_lambda(k, l, &k_1, &k_2, &props);
        (k_1 * lambda, k_2 * lambda)
    }
}

impl Deformer<f64> {
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

    /// Derivative of h_theta wrt ln t
    #[inline]
    fn h_theta_ln_grad(t: f64, m: f64) -> f64 {
        if t <= 0.0 {
            0.0
        } else {
            m / t / (t + m)
        }
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
                    * (-2.0
                        * Deformer::h_theta_ln_grad(
                            -2.0 * (mom - &self.qs[l]).dot(&(mom - &self.qs[i])),
                            self.m1_sq,
                        )),
            )
        }
    }

    fn deform_int(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, Complex), &'static str> {
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
        let k_ext = LorentzVector::from_args(
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
        assert!(n <= Q_UB);
        let f = Deformer::g(&k_centre, self.gamma1, self.m2_sq);
        let mut c = [f; Q_UB];
        let f_grad = Deformer::g_ln_grad(&k_centre, self.m2_sq);
        let mut c_grad = [f_grad; Q_UB];

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
        let n1: f64 = c[..n].iter().sum();
        if 1.0 / (4.0 * n1) < lambda_cycle {
            lambda_cycle = 1.0 / (4.0 * n1);

            lambda_grad = LorentzVector::new();
            for x in &c_grad[..n] {
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

        // for testing
        if false {
            let (k_new, jac_dual) = self.jacobian_using_dual(mom, false);

            if (k0 * lambda_cycle - k_new).square().abs() > 0.0001 {
                println!("Different deformation: k={}, k_new={}", k0, k_new);
            }

            if (jac_dual.norm() - jac.norm()).abs() > 0.000001 {
                println!(
                    "jac={}, jac_dual={}, k={}, k_new={}",
                    jac, jac_dual, k, k_new
                );
            }
        }

        Ok((k, jac))
    }

    fn deform_ext(
        &self,
        mom: &LorentzVector<f64>,
    ) -> Result<(LorentzVector<Complex>, Complex), &'static str> {
        let k = LorentzVector::from_args(
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
    ) -> Result<(LorentzVector<Complex>, Complex), &'static str> {
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
    ) -> Result<(LorentzVector<Complex>, LorentzVector<Complex>), &'static str> {
        let (kd, ld) = match id {
            DOUBLE_BOX_ID => self.deform_doublebox(k, l),
            DOUBLE_TRIANGLE_ID => self.deform_doubletriangle(k, l),
            TRIANGLE_BOX_ID => self.deform_trianglebox(k, l),
            TRIANGLE_BOX_ALTERNATIVE_ID => self.deform_trianglebox_alternative(k, l),
            DIAGONAL_BOX_ID => self.deform_diagonalbox(k, l),
            DOUBLE_BOX_SB_ID => self.deform_doublebox_sb(k, l),
            _ => return Err("Unknown id"),
        };

        Ok((&kd.to_complex(false) + k, &ld.to_complex(false) + l))
    }

    pub fn jacobian_using_dual(
        &self,
        k: &LorentzVector<f64>,
        use_cij: bool,
    ) -> (LorentzVector<f64>, Complex) {
        let mut dual_deformer: Deformer<Dual<f64>> = Deformer::new(
            self.e_cm_sq,
            self.mu_sq.im,
            self.m1_fac,
            self.m2_fac,
            self.m3_fac,
            self.m4_fac,
            self.gamma1,
            self.gamma2,
            self.soft_fac,
            self.region,
            &self.masses,
        )
        .unwrap();
        dual_deformer.set_qs_iter(self.qs.iter().map(|x| LorentzVector::from_f64(*x)));

        let mut grad = [[Complex::new(0., 0.); 4]; 4];

        let mut k_new = LorentzVector::new();
        for i in 0..4 {
            grad[i][i] += Complex::new(1., 0.); // for the real part

            // we need to decide where to place the 1 for the ep
            let mut k_ep: LorentzVector<Dual<f64>> = LorentzVector::from_f64(*k);
            k_ep[i] = Dual::new(k[i], 1.0);

            // disable ca, like in the one-loop code
            let res_k = dual_deformer.deform_int_no_jacobian(&k_ep, use_cij, false, true);

            k_new = LorentzVector::from_args(
                res_k[0].real(),
                res_k[1].real(),
                res_k[2].real(),
                res_k[3].real(),
            );

            for j in 0..4 {
                grad[i][j] += Complex::new(0.0, res_k[j].dual());
            }
        }

        (k_new, determinant(&grad))
    }

    pub fn jacobian_using_dual9_two_loops(
        &mut self,
        id: usize,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
    ) -> (LorentzVector<Complex>, LorentzVector<Complex>, Complex) {
        // compute with dual9
        let mut dual9_deformer: Deformer<Dual9<f64>> = Deformer::new(
            self.e_cm_sq,
            self.mu_sq.im,
            self.m1_fac,
            self.m2_fac,
            self.m3_fac,
            self.m4_fac,
            self.gamma1,
            self.gamma2,
            self.soft_fac,
            self.region,
            &self.masses,
        ).unwrap();
        dual9_deformer
            .set_external_momenta_iter(self.ext.iter().map(|x| LorentzVector::from_f64(*x)));

        let mut k_dual9: LorentzVector<Dual9<f64>> = LorentzVector::from_f64(*k);
        let mut l_dual9: LorentzVector<Dual9<f64>> = LorentzVector::from_f64(*l);

        for i in 0..8 {
            if i < 4 {
                k_dual9[i][i + 1] = 1.0;
            } else {
                l_dual9[i - 4][i + 1] = 1.0;
            }
        }

        let (k_res, l_res) = dual9_deformer
            .deform_two_loops_no_jac(id, &k_dual9, &l_dual9)
            .unwrap();

        let mut grad = [[Complex::new(0., 0.); 8]; 8];
        for i in 0..8 {
            grad[i][i] += Complex::new(1., 0.);
            for j in 0..8 {
                if i < 4 {
                    grad[i][j] += Complex::new(0.0, k_res[i][j + 1]);
                } else {
                    grad[i][j] += Complex::new(0.0, l_res[i - 4][j + 1]);
                }
            }
        }

        (
            &k_res.map(|x| x.real()).to_complex(false) + k,
            &l_res.map(|x| x.real()).to_complex(false) + l,
            determinant8(&grad),
        )
    }

    pub fn jacobian_using_dual_two_loops(
        &mut self,
        id: usize,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
    ) -> (LorentzVector<Complex>, LorentzVector<Complex>, Complex) {
        // create a deformer that works with dual numbers
        let mut dual_deformer: Deformer<Dual<f64>> = Deformer::new(
            self.e_cm_sq,
            self.mu_sq.im,
            self.m1_fac,
            self.m2_fac,
            self.m3_fac,
            self.m4_fac,
            self.gamma1,
            self.gamma2,
            self.soft_fac,
            self.region,
            &self.masses,
        )
        .unwrap();
        dual_deformer
            .set_external_momenta_iter(self.ext.iter().map(|x| LorentzVector::from_f64(*x)));

        let mut grad = [[Complex::new(0., 0.); 8]; 8];

        let mut k_dual: LorentzVector<Dual<f64>> = LorentzVector::from_f64(*k);
        let mut l_dual: LorentzVector<Dual<f64>> = LorentzVector::from_f64(*l);

        let mut k_res: LorentzVector<Dual<f64>> = LorentzVector::new();
        let mut l_res: LorentzVector<Dual<f64>> = LorentzVector::new();
        for i in 0..8 {
            grad[i][i] += Complex::new(1., 0.); // for the real part

            // we need to decide where to place the 1 for the ep
            if i < 4 {
                k_dual[i] = Dual::new(k[i], 1.0);
            } else {
                l_dual[i - 4] = Dual::new(l[i - 4], 1.0);
            }

            let (kr, lr) = dual_deformer
                .deform_two_loops_no_jac(id, &k_dual, &l_dual)
                .unwrap();
            k_res = kr;
            l_res = lr;

            if i < 4 {
                k_dual[i] = Dual::new(k[i], 0.0);
            } else {
                l_dual[i - 4] = Dual::new(l[i - 4], 0.0);
            }

            for j in 0..4 {
                grad[i][j] += Complex::new(0.0, k_res[j].dual());
                grad[i][j + 4] += Complex::new(0.0, l_res[j].dual());
            }
        }

        (
            &k_res.map(|x| x.real()).to_complex(false) + k,
            &l_res.map(|x| x.real()).to_complex(false) + l,
            determinant8(&grad),
        )
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
        let a = eps * 2.0; // cut-off for when the eps should not be scaled with the value of k

        for i in 0..8 {
            let mut ep = eps;
            let mut k_n = k.clone();
            let mut l_n = l.clone();
            if i < 4 {
                if k_n[i].abs() > a {
                    ep = k_n[i].abs() * eps;
                }
                k_n = k_n + &DIRECTIONS[i] * ep;
            } else {
                if l_n[i - 4].abs() > a {
                    ep = l_n[i - 4].abs() * eps;
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
        let a = eps * 2.0; // cut-off for when the eps should not be scaled with the value of k

        for i in 0..8 {
            let mut ep = eps;

            let (res_n_k, res_n_l, res_p_k, res_p_l) = if i < 4 {
                if k[i].abs() > a {
                    ep = k[i].abs() * eps;
                }
                let k_n = k - &DIRECTIONS[i] * ep;
                let k_p = k + &DIRECTIONS[i] * ep;

                let (res_n_k, res_n_l) = self.deform_two_loops(id, &k_n, l).unwrap();
                let (res_p_k, res_p_l) = self.deform_two_loops(id, &k_p, l).unwrap();
                (res_n_k, res_n_l, res_p_k, res_p_l)
            } else {
                if l[i - 4].abs() > a {
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
