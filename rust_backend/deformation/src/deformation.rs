use num;
use std::mem;
use vector::LorentzVector;

type Complex = num::Complex<f64>;

pub struct Deformer {
    pub qs: Vec<LorentzVector<f64>>,
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
}

impl Deformer {
    pub fn new(e_cm_sq: f64) -> Result<Deformer, &'static str> {
        Ok(Deformer {
            qs: Vec::with_capacity(4),
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
        })
    }

    /// Set new qs and update all parameters accordingly.
    pub fn set_qs(&mut self, qs: Vec<LorentzVector<f64>>) {
        self.qs = qs;

        self.p_plus = self.compute_p(true);
        self.p_min = self.compute_p(false);
        self.mu_p = (&self.p_min - &self.p_plus).square();

        self.m1_sq = self.m1_fac * self.m1_fac * self.e_cm_sq;
        self.m2_sq = self.m1_fac * self.m1_fac * (if self.mu_p > self.e_cm_sq { self.mu_p } else { self.e_cm_sq });
        self.m3_sq = self.m1_fac * self.m1_fac * (if self.mu_p > self.e_cm_sq { self.mu_p } else { self.e_cm_sq });
        self.m4_sq = self.m1_fac * self.m1_fac * self.e_cm_sq;
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

    pub fn deform(&self, mom: &LorentzVector<f64>) -> Result<LorentzVector<Complex>, &str> {
        Err("not implemented")
    }
}
