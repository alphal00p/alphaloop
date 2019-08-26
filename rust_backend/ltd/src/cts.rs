use amplitude::eeAA;
use arrayvec::ArrayVec;
use num::Complex;
use topologies::{LTDCache, Topology};
use utils;
use vector::LorentzVector;
use FloatLike;

impl Topology {
    #[inline]
    ///Compute the xi's for the collinear limits, here mom should always by on-shell
    fn collinear_x<T: FloatLike>(
        &self,
        loopmom: &LorentzVector<Complex<T>>,
        mom: &LorentzVector<f64>,
    ) -> Complex<T> {
        let eta = mom.dual();
        (loopmom[0] * T::from_f64(eta[0]).unwrap()
            - loopmom[1] * T::from_f64(eta[1]).unwrap()
            - loopmom[2] * T::from_f64(eta[2]).unwrap()
            - loopmom[3] * T::from_f64(eta[3]).unwrap())
            / T::from_f64(eta.dot(mom)).unwrap()
    }

    #[inline]
    fn xij<T: FloatLike>(&self, i: usize, j: usize) -> T {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        T::one() + T::from_f64(qji.square()).unwrap() / (T::from_f64(2.0 * qji.dot(&pi)).unwrap())
    }

    #[inline]
    fn tij<T: FloatLike>(&self, i: usize, j: usize, x: T) -> T {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        T::one()
            / ((T::one() - x) * T::from_f64(2.0 * pi.dot(&qji)).unwrap()
                + T::from_f64(qji.square()).unwrap())
    }

    fn aij<T: FloatLike>(&self, i: usize, j: usize) -> T {
        let mut result = T::one();
        let n_prop = self.loop_lines[0].propagators.len();
        let im1 = if i == 0 { n_prop - 1 } else { i - 1 };
        let x0 = self.xij(i, j);

        for k in 0..n_prop {
            if k == i || k == im1 || k == j {
                continue;
            } else {
                result *= self.tij(i, k, x0);
            }
        }
        result
    }

    fn bi<T: FloatLike>(&self, i: usize) -> T {
        let pi = self.external_kinematics[i];
        let n_prop = self.loop_lines[0].propagators.len();
        let tau = 1e-10 * self.e_cm_squared;
        let ip1 = (i + 1) % n_prop;
        let im1 = if i == 0 { n_prop - 1 } else { i - 1 };

        if pi.square().abs() < tau {
            self.aij(i, ip1)
        } else {
            self.aij(ip1, im1)
        }
    }
    fn eval_porp_limit<T: FloatLike>(
        &self,
        eval_n: usize,
        coll_n: usize,
        x: &Complex<T>,
        p: &LorentzVector<f64>,
    ) -> LorentzVector<Complex<T>> {
        let qij =
            self.loop_lines[0].propagators[eval_n].q - self.loop_lines[0].propagators[coll_n].q;
        LorentzVector {
            t: x * T::from_f64(p.t).unwrap() + T::from_f64(qij.t).unwrap(),
            x: x * T::from_f64(p.x).unwrap() + T::from_f64(qij.x).unwrap(),
            y: x * T::from_f64(p.y).unwrap() + T::from_f64(qij.y).unwrap(),
            z: x * T::from_f64(p.z).unwrap() + T::from_f64(qij.z).unwrap(),
        }
    }

    pub fn counterterm<T: FloatLike>(
        &self,
        k_def: &[LorentzVector<Complex<T>>],
        cut_2energy: Complex<T>,
        cut_id: usize,
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        //Define mu_uv_sq for the collinear counterterms
        let mu_uv_sq = Complex::new(
            T::from_f64(self.settings.general.mu_uv_sq_re_im[0]).unwrap(),
            T::from_f64(self.settings.general.mu_uv_sq_re_im[1]).unwrap(),
        );
        match self.n_loops {
            1 => {
                if self.on_shell_flag == 0 {
                    return Complex::default();
                }

                let ll = &self.loop_lines[0];

                //Propagators (MAX 10 )
                //TODO: Make it flexible
                let props: ArrayVec<[num::Complex<T>; 10]> = ll
                    .propagators
                    .iter()
                    .map(|p| {
                        utils::powi(k_def[0].t + T::from_f64(p.q.t).unwrap(), 2)
                            - cache.complex_prop_spatial[p.id]
                    })
                    .collect();

                let kqs: ArrayVec<[LorentzVector<Complex<T>>; 10]> = ll
                    .propagators
                    .iter()
                    .map(|p| LorentzVector {
                        t: k_def[0].t + T::from_f64(p.q.t).unwrap(),
                        x: k_def[0].x + T::from_f64(p.q.x).unwrap(),
                        y: k_def[0].y + T::from_f64(p.q.y).unwrap(),
                        z: k_def[0].z + T::from_f64(p.q.z).unwrap(),
                    })
                    .collect();

                //Compute CT
                //Use specific CT whenever possible
                if self.settings.general.use_collinear_ct
                    && props.len() == 4
                    && self.on_shell_flag == 8
                {
                    let mut ct_numerator: Complex<T> = Complex::default();
                    //Incoming external: p[n]=q[n]-q[n-1]
                    let x4 =
                        self.collinear_x(&kqs[3], &(ll.propagators[3].q - ll.propagators[2].q)); //)self.external_kinematics[3]);

                    //let ct_uv = T::one();
                    let ct_uv = mu_uv_sq
                        * mu_uv_sq
                        * utils::finv(mu_uv_sq - props[2])
                        * utils::finv(mu_uv_sq - props[3]);

                    // Use explicit expression or limits to compute the coutnerterms
                    ct_numerator -= if false {
                        let s = (kqs[1] - kqs[3]).square();
                        let t = (kqs[0] - kqs[2]).square();

                        let m1 = T::from_f64((self.external_kinematics[0]).square()).unwrap();
                        let m3 = T::from_f64((self.external_kinematics[2]).square()).unwrap();

                        let x4b: Complex<T> = Complex::new(T::one(), T::zero()) - x4;

                        utils::finv((x4 * t + x4b * m1) * (x4b * s + x4 * m3))
                            * (props[0] * props[1])
                            * ct_uv
                    } else {
                        let prop0_c4 = self
                            .eval_porp_limit(
                                0,
                                3,
                                &x4,
                                &(ll.propagators[3].q - ll.propagators[2].q),
                            )
                            .square();
                        let prop1_c4 = self
                            .eval_porp_limit(
                                1,
                                3,
                                &x4,
                                &(ll.propagators[3].q - ll.propagators[2].q),
                            )
                            .square();
                        utils::finv(prop0_c4 * prop1_c4) * (props[0] * props[1]) * ct_uv
                    };
                    ct_numerator
                } else {
                    /*
                    Generic one-loop CounterTerms
                    */

                    let mut ct_numerator = Complex::default();

                    //Loop Over all the aij
                    for i in 0..props.len() {
                        for j in 0..props.len() {
                            let ip1 = (i + 1) % props.len();
                            let im1 = if i == 0 { props.len() - 1 } else { i - 1 };
                            let im2 = if i <= 1 { props.len() - 2 + i } else { i - 2 };

                            if j == i || j == im1 {
                                continue;
                            } else {
                                //Define triangle propagators
                                let mut tri_prod = Complex::new(T::one(), T::zero());
                                for k in 0..props.len() {
                                    if k == i || k == im1 || k == j {
                                        continue;
                                    } else {
                                        tri_prod *= props[k];
                                    }
                                }
                                //Get the right factor
                                ct_numerator -= tri_prod
                                    * if j == ip1 {
                                        self.bi(i)
                                    } else if j == im2 {
                                        //Avoid double counting
                                        T::zero()
                                    } else {
                                        self.aij(i, j)
                                    };
                            }
                        }
                    }
                    ct_numerator
                }
            }
            _ => Complex::default(),
        }
    }
}
