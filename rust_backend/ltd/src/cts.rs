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

    pub fn counterterm<T: FloatLike>(
        &self,
        k_def: &[LorentzVector<Complex<T>>],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
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

                //Use specific CT whenever possible
                let use_collinear_ct = true;
                //Compute CT
                if use_collinear_ct && props.len() == 4 && self.on_shell_flag == 8 {
                    let mut ct_numerator: Complex<T> = Complex::default();

                    let x4 = self.collinear_x(&kqs[2], &self.external_kinematics[3]);
                    let x4b: Complex<T> = Complex::new(T::one(), T::zero()) - x4;

                    let s = (kqs[1] - kqs[3]).square();
                    let t = (kqs[0] - kqs[2]).square();

                    let m1 = T::from_f64((self.external_kinematics[0]).square()).unwrap();
                    let m3 = T::from_f64((self.external_kinematics[2]).square()).unwrap();

                    let mu_sq = Complex::new(T::zero(), T::from_f64(1e9).unwrap());

                    ct_numerator += utils::finv((x4b * t + x4 * m1) * (x4 * s + x4b * m3))
                        * (props[0] * props[1])
                        * mu_sq
                        * (mu_sq - props[2] - props[3])
                        * utils::finv(mu_sq - props[2])
                        * utils::finv(mu_sq - props[3]);

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
            _ => {
                //println!("No CounterTerms");
                Complex::default()
            }
        }
    }
}
