use arrayvec::ArrayVec;
use float;
use num::Complex;
use num_traits::{Float, FromPrimitive, One, Zero};
use topologies::Topology;
use vector::{LorentzVector, RealNumberLike};

impl Topology {
    /*    1-LOOP    */
    #[inline]
    fn xij<
        T: From<float>
            + FromPrimitive
            + Float
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + 'static,
    >(
        &self,
        i: usize,
        j: usize,
    ) -> T {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        T::one() + T::from_f64(qji.square()).unwrap() / (T::from_f64(2.0 * qji.dot(&pi)).unwrap())
    }

    #[inline]
    fn tij<
        T: From<float>
            + FromPrimitive
            + Float
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + 'static,
    >(
        &self,
        i: usize,
        j: usize,
        x: T,
    ) -> T {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        T::one()
            / ((T::one() - x) * T::from_f64(2.0 * pi.dot(&qji)).unwrap()
                + T::from_f64(qji.square()).unwrap())
    }

    fn aij<
        T: From<float>
            + FromPrimitive
            + Float
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + 'static,
    >(
        &self,
        i: usize,
        j: usize,
    ) -> T {
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

    fn bi<
        T: From<float>
            + FromPrimitive
            + Float
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + 'static,
    >(
        &self,
        i: usize,
    ) -> T {
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

    pub fn counterterm<
        T: From<float>
            + FromPrimitive
            + Float
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + 'static,
    >(
        &self,
        cut_loop_mom: &[LorentzVector<Complex<T>>],
    ) -> Complex<T> {
        match self.n_loops {
            1 => {
                if self.on_shell_flag == 0 {
                    return Complex::default();
                }

                //println!("One Loop CounterTerms");
                let mom = &cut_loop_mom[0];

                //Propagators (MAX 10 )
                //TODO: Make it flexible
                let props: ArrayVec<[num::Complex<T>; 10]> = self.loop_lines[0]
                    .propagators
                    .iter()
                    .map(|x| {
                        let q: LorentzVector<Complex<T>> = x.q.cast().to_complex(true);
                        (mom + q).square() - T::from_f64(x.m_squared).unwrap()
                    })
                    .collect();

                //Use specific CT whenever possible
                let use_special_ct = false;
                //Compute CT
                if use_special_ct && props.len() == 4 {
                    // TODO: convert qs to float
                    let s11 = (self.loop_lines[0].propagators[3].q
                        - self.loop_lines[0].propagators[0].q)
                        .cast::<T>()
                        .square();
                    let s22 = (self.loop_lines[0].propagators[0].q
                        - self.loop_lines[0].propagators[1].q)
                        .cast::<T>()
                        .square();
                    let s33 = (self.loop_lines[0].propagators[1].q
                        - self.loop_lines[0].propagators[2].q)
                        .cast::<T>()
                        .square();
                    let s44 = (self.loop_lines[0].propagators[2].q
                        - self.loop_lines[0].propagators[3].q)
                        .cast::<T>()
                        .square();
                    let s12 = (self.loop_lines[0].propagators[1].q
                        - self.loop_lines[0].propagators[3].q)
                        .cast::<T>()
                        .square();
                    let s23 = (self.loop_lines[0].propagators[0].q
                        - self.loop_lines[0].propagators[2].q)
                        .cast::<T>()
                        .square();

                    let ai = vec![
                        (T::one() - (s44 + s33) / s12) / (s23 - (s11 * s33 + s22 * s44) / s12),
                        (T::one() - (s11 + s44) / s23) / (s12 - (s22 * s44 + s33 * s11) / s23),
                        (T::one() - (s22 + s11) / s12) / (s23 - (s33 * s11 + s44 * s22) / s12),
                        (T::one() - (s33 + s22) / s23) / (s12 - (s44 * s22 + s11 * s33) / s23),
                    ];
                    let ct_numerator =
                        -props[0] * ai[0] - props[1] * ai[1] - props[2] * ai[2] - props[3] * ai[3];
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
