use arrayvec::ArrayVec;
use float;
use num_traits::{FromPrimitive, One, Zero};
use topologies::Topology;
use vector::LorentzVector;
use Complex;

impl Topology {
    /*    1-LOOP    */
    #[inline]
    fn xij(&self, i: usize, j: usize) -> float {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        float::one()
            + float::from_f64(qji.square()).unwrap()
                / (float::from_f64(2.0 * qji.dot(&pi)).unwrap())
    }

    #[inline]
    fn tij(&self, i: usize, j: usize, x: float) -> float {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        float::one()
            / ((float::one() - x) * float::from_f64(2.0 * pi.dot(&qji)).unwrap()
                + float::from_f64(qji.square()).unwrap())
    }

    fn aij(&self, i: usize, j: usize) -> float {
        let mut result = float::one();
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

    fn bi(&self, i: usize) -> float {
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

    pub fn counterterm(&self, cut_loop_mom: &[LorentzVector<Complex>]) -> Complex {
        match self.n_loops {
            1 => {
                if self.on_shell_flag == 0 {
                    return Complex::default();
                }

                //println!("One Loop CounterTerms");
                let mom = &cut_loop_mom[0];

                //Propagators (MAX 10 )
                //TODO: Make it flexible
                let props: ArrayVec<[num::Complex<float>; 10]> = self.loop_lines[0]
                    .propagators
                    .iter()
                    .map(|x| {
                        let q: LorentzVector<Complex> = x.q.cast().to_complex(true);
                        (mom + q).square() - float::from_f64(x.m_squared).unwrap()
                    })
                    .collect();

                //Use specific CT whenever possible
                let use_special_ct = true;
                //Compute CT
                if use_special_ct && props.len() == 4 {
                    // TODO: convert qs to float
                    let s11 = (self.loop_lines[0].propagators[3].q
                        - self.loop_lines[0].propagators[0].q)
                        .cast::<float>()
                        .square();
                    let s22 = (self.loop_lines[0].propagators[0].q
                        - self.loop_lines[0].propagators[1].q)
                        .cast::<float>()
                        .square();
                    let s33 = (self.loop_lines[0].propagators[1].q
                        - self.loop_lines[0].propagators[2].q)
                        .cast::<float>()
                        .square();
                    let s44 = (self.loop_lines[0].propagators[2].q
                        - self.loop_lines[0].propagators[3].q)
                        .cast::<float>()
                        .square();
                    let s12 = (self.loop_lines[0].propagators[1].q
                        - self.loop_lines[0].propagators[3].q)
                        .cast::<float>()
                        .square();
                    let s23 = (self.loop_lines[0].propagators[0].q
                        - self.loop_lines[0].propagators[2].q)
                        .cast::<float>()
                        .square();
                    //println!("{:?}::{:?}::{:?}::{:?}", s44, s33, s22, s11);
                    let s34 = s12;
                    let s14 = s23;
                    let ai = vec![
                        (float::one() - (s44 + s33) / s12) / (s23 - (s11 * s33 + s22 * s44) / s12),
                        (float::one() - (s11 + s44) / s23) / (s12 - (s22 * s44 + s33 * s11) / s23),
                        (float::one() - (s22 + s11) / s12) / (s23 - (s33 * s11 + s44 * s22) / s12),
                        (float::one() - (s33 + s22) / s23) / (s12 - (s44 * s22 + s11 * s33) / s23),
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
                                let mut tri_prod = Complex::new(float::one(), float::zero());
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
                                        float::zero()
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
