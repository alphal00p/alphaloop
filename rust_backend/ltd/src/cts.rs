use arrayvec::ArrayVec;
use float;
use num_traits::One;
use num_traits::Zero;
use topologies::Topology;
use vector::LorentzVector;
use Complex;

impl Topology {
    /*    1-LOOP    */
    #[inline]
    fn xij(&self, i: usize, j: usize) -> float {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        1.0 + qji.square() / (2.0 * qji.dot(&pi))
    }

    #[inline]
    fn tij(&self, i: usize, j: usize, x: float) -> float {
        let qji = self.loop_lines[0].propagators[j].q - self.loop_lines[0].propagators[i].q;
        let pi = self.external_kinematics[i];
        1.0 / ((1.0 - x) * 2.0 * pi.dot(&qji) + qji.square())
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
        //TODO: introduce e_cm_sq
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
                    .map(|x| (mom + &x.q).square() - x.m_squared)
                    .collect();

                //Use specific CT whenever possible
                let use_special_ct = true;
                //Compute CT
                if use_special_ct && props.len() == 4 {
                    let s11 = (self.loop_lines[0].propagators[3].q
                        - self.loop_lines[0].propagators[0].q)
                        .square();
                    let s22 = (self.loop_lines[0].propagators[0].q
                        - self.loop_lines[0].propagators[1].q)
                        .square();
                    let s33 = (self.loop_lines[0].propagators[1].q
                        - self.loop_lines[0].propagators[2].q)
                        .square();
                    let s44 = (self.loop_lines[0].propagators[2].q
                        - self.loop_lines[0].propagators[3].q)
                        .square();
                    let s12 = (self.loop_lines[0].propagators[1].q
                        - self.loop_lines[0].propagators[3].q)
                        .square();
                    let s23 = (self.loop_lines[0].propagators[0].q
                        - self.loop_lines[0].propagators[2].q)
                        .square();
                    //println!("{:?}::{:?}::{:?}::{:?}", s44, s33, s22, s11);
                    let s34 = s12;
                    let s14 = s23;
                    let ai = vec![
                        (1.0 - (s44 + s33) / s12) / (s23 - (s11 * s33 + s22 * s44) / s12),
                        (1.0 - (s11 + s44) / s23) / (s12 - (s22 * s44 + s33 * s11) / s23),
                        (1.0 - (s22 + s11) / s12) / (s23 - (s33 * s11 + s44 * s22) / s12),
                        (1.0 - (s33 + s22) / s23) / (s12 - (s44 * s22 + s11 * s33) / s23),
                    ];
                    let ct_numerator =
                        -ai[0] * props[0] - ai[1] * props[1] - ai[2] * props[2] - ai[3] * props[3];
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
                                ct_numerator -= if j == ip1 {
                                    self.bi(i)
                                } else if j == im2 {
                                    //Avoid double counting
                                    0.0
                                } else {
                                    self.aij(i, j)
                                } * tri_prod;
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
