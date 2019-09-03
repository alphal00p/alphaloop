use arrayvec::ArrayVec;
use gamma_chain::GammaChain;
use num::Complex;
use topologies::{Cut, LTDCache, Topology};
use utils;
use vector::LorentzVector;
use FloatLike;
use MAX_LOOP;

#[allow(
    non_snake_case,
    non_upper_case_globals,
    non_camel_case_types,
    dead_code
)]
pub mod Parameters {
    pub const g_f: f64 = 1.166390e-05;
    pub const alpha_s: f64 = 1.180000e-01;
    pub const alpha_ew: f64 = 1. / 1.325070e+02;
    pub const C_F: f64 = 4. / 3.;
}

#[allow(non_camel_case_types, dead_code)]
enum Polarizations {
    UPlus,
    UMinus,
    UBarPlus,
    UBarMinus,
    VPlus,
    VMinus,
    VBarPlus,
    VBarMinus,
    APlus,
    AMinus,
}

#[derive(Debug)]
struct Diagram<'a> {
    alias: &'static str,
    cuts: &'a [usize],
    ct: bool,
}

impl<'a> Diagram<'a> {
    pub fn new(
        diag_name: &'static str,
        cuts: &'a [usize],
        ct: bool,
    ) -> Result<Diagram<'a>, &'static str> {
        Ok(Diagram {
            alias: diag_name,
            cuts: cuts,
            ct: ct,
        })
    }
}

impl Topology {
    pub fn evaluate_amplitude_cut<T: FloatLike>(
        &self,
        mut k_def: &mut ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &'static str> {
        match self.n_loops {
            1 => {
                match &self.name {
                    v if v.contains("eeAA_amplitude") => {
                        self.set_loop_momentum_energies(k_def, cut, mat, cache);
                        let ll = &self.loop_lines[0];
                        // compute propagators
                        let props: ArrayVec<[num::Complex<T>; 10]> = ll
                            .propagators
                            .iter()
                            .map(|p| {
                                utils::powi(k_def[0].t + T::from_f64(p.q.t).unwrap(), 2)
                                    - cache.complex_prop_spatial[p.id]
                            })
                            .collect();
                        // compute residue energy
                        let mut cut_2energy = Complex::new(T::one(), T::zero());
                        let mut cut_id = 0;
                        for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
                            // get the loop line result from the cache
                            cut_2energy *= match ll_cut {
                                Cut::PositiveCut(j) | Cut::NegativeCut(j) => {
                                    cut_id = ll.propagators[*j].id;
                                    cache.complex_cut_energies[cut_id] * Into::<T>::into(2.)
                                }
                                _ => Complex::new(T::one(), T::zero()),
                            };
                        }
                        //Rotate back PS for numerator
                        let rot_matrix = self.rotation_matrix;
                        let mut external_kinematics = self.external_kinematics.clone();

                        for e in &mut external_kinematics {
                            let old_x = e.x;
                            let old_y = e.y;
                            let old_z = e.z;
                            e.x = rot_matrix[0][0] * old_x
                                + rot_matrix[1][0] * old_y
                                + rot_matrix[2][0] * old_z;
                            e.y = rot_matrix[0][1] * old_x
                                + rot_matrix[1][1] * old_y
                                + rot_matrix[2][1] * old_z;
                            e.z = rot_matrix[0][2] * old_x
                                + rot_matrix[1][2] * old_y
                                + rot_matrix[2][2] * old_z;
                        }
                        let mut l = k_def[0];
                        let old_x = l.x;
                        let old_y = l.y;
                        let old_z = l.z;
                        l.x = old_x * T::from_f64(rot_matrix[0][0]).unwrap()
                            + old_y * T::from_f64(rot_matrix[1][0]).unwrap()
                            + old_z * T::from_f64(rot_matrix[2][0]).unwrap();
                        l.y = old_x * T::from_f64(rot_matrix[0][1]).unwrap()
                            + old_y * T::from_f64(rot_matrix[1][1]).unwrap()
                            + old_z * T::from_f64(rot_matrix[2][1]).unwrap();
                        l.z = old_x * T::from_f64(rot_matrix[0][2]).unwrap()
                            + old_y * T::from_f64(rot_matrix[1][2]).unwrap()
                            + old_z * T::from_f64(rot_matrix[2][2]).unwrap();
                        //Create Amplitude
                        let mut amp =
                            eeAA::new(external_kinematics, props.to_vec(), l, cut_2energy, cut_id)
                                .unwrap();
                        return amp.compute_amplitude();
                    }
                    _ => {
                        let v = self.evaluate_cut(&mut k_def, cut, mat, cache)?;
                        // Assuming that there is no need for the residue energy or the cut_id
                        let ct =
                            self.counterterm(&k_def[..self.n_loops], Complex::default(), 0, cache);
                        Ok(v * (ct + T::one()))
                    }
                }
            }
            _ => {
                let v = self.evaluate_cut(&mut k_def, cut, mat, cache)?;
                // Assuming that there is no need for the residue energy or the cut_id
                let ct = self.counterterm(&k_def[..self.n_loops], Complex::default(), 0, cache);
                Ok(v * (ct + T::one()))
            }
        }
    }
}

#[allow(unused_variables)]
fn compute_polarization<T: FloatLike>(
    p: LorentzVector<f64>,
    polarization: Polarizations,
) -> Result<ArrayVec<[Complex<T>; 4]>, &'static str> {
    //Only for massless case
    //Spherical coordinates of the spatial part of p
    let rho = p.spatial_squared().sqrt();
    let theta = (p[3] / rho).acos();
    let phi = if p[2].abs() < 1e-10 && p[1].abs() < 1e-10 && p[3] > 0.0 {
        0.0
    } else if p[2].abs() < 1e-10 && p[1].abs() < 1e-10 && p[3] < 0.0 {
        std::f64::consts::PI
    } else {
        //p.arctan2(p[2], p[1])
        p[2].atan2(p[1])
    };
    let phase = (Complex::new(0.0, 1.0) * phi).exp();

    // define gamma matrices in Weyl representation,
    let gamma0 = [
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
    ];
    let gamma5 = [
        [-1., 0., 0., 0.],
        [0., -1., 0., 0.],
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
    ];

    //HELAS convention for V_+- to U_-+ conversion
    // E > 0 : v_plus = -u_minus, v_minus = -u_plus
    // E < 0 : v_plus = u_minus, v_minus = u_plus
    let factor = match polarization {
        Polarizations::VPlus
        | Polarizations::VBarPlus
        | Polarizations::VMinus
        | Polarizations::VBarMinus
            if p[0] > 0. =>
        {
            -1.0
        }
        _ => 1.0,
    };

    let pol = match polarization {
        Polarizations::UPlus
        | Polarizations::UBarPlus
        | Polarizations::VMinus
        | Polarizations::VBarMinus => {
            // spinors as documented in HELAS reference, rewritten in Dirac representation and in polar coordinates
            let mut u_plus: Vec<Complex<f64>> = if p[0] > 0.0 {
                vec![
                    Complex::default(),
                    Complex::default(),
                    Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                    (1.0 - (theta).cos()).sqrt() * phase,
                ]
            } else {
                vec![
                    Complex::new(0.0, (1.0 + (theta).cos()).sqrt()),
                    Complex::new(0.0, 1.0) * (1.0 - (theta).cos()).sqrt() * phase,
                    Complex::default(),
                    Complex::default(),
                ]
            }
            .iter()
            .map(|x| factor * x * (rho).sqrt())
            .collect();
            // this line adjust sign to HELAS convention
            if p[2] == 0.0 && p[1] == 0.0 && p[3] < 0.0 {
                u_plus = u_plus.iter().map(|x| -x).collect();
            }
            match polarization {
                Polarizations::UPlus | Polarizations::VMinus => u_plus,
                _ => {
                    let mut u_bar_plus = Vec::new();
                    for g in gamma0.iter() {
                        let mut res = Complex::default();
                        for (i, x) in u_plus.iter().enumerate() {
                            res += x.conj() * g[i];
                        }
                        u_bar_plus.push(res);
                    }
                    u_bar_plus
                }
            }
        }
        Polarizations::UMinus
        | Polarizations::UBarMinus
        | Polarizations::VPlus
        | Polarizations::VBarPlus => {
            // spinors as documented in HELAS reference, rewritten in Dirac representation and in polar coordinates
            // according to HELAS convention, v_plus = -u_minus, v_minus = -u_plus
            let mut u_minus: Vec<Complex<f64>> = if p[0] > 0.0 {
                vec![
                    -(1.0 - (theta).cos()).sqrt() / phase,
                    Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                    Complex::default(),
                    Complex::default(),
                ]
            } else {
                vec![
                    Complex::default(),
                    Complex::default(),
                    -Complex::new(0.0, 1.0) * (1.0 - (theta).cos()).sqrt() / phase,
                    Complex::new(0.0, (1.0 + (theta).cos()).sqrt()),
                ]
            }
            .iter()
            .map(|x| factor * x * (rho).sqrt())
            .collect();
            // this line adjust sign to HELAS convention
            if p[2] == 0.0 && p[1] == 0.0 && p[3] < 0.0 {
                u_minus = u_minus.iter().map(|x| -x).collect();
            }
            match polarization {
                Polarizations::UMinus | Polarizations::VPlus => u_minus,
                _ => {
                    let mut u_bar_minus = Vec::new();
                    for g in gamma0.iter() {
                        let mut res = Complex::default();
                        for (i, x) in u_minus.iter().enumerate() {
                            res += x.conj() * g[i];
                        }
                        u_bar_minus.push(res);
                    }
                    u_bar_minus
                }
            }
        }
        Polarizations::APlus => {
            // photon polarization vectors
            let a1 = vec![
                0.,
                (theta).cos() * (phi).cos(),
                (theta).cos() * (phi).sin(),
                -(theta).sin(),
            ];
            let a2 = vec![0., -(phi).sin(), (phi).cos(), 0.];

            //helicity eigenstates of photons
            let ii = Complex::new(0.0, 1.0);
            let norm = 1. / (2.0 as f64).sqrt();
            let mut a_plus = Vec::new();
            for i in 0..a1.len() {
                a_plus.push(norm * (-a1[i] - ii * a2[i]));
            }

            a_plus
        }
        Polarizations::AMinus => {
            // photon polarization vectors
            let a1 = vec![
                0.,
                (theta).cos() * (phi).cos(),
                (theta).cos() * (phi).sin(),
                -(theta).sin(),
            ];
            let a2 = vec![0., -(phi).sin(), (phi).cos(), 0.];
            //helicity eigenstates of photons
            let ii = Complex::new(0.0, 1.0);
            let norm = 1. / (2.0 as f64).sqrt();
            let mut a_minus = Vec::new();
            for i in 0..a1.len() {
                a_minus.push(norm * (a1[i] - ii * a2[i]));
            }

            a_minus
        }
    };
    Ok(pol
        .iter()
        .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
        .collect())
}
//END compute_polarization
#[derive(Debug)]
#[allow(unused_variables, non_camel_case_types)]
pub struct eeAA<'a, T: FloatLike> {
    external_kinematics: Vec<LorentzVector<Complex<T>>>,
    den: Vec<Complex<T>>,
    vectors: [LorentzVector<Complex<T>>; 8],
    vbar: ArrayVec<[Complex<T>; 4]>,
    u: ArrayVec<[Complex<T>; 4]>,
    a_mu: LorentzVector<Complex<T>>,
    a_nu: LorentzVector<Complex<T>>,
    cut_2energy: Complex<T>,
    diags_and_cuts: [Diagram<'a>; 14],
    cut_id: usize,
}

#[allow(unused_variables, non_camel_case_types)]
impl<'a, T: FloatLike> eeAA<'a, T> {
    pub fn new(
        external_kinematics: Vec<LorentzVector<f64>>,
        propagators: Vec<Complex<T>>,
        loop_momentum: LorentzVector<Complex<T>>,
        cut_2energy: Complex<T>,
        cut_id: usize,
    ) -> Result<eeAA<'a, T>, &'static str> {
        //For on-shell photons
        let a_mu = LorentzVector::from_slice(
            &compute_polarization(external_kinematics[3], Polarizations::APlus).unwrap(),
        );
        let a_nu = LorentzVector::from_slice(
            &compute_polarization(external_kinematics[4], Polarizations::APlus).unwrap(),
        );
        //Spinors
        //This agrees with HELAS convention
        let u = if external_kinematics[0][0] > 0. {
            compute_polarization(external_kinematics[0], Polarizations::UPlus).unwrap()
        } else {
            compute_polarization(-external_kinematics[0], Polarizations::UPlus).unwrap()
        };
        let vbar = if external_kinematics[2][0] > 0. {
            compute_polarization(external_kinematics[2], Polarizations::UBarPlus).unwrap()
        } else {
            compute_polarization(external_kinematics[2], Polarizations::VBarPlus).unwrap()
        };

        let gamma_0 = LorentzVector::from_args(
            Complex::new(T::one(), T::zero()),
            Complex::new(T::zero(), T::zero()),
            Complex::new(T::zero(), T::zero()),
            Complex::new(T::zero(), T::zero()),
        );

        let k = loop_momentum;
        let mut ps: Vec<LorentzVector<Complex<T>>> = external_kinematics
            .iter()
            .map(|mom| (*mom).map(|x| Complex::new(T::from_f64(x).unwrap(), T::zero())))
            .collect();

        let mut den = propagators;
        den[cut_id] = cut_2energy;
        // When considering also the UV propagators
        ps.remove(1);
        let vectors = [
            -k - ps[0],
            -k - ps[0] - ps[1],
            -k + ps[3],
            -k,
            -ps[1] - ps[2],
            a_mu,
            a_nu,
            gamma_0,
        ];
        let diags_and_cuts = [
            //Diagrams
            Diagram::new("D1", &[0, 3, 4], false).unwrap(),
            Diagram::new("D2", &[0, 2, 3], false).unwrap(),
            Diagram::new("D3", &[0, 2, 3, 4], false).unwrap(),
            Diagram::new("D4", &[0, 3], false).unwrap(),
            //Counterterms
            Diagram::new("IR", &[0, 2, 4], true).unwrap(),
            Diagram::new("UV1", &[1], true).unwrap(),
            Diagram::new("UV2", &[1], true).unwrap(),
            Diagram::new("UV4", &[1], true).unwrap(),
            Diagram::new("UVIR", &[1], true).unwrap(),
            //UV approximations of relevant contributions
            Diagram::new("D1UV_LO", &[0], false).unwrap(),
            Diagram::new("D2UV_LO", &[0], false).unwrap(),
            Diagram::new("D3UV_LO", &[0], false).unwrap(),
            Diagram::new("D4UV_LO", &[0], false).unwrap(),
            Diagram::new("IRUV_LO", &[0], false).unwrap(),
        ];

        return Ok(eeAA {
            external_kinematics: ps,
            den: den,
            vectors: vectors,
            u: u,
            vbar: vbar,
            a_mu: a_mu,
            a_nu: a_nu,
            diags_and_cuts: diags_and_cuts,
            cut_2energy: cut_2energy,
            cut_id: cut_id,
        });
    }

    pub fn compute_chain(&self, indices: &[i8]) -> Result<Complex<T>, &'static str> {
        GammaChain::new(
            self.vbar.as_slice(),
            self.u.as_slice(),
            indices,
            &self.vectors,
        )
        .unwrap()
        .compute_chain()
    }

    pub fn compute_amplitude(&mut self) -> Result<Complex<T>, &'static str> {
        let mut res: Complex<T> = Complex::default();
        let diaglist = if self.cut_2energy.norm() > T::from_f64(1e4).unwrap() {
            [
                "D1UV_LO", "D2UV_LO", "D3UV_LO", "D4UV_LO", "IRUV_LO", "UV1", "UV2", "UV4", "UVIR",
            ]
        } else {
            ["D1", "D2", "D3", "D4", "IR", "UV1", "UV2", "UV4", "UVIR"]
        };
        //Use diagram own denominator
        for diag_and_cut in self.diags_and_cuts.iter() {
            res += if diag_and_cut.cuts.iter().any(|v| v == &self.cut_id)
                & diaglist.iter().any(|v| v == &diag_and_cut.alias)
            {
                //Compute denominator
                let mut diag_den = Complex::new(T::one(), T::zero());
                for cut in diag_and_cut.cuts.iter() {
                    diag_den *= self.den[*cut];
                }
                //Check if it is a counter term and return result
                if diag_and_cut.ct {
                    -self.numerator(diag_and_cut.alias).unwrap() * utils::finv(diag_den)
                } else {
                    self.numerator(diag_and_cut.alias).unwrap() * utils::finv(diag_den)
                }
            } else {
                Complex::default()
            };
        }

        return Ok(res);
    }

    pub fn numerator(&self, name: &'static str) -> Result<Complex<T>, &'static str> {
        //invariants
        let s23_inv =
            utils::finv((self.external_kinematics[1] + self.external_kinematics[2]).square());
        let factor = Complex::new(
            T::from_f64(
                Parameters::C_F * Parameters::alpha_s * Parameters::alpha_ew / 9.0
                    * (4.0 * std::f64::consts::PI).powf(2.0),
            )
            .unwrap(),
            T::zero(),
        );

        let numerator = match name {
            "D1" => (-factor * s23_inv) * self.compute_chain(&[6, 5, -1, 3, 7, 4, -1]).unwrap(),
            "D2" => (-factor * s23_inv) * self.compute_chain(&[-1, 2, 6, 3, -1, 5, 7]).unwrap(),
            "D3" => (-factor) * self.compute_chain(&[-1, 2, 6, 3, 7, 4, -1]).unwrap(),
            "D4" => {
                (-factor * s23_inv * s23_inv)
                    * self.compute_chain(&[6, 5, -1, 3, -1, 5, 7]).unwrap()
            }
            "IR" => (-factor * s23_inv) * self.compute_chain(&[-1, 2, 6, 5, 7, 4, -1]).unwrap(),
            "UV1" => {
                //First contribution of the double derivative
                let fdd = self.compute_chain(&[6, 5, -1, 1, 7, 1, -1]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                //Second contribution of the double derivative
                let dfd = -(self.compute_chain(&[6, 5, -1, 8, 7, 1, -1]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 1, 7, 8, -1]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                //Third contribution of the double derivative
                let ddf = self.compute_chain(&[6, 5, -1, 8, 7, 8, -1]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));

                (-factor * s23_inv / T::from_f64(2.0).unwrap()) * (fdd + dfd + ddf)
            }
            "UV2" => {
                let fdd = self.compute_chain(&[-1, 1, 6, 1, -1, 5, 7]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[-1, 8, 6, 1, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[-1, 1, 6, 8, -1, 5, 7]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = self.compute_chain(&[-1, 8, 6, 8, -1, 5, 7]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));
                (-factor * s23_inv / T::from_f64(2.0).unwrap()) * (fdd + dfd + ddf)
            }
            "UV4" => {
                //LEADING DIVERGENCE
                let fd = self.compute_chain(&[6, 5, -1, 1, -1, 5, 7]).unwrap()
                    * T::from_f64(-2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));

                let df = -self.compute_chain(&[6, 5, -1, 8, -1, 5, 7]).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 1));

                //SUBLEADING DIVERGENCE
                let fdd = self.compute_chain(&[6, 5, -1, 1, 5, 1, -1, 5, 7]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[6, 5, -1, 8, 5, 1, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 1, 5, 8, -1, 5, 7]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = self.compute_chain(&[6, 5, -1, 8, 5, 8, -1, 5, 7]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));

                (-factor * s23_inv * s23_inv)
                    * ((fd + df) - (fdd + dfd + ddf) / T::from_f64(2.0).unwrap())
            }
            "UVIR" => {
                let fdd = self.compute_chain(&[-1, 1, 6, 5, 7, 1, -1]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[-1, 8, 6, 5, 7, 1, -1]).unwrap()
                    + self.compute_chain(&[-1, 1, 6, 5, 7, 8, -1]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = self.compute_chain(&[-1, 8, 6, 5, 7, 8, -1]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));
                (factor * s23_inv / T::from_f64(2.0).unwrap()) * (fdd + dfd + ddf)
            }
            "D1UV_LO" => {
                let fdd = self.compute_chain(&[6, 5, -1, 3, 7, 4, -1]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[6, 5, -1, 8, 7, 4, -1]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 3, 7, 8, -1]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = self.compute_chain(&[6, 5, -1, 8, 7, 8, -1]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));
                (-factor * s23_inv / T::from_f64(2.0).unwrap()) * (fdd + dfd + ddf)
            }
            "D2UV_LO" => {
                let fdd = self.compute_chain(&[-1, 2, 6, 3, -1, 5, 7]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[-1, 8, 6, 3, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[-1, 2, 6, 8, -1, 5, 7]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = self.compute_chain(&[-1, 8, 6, 8, -1, 5, 7]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));
                (-factor * s23_inv / T::from_f64(2.0).unwrap()) * (fdd + dfd + ddf)
            }
            "D3UV_LO" => {
                let fddd = self.compute_chain(&[-1, 2, 6, 3, 7, 4, -1]).unwrap()
                    * T::from_f64(-120.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 6));
                let dfdd = -(self.compute_chain(&[-1, 8, 6, 3, 7, 4, -1]).unwrap()
                    + self.compute_chain(&[-1, 2, 6, 8, 7, 4, -1]).unwrap()
                    + self.compute_chain(&[-1, 2, 6, 3, 7, 8, -1]).unwrap())
                    * T::from_f64(60.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 5));
                let ddfd = (self.compute_chain(&[-1, 8, 6, 8, 7, 4, -1]).unwrap()
                    + self.compute_chain(&[-1, 8, 6, 3, 7, 8, -1]).unwrap()
                    + self.compute_chain(&[-1, 2, 6, 8, 7, 8, -1]).unwrap())
                    * T::from_f64(-24.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dddf = -self.compute_chain(&[-1, 8, 6, 8, 7, 8, -1]).unwrap()
                    * T::from_f64(6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                (-factor / T::from_f64(6.0).unwrap()) * (fddd + dfdd + ddfd + dddf)
            }
            "D4UV_LO" => {
                //LEADING DIVERGENCE
                let fd = self.compute_chain(&[6, 5, -1, 3, -1, 5, 7]).unwrap()
                    * T::from_f64(-2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));

                let df = -self.compute_chain(&[6, 5, -1, 8, -1, 5, 7]).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 1));
                //SUBLEADING DIVERGENCE
                let fdd = (self.compute_chain(&[6, 5, -1, 3, 1, 5, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 3, 5, 1, -1, 5, 7]).unwrap())
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[6, 5, -1, 8, 1, 5, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 8, 5, 1, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 3, 8, 5, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 3, 5, 8, -1, 5, 7]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = (self.compute_chain(&[6, 5, -1, 8, 5, 8, -1, 5, 7]).unwrap()
                    + self.compute_chain(&[6, 5, -1, 8, 8, 5, -1, 5, 7]).unwrap())
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));
                (-factor * s23_inv * s23_inv)
                    * ((fd + df) - (fdd + dfd + ddf) / T::from_f64(2.0).unwrap())
            }
            "IRUV_LO" => {
                let fdd = self.compute_chain(&[-1, 2, 6, 5, 7, 4, -1]).unwrap()
                    * T::from_f64(12.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 4));
                let dfd = -(self.compute_chain(&[-1, 8, 6, 5, 7, 4, -1]).unwrap()
                    + self.compute_chain(&[-1, 2, 6, 5, 7, 8, -1]).unwrap())
                    * T::from_f64(-6.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 3));
                let ddf = self.compute_chain(&[-1, 8, 6, 5, 7, 8, -1]).unwrap()
                    * T::from_f64(2.).unwrap()
                    * utils::finv(utils::powi(self.cut_2energy, 2));
                (factor * s23_inv / T::from_f64(2.0).unwrap()) * (fdd + dfd + ddf)
            }
            _ => return panic!("Unknown diagram/counterterm {} for eeAA", name),
        };
        return Ok(numerator);
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use rand::Rng;
    use topologies::Topology;
    use Settings;

    #[allow(non_snake_case, dead_code)]
    #[test]
    fn test_ee_aa() {
        let my_top = "manual_eeAA_amplitude_E";
        let settings = Settings::from_file("../../LTD/hyperparameters.yaml");
        let mut topologies = Topology::from_file("../../LTD/topologies.yaml", &settings);
        let topo = topologies.get_mut(my_top).expect("Unknown topology");
        topo.process();

        //Tree Level
        let mut ps = topo.external_kinematics.clone();
        ps.remove(1);

        //For on-shell photons
        let a_mu =
            LorentzVector::from_slice(&compute_polarization(ps[2], Polarizations::APlus).unwrap());
        let a_nu =
            LorentzVector::from_slice(&compute_polarization(ps[3], Polarizations::APlus).unwrap());

        //Spinors
        //This agrees with HELAS convention
        let u = compute_polarization(ps[0], Polarizations::UPlus).unwrap();
        let vbar = compute_polarization(ps[1], Polarizations::VBarPlus).unwrap();
        //let vbar = compute_polarization(-ps[1], Polarizations::UBarPlus).unwrap();
        let vbar0 = compute_polarization::<f64>(ps[1], Polarizations::VBarPlus).unwrap();
        println!("Diff: [");
        for (x0, x) in vbar0.iter().zip(vbar.iter()) {
            println!("{:.5e}", x0 + x);
        }
        println!("]");
        //let vbar = compute_polarization(ps[1], Polarizations::UBarPlus).unwrap();
        //let box_den = props[0] * props[1] * props[2] * props[3] * props[4];
        let gamma_0 = LorentzVector::from_args(
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        );
        let vectors = [
            (-ps[1] - ps[2]).map(|x| Complex::new(x, 0.0)),
            a_mu,
            a_nu,
            gamma_0,
        ];

        let mut num = GammaChain::new(&vbar, &u, &[2, 1, 3], &vectors).unwrap();
        //let mut num = GammaChain::new(vbar.as_slice(), u.as_slice(), &[2, 1, 3], &vectors).unwrap();
        let factor = -Parameters::alpha_ew * 4.0 * std::f64::consts::PI / 9.0;
        let result = num.compute_chain().unwrap() * factor / (ps[1] + ps[2]).square();

        //println!("chain: {:.6e}", num.compute_chain().unwrap());
        //println!("vbar: {:.6e}", LorentzVector::from_slice(&vbar));
        //println!("u: {:.6e}", LorentzVector::from_slice(&u));
        println!("u: [");
        for v in u.iter() {
            println!("{:.5e}", v);
        }
        println!("]");

        println!("vbar: [");
        for v in vbar.iter() {
            println!("{:.5e}", v);
        }
        println!("]");

        println!(
            "a_mu = [{:.5e}, {:.5e}, {:.5e}, {:.5e}]",
            a_mu.t, a_mu.x, a_mu.y, a_mu.z,
        );
        println!(
            "a_nu = [{:.5e}, {:.5e}, {:.5e}, {:.5e}]",
            a_nu.t, a_nu.x, a_nu.y, a_nu.z,
        );

        let expected = Complex::new(0.0, 1.49019647202147543e-2);
        println!("res: {:.5e}", result);
        println!("exp: {:.5e}", expected);

        let vectors = [
            ps[0].map(|x| Complex::new(x, 0.0)),
            ps[1].map(|x| Complex::new(x, 0.0)),
            ps[2].map(|x| Complex::new(x, 0.0)),
            ps[3].map(|x| Complex::new(x, 0.0)),
            (-ps[1] - ps[2]).map(|x| Complex::new(x, 0.0)),
            a_mu,
            a_nu,
        ];
        let mu_r = 91.188000000000002;
        let M0 = GammaChain::new(&vbar, &u, &[6, 5, 7], &vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        let M1 = GammaChain::new(&vbar, &u, &[-1, -2, 6, 5, 7, -2, -1], &vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        let M2 = GammaChain::new(&vbar, &u, &[-1, 1, 6, 5, 7, 2, -1], &vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        let M3 = GammaChain::new(&vbar, &u, &[6, 4, 1, 4, 7], &vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        let M4 = GammaChain::new(&vbar, &u, &[6, 3, 2, 3, 7], &vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        let M5 = GammaChain::new(&vbar, &u, &[6, 1, 4, 1, 7], &vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        println!("M0 = {:?}", M0);
        println!("M1 = {:?}", M1);
        println!("M2 = {:?}", M2);
        println!("M3 = {:?}", M3);
        println!("M4 = {:?}", M4);
        println!("M5 = {:?}", M5);
        println!("M1/M0 = {:?}", M1 / M0);
        println!(
            "s = {:?}, t = {:?}",
            (ps[0] + ps[1]).square(),
            (ps[1] + ps[2]).square(),
        );

        assert!((result - expected).norm() < 1e-10);
    }

    #[allow(non_snake_case, dead_code)]
    #[test]
    fn test_rot() {
        let my_top = "manual_eeAA_amplitude_E";
        let settings = Settings::from_file("../../LTD/hyperparameters.yaml");
        let mut topologies = Topology::from_file("../../LTD/topologies.yaml", &settings);
        let topo = topologies.get_mut(my_top).expect("Unknown topology");
        topo.process();
        //Rotate
        let mut topos = vec![topo.clone()];
        let mut rng = rand::thread_rng();
        let angle = rng.gen::<f64>() * 2. * std::f64::consts::PI;
        let mut rv = (rng.gen::<f64>(), rng.gen::<f64>(), rng.gen::<f64>()); // rotation axis
        let inv_norm = 1.0 / (rv.0 * rv.0 + rv.1 * rv.1 + rv.2 * rv.2).sqrt();
        rv = (rv.0 * inv_norm, rv.1 * inv_norm, rv.2 * inv_norm);
        topos.push(topo.rotate(angle, rv));

        for topo in topos.iter() {
            println!("{:?}", topo.rotation_matrix);
            //Tree Level
            let mut ps = topo.external_kinematics.clone();
            ps.remove(1);
            let rot_matrix = topo.rotation_matrix;
            for e in &mut ps {
                let old_x = e.x;
                let old_y = e.y;
                let old_z = e.z;
                e.x =
                    rot_matrix[0][0] * old_x + rot_matrix[1][0] * old_y + rot_matrix[2][0] * old_z;
                e.y =
                    rot_matrix[0][1] * old_x + rot_matrix[1][1] * old_y + rot_matrix[2][1] * old_z;
                e.z =
                    rot_matrix[0][2] * old_x + rot_matrix[1][2] * old_y + rot_matrix[2][2] * old_z;
            }

            println!("{:?}", ps);

            //For on-shell photons
            let a_mu = LorentzVector::from_slice(
                &compute_polarization(ps[2], Polarizations::APlus).unwrap(),
            );
            let a_nu = LorentzVector::from_slice(
                &compute_polarization(ps[3], Polarizations::APlus).unwrap(),
            );

            //Spinors
            //This agrees with HELAS convention
            let u = compute_polarization(ps[0], Polarizations::UPlus).unwrap();
            let vbar = compute_polarization(ps[1], Polarizations::VBarPlus).unwrap();
            println!("u: [");
            for v in u.iter() {
                println!("{:.5e}", v);
            }
            println!("]");

            println!("vbar: [");
            for v in vbar.iter() {
                println!("{:.5e}", v);
            }
            println!("]");

            println!(
                "a_mu = [{:.5e}, {:.5e}, {:.5e}, {:.5e}]",
                a_mu.t, a_mu.x, a_mu.y, a_mu.z,
            );
            println!(
                "a_nu = [{:.5e}, {:.5e}, {:.5e}, {:.5e}]",
                a_nu.t, a_nu.x, a_nu.y, a_nu.z,
            );

            let vectors = [
                ps[0].map(|x| Complex::new(x, 0.0)),
                ps[1].map(|x| Complex::new(x, 0.0)),
                ps[2].map(|x| Complex::new(x, 0.0)),
                ps[3].map(|x| Complex::new(x, 0.0)),
                (-ps[1] - ps[2]).map(|x| Complex::new(x, 0.0)),
                a_mu,
                a_nu,
            ];
            let M0 = GammaChain::new(&vbar, &u, &[6, 5, 7], &vectors)
                .unwrap()
                .compute_chain()
                .unwrap();
            let M1 = GammaChain::new(&vbar, &u, &[-1, -2, 6, 5, 7, -2, -1], &vectors)
                .unwrap()
                .compute_chain()
                .unwrap();
            let M2 = GammaChain::new(&vbar, &u, &[-1, 1, 6, 5, 7, 2, -1], &vectors)
                .unwrap()
                .compute_chain()
                .unwrap();
            let M3 = GammaChain::new(&vbar, &u, &[6, 4, 1, 4, 7], &vectors)
                .unwrap()
                .compute_chain()
                .unwrap();
            let M4 = GammaChain::new(&vbar, &u, &[6, 3, 2, 3, 7], &vectors)
                .unwrap()
                .compute_chain()
                .unwrap();
            let M5 = GammaChain::new(&vbar, &u, &[6, 1, 4, 1, 7], &vectors)
                .unwrap()
                .compute_chain()
                .unwrap();
            println!("M0 = {:?}", M0);
            println!("M1 = {:?}", M1);
            println!("M2 = {:?}", M2);
            println!("M3 = {:?}", M3);
            println!("M4 = {:?}", M4);
            println!("M5 = {:?}", M5);
            println!("M1/M0 = {:?}", M1 / M0);
            println!(
                "s = {:?}, t = {:?}",
                (ps[0] + ps[1]).square(),
                (ps[1] + ps[2]).square(),
            );
        }
        assert!((1.0 - 2.0) < 1e-10);
    }
}
