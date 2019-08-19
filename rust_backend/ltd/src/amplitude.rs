use gamma_chain::GammaChain;
use num::Complex;
use utils;
use vector::LorentzVector;
use FloatLike;

#[allow(non_camel_case_types, dead_code)]
enum Parameters {
    g_f,
    alpha_s,
    alpha_ew,
    C_F,
}

fn parameters(parameter: Parameters) -> f64 {
    match parameter {
        Parameters::g_f => 1.166390e-05,
        Parameters::alpha_s => 1.180000e-01,
        Parameters::alpha_ew => 1. / 1.325070e+02,
        Parameters::C_F => 4. / 3.,
    }
}

#[allow(non_camel_case_types, dead_code)]
enum Polarizations {
    UPlus,
    UMinus,
    UBarPlus,
    UBarMinus,
    APlus,
    AMinus,
}

#[allow(unused_variables)]
fn compute_polarization(
    p: LorentzVector<f64>,
    polarization: Polarizations,
) -> Result<Vec<Complex<f64>>, &'static str> {
    //Only for massless case
    //Spherical coordinates of the spatial part of p
    let rho = p.spatial_squared();
    let theta = (p[3] / rho).acos();
    let phi = if p[2] == 0.0 && p[1] == 0.0 && p[3] > 0.0 {
        0.0
    } else if p[2] == 0.0 && p[1] == 0.0 && p[3] < 0.0 {
        std::f64::consts::PI
    } else {
        //p.arctan2(p[2], p[1])
        (p[2] / p[1]).atan()
    };

    let phase = (Complex::new(0.0, 1.0) * phi).exp();
    // define gamma matrices in Dirac representation, only for the use of computing polarization vectors and possibly change of basis between Dirac and Weyl representations
    let gamma0 = [
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., -1., 0.],
        [0., 0., 0., -1.],
    ];
    let gamma5 = [
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
    ];

    match polarization {
        Polarizations::UPlus | Polarizations::UBarPlus => {
            // spinors as documented in HELAS reference, rewritten in Dirac representation and in polar coordinates
            // according to HELAS convention, v_plus = -u_minus, v_minus = -u_plus
            let mut u_plus: Vec<Complex<f64>> = vec![
                Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                (1.0 - (theta).cos()).sqrt() * phase,
                Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                (1.0 - (theta).cos()).sqrt() * phase,
            ]
            .iter()
            .map(|x| x * (rho / 2.).sqrt())
            .collect();
            // this line adjust sign to HELAS convention
            if p[2] == 0.0 && p[1] == 0.0 && p[3] < 0.0 {
                u_plus = u_plus.iter().map(|x| -x).collect();
            }
            match polarization {
                Polarizations::UPlus => {
                    return Ok(u_plus);
                }
                _ => {
                    let mut u_bar_plus = Vec::new();
                    for g in gamma0.iter() {
                        let mut res = Complex::default();
                        for (i, x) in u_plus.iter().enumerate() {
                            res += x * g[i];
                        }
                        u_bar_plus.push(res);
                    }
                    return Ok(u_bar_plus);
                }
            }
        }
        Polarizations::UMinus | Polarizations::UBarMinus => {
            // spinors as documented in HELAS reference, rewritten in Dirac representation and in polar coordinates
            // according to HELAS convention, v_plus = -u_minus, v_minus = -u_plus
            let mut u_minus: Vec<Complex<f64>> = vec![
                -(1.0 - (theta).cos()).sqrt() / phase,
                Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                (1.0 - (theta).cos()).sqrt() / phase,
                Complex::new(-(1.0 + (theta).cos()).sqrt(), 0.0),
            ]
            .iter()
            .map(|x| x * (rho / 2.).sqrt())
            .collect();
            // this line adjust sign to HELAS convention
            if p[2] == 0.0 && p[1] == 0.0 && p[3] < 0.0 {
                u_minus = u_minus.iter().map(|x| -x).collect();
            }
            match polarization {
                Polarizations::UMinus => {
                    return Ok(u_minus);
                }
                _ => {
                    let mut u_bar_minus = Vec::new();
                    for g in gamma0.iter() {
                        let mut res = Complex::default();
                        for (i, x) in u_minus.iter().enumerate() {
                            res += x * g[i];
                        }
                        u_bar_minus.push(res);
                    }
                    return Ok(u_bar_minus);
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

            return Ok(a_plus);
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

            return Ok(a_minus);
        }
    }
}
//END compute_polarization
#[derive(Debug)]
#[allow(unused_variables, non_camel_case_types)]
pub struct eeAA<T: FloatLike> {
    external_kinematics: Vec<LorentzVector<f64>>,
    den: Vec<Complex<T>>,
    vectors: Vec<LorentzVector<Complex<T>>>,
    vbar: Vec<Complex<T>>,
    u: Vec<Complex<T>>,
    a_mu: LorentzVector<Complex<T>>,
    a_nu: LorentzVector<Complex<T>>,
    mu_uv: Complex<T>,
}

#[allow(unused_variables, non_camel_case_types)]
impl<T: FloatLike> eeAA<T> {
    pub fn new(
        external_kinematics: Vec<LorentzVector<f64>>,
        propagators: Vec<Complex<T>>,
        loop_momentum: LorentzVector<Complex<T>>,
        mu_uv: Complex<T>,
    ) -> Result<eeAA<T>, &'static str> {
        //For positive energies
        //This agrees with HELAS convention
        let u = compute_polarization(external_kinematics[0], Polarizations::UPlus)
            .unwrap()
            .iter()
            .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
            .collect();
        let vbar = compute_polarization(external_kinematics[1], Polarizations::UBarPlus)
            .unwrap()
            .iter()
            .map(|x| -Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
            .collect();

        //random contraction for the offshell photons
        let a_mu = LorentzVector::from_args(
            Complex::new(T::from_f64(1.).unwrap(), T::zero()),
            Complex::new(T::from_f64(2.).unwrap(), T::zero()),
            Complex::new(T::from_f64(3.).unwrap(), T::zero()),
            Complex::new(T::from_f64(4.).unwrap(), T::zero()),
        );

        let a_nu = LorentzVector::from_args(
            Complex::new(T::from_f64(2.).unwrap(), T::zero()),
            Complex::new(T::from_f64(3.).unwrap(), T::zero()),
            Complex::new(T::from_f64(5.).unwrap(), T::zero()),
            Complex::new(T::from_f64(7.).unwrap(), T::zero()),
        );

        let k = loop_momentum;
        let ps: Vec<LorentzVector<Complex<T>>> = external_kinematics
            .iter()
            .map(|mom| (*mom).map(|x| Complex::new(T::from_f64(x).unwrap(), T::zero())))
            .collect();
        let vectors = vec![
            k + ps[0],
            -k - ps[0] - ps[1],
            -k + ps[3],
            -k,
            -ps[1] - ps[2],
            a_mu,
            a_nu,
        ];
        return Ok(eeAA {
            external_kinematics: external_kinematics,
            den: propagators,
            vectors: vectors,
            u: u,
            vbar: vbar,
            a_mu: a_mu,
            a_nu: a_nu,
            mu_uv: mu_uv,
        });
    }

    pub fn compute_chain(&mut self, indices: Vec<i8>) -> Result<Complex<T>, &'static str> {
        GammaChain::new(
            self.vbar.clone(),
            self.u.clone(),
            indices,
            self.vectors.clone(),
        )
        .unwrap()
        .compute_chain()
    }

    pub fn compute_amplitude(&mut self) -> Result<Complex<T>, &'static str> {
        let mut numerator: Complex<T> = Complex::new(-T::one(), T::zero());

        //invariants
        let s12 = T::from_f64((self.external_kinematics[0] + self.external_kinematics[1]).square())
            .unwrap();
        let s23 = T::from_f64((self.external_kinematics[1] + self.external_kinematics[2]).square())
            .unwrap();
        let factor =
            T::from_f64((parameters(Parameters::alpha_ew) * 4.0 * std::f64::consts::PI).powf(2.0))
                .unwrap();

        /* ====== DIAGRAMS ====== */
        //D1
        let den1 = self.den[1];

        let indices = vec![6, 5, -1, 3, 7, 4, -1];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator += chain_res * (-factor / s23) * den1;

        //D2
        let indices = vec![-1, 2, 6, 3, -1, 5, 7];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator += chain_res * (-factor / s23) * &self.den[3];

        //D3
        let indices = vec![-1, 2, 6, 3, 7, 4, -1];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator += chain_res * (-factor);

        //D4
        let indices = vec![6, 5, -1, 3, -1, 5, 7];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator += chain_res * (-factor / s23 / s23) * (self.den[2] + self.den[4]);

        /* ====== COUNTERTERMS ====== */
        let den_uv = self.den[0] - self.mu_uv;
        let den_uv_sq = den_uv * den_uv;

        let box_den = self.den[0] * den1 * self.den[2] * self.den[3];
        //IR
        let indices = vec![-1, 2, 6, 5, 7, 4, -1];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator -= chain_res * (-factor / s23 / s23) * self.den[2];

        //UV1
        let indices = vec![6, 5, -1, 1, 7, 1, -1];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator -= chain_res * (-factor / s23 / s23) * box_den * utils::finv(den_uv * den_uv_sq);

        //UV2
        let indices = vec![-1, 1, 6, 1, -1, 5, 7];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator -= chain_res * (-factor / s23 / s23) * box_den * utils::finv(den_uv * den_uv_sq);

        //UV4
        let indices = vec![6, 5, -1, 1, -1, 5, 7];
        let mut chain_res = self.compute_chain(indices).unwrap() * utils::finv(den_uv_sq);
        let indices = vec![6, 5, -1, 1, 5, 1, -1, 5, 7];
        chain_res -= self.compute_chain(indices).unwrap() * utils::finv(den_uv * den_uv_sq);
        numerator -= chain_res * (-factor / s23 / s23) * box_den;

        //UVIR
        let indices = vec![-1, 1, 6, 5, 7, 1, -1];
        let chain_res = self.compute_chain(indices).unwrap();
        numerator -= chain_res * (factor / s23 / s23) * box_den * utils::finv(den_uv * den_uv_sq);

        return Ok(numerator);
    }
}
