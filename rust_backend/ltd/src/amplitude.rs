use gamma_chain::GammaChain;
use num::Complex;
use utils;
use vector::LorentzVector;
use FloatLike;

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

enum Polarizations {
    UPlus,
    UMinus,
    UBarPlus,
    UBarMinus,
    APlus,
    AMinus,
}

pub fn compute_polarization(
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
                    let u_bar_plus = Vec::new();
                    for g in gamma0.iter() {
                        let res = Complex::default();
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
                    let u_bar_minus = Vec::new();
                    for g in gamma0.iter() {
                        let res = Complex::default();
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
            let a_plus = Vec::new();
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
            let a_minus = Vec::new();
            for i in 0..a1.len() {
                a_minus.push(norm * (a1[i] - ii * a2[i]));
            }

            return Ok(a_minus);
        }
    }
}
pub struct eeAA<T: FloatLike> {
    external_kinematics: Vec<LorentzVector<f64>>,
    den: Vec<Complex<T>>,
    vectors: Vec<LorentzVector<Complex<T>>>,
    vbar: Vec<Complex<T>>,
    u: Vec<Complex<T>>,
    a_mu: Vec<Complex<T>>,
    a_nu: Vec<Complex<T>>,
    mu_uv: Complex<T>,
}

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
        let a_mu = LorentzVector::new(1., 2., 3., 4.)
            .iter()
            .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
            .collect();

        let a_nu = LorentzVector::new(2., 3., 5., 7.)
            .iter()
            .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
            .collect();

        let k = loop_momentum;
        let ps = external_kinematics;
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

    pub fn compute_amplitude(&mut self) -> Result<Complex<T>, &'static str> {
        let numerator = Complex::default() - 1;

        //invariants
        let s12 = (self.external_kinematics[0] + self.external_kinematics[1]).square();
        let s23 = (self.external_kinematics[1] + self.external_kinematics[2]).square();
        let factor = (parameters(Parameters::alpha_ew) * 4.0 * std::f64::consts::PI).pow(2);

        /* ====== DIAGRAMS ====== */
        //D1
        let indices = vec![6, 5, -1, 3, 7, 4, -1];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator += -factor / s23 * chain_res * self.den[1];

        //D2
        let indices = vec![-1, 2, 6, 3, -1, 5, 7];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator += -factor / s23 * chain_res * self.den[3];

        //D3
        let indices = vec![-1, 2, 6, 3, 7, 4, -1];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator += -factor * chain_res;

        //D4
        let indices = vec![6, 5, -1, 3, -1, 5, 7];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator += -factor / s23.pow(2) * chain_res * (self.den[2] + self.den[4]);

        /* ====== COUNTERTERMS ====== */
        let den_uv = self.den[0] - self.mu_uv;
        let box_den = self.den[0] * self.den[1] * self.den[2] * self.den[3];
        //IR
        let indices = vec![-1, 2, 6, 5, 7, 4, -1];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator -= -factor / s23.pow(2) * chain_res * self.den[2];

        //UV1
        let indices = vec![6, 5, -1, 1, 7, 1, -1];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator -= -factor / s23.pow(2) * chain_res * box_den * utils::finv(den_uv.pow(3));

        //UV2
        let indices = vec![-1, 1, 6, 1, -1, 5, 7];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator -= -factor / s23.pow(2) * chain_res * box_den * utils::finv(den_uv.pow(3));

        //UV4
        let indices = vec![6, 5, -1, 1, -1, 5, 7];
        let mut chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap()
            * utils::finv(den_uv.pow(2));
        let indices = vec![6, 5, -1, 1, 5, 1, -1, 5, 7];
        chain_res -= GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap()
            * utils::finv(den_uv.pow(3));
        numerator -= -factor / s23.pow(2) * chain_res * box_den;

        //UVIR
        let indices = vec![-1, 1, 6, 5, 7, 1, -1];
        let chain_res = GammaChain::new(self.vbar, self.v, indices, self.vectors)
            .unwrap()
            .compute_chain()
            .unwrap();
        numerator -= factor / s23.pow(2) * chain_res * box_den * utils::finv(den_uv.pow(3));

        return Ok(numerator);
    }
}
