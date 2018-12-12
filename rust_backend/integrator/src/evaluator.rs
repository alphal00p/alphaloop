use aggregator::Settings;
use deformation::deformation::Deformer;
use deformation::deformation::{
    DIAGONAL_BOX_ID, DOUBLE_BOX_ID, DOUBLE_BOX_SB_ID, DOUBLE_TRIANGLE_ID,
    TRIANGLE_BOX_ALTERNATIVE_ID, TRIANGLE_BOX_ID,
};
use deformation::parameterization::Parameterizer;
use integrand::integrands::Integrand;
use vector::LorentzVector;
use Complex;

const ONE_LOOP_ID: usize = 1000;

#[derive(Debug, Deserialize)]
pub struct Topology {
    pub loops: usize,
    e_cm_sq: f64,
    external_momenta: Vec<Vec<f64>>,
    #[serde(default)]
    id: usize,
    #[serde(default)]
    name: String,
    #[serde(default)]
    on_shell_flag: usize,
}

impl Topology {
    pub fn build_evaluator(&self, settings: &Settings) -> Evaluator {
        let ext = self
            .external_momenta
            .iter()
            .map(|p| LorentzVector::from_slice(&p))
            .collect();
        match self.loops {
            1 => Evaluator::new(&self.name, self.e_cm_sq, settings, ext, self.on_shell_flag),
            2 => Evaluator::new_two_loop(self.id, self.e_cm_sq, settings, ext, self.on_shell_flag),
            _ => unreachable!("Invalid number of loops"),
        }
    }
}

pub struct UserData {
    pub evaluator: Vec<Evaluator>, // one evaluator per core
    pub running_max: f64,
}

#[inline]
pub fn integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    _nvec: usize,
    core: i32,
) -> Result<(), &'static str> {
    // master is -1
    let r = user_data.evaluator[(core + 1) as usize].evaluate(x)?;

    if r.is_finite() {
        if r.re.abs() > user_data.running_max {
            user_data.running_max = r.re.abs();
            //println!("max {:e}", r.re);
        }

        //println!("{:e}", r.re);
        f[0] = r.re;
    //f[1] = r.im;
    } else {
        if x.iter().any(|v| *v == 0.0 as f64 || *v == 1.0 as f64) {
            println!("Hitting the hypercube side");
        } else {
            println!("Bad point: {}", r);
        }
        f[0] = 0.;
        //f[1] = 0.;
    }

    Ok(())
}

#[derive(Clone)]
pub struct Evaluator {
    pub id: usize,
    pub parameterizer: Parameterizer,
    pub deformer: Deformer<f64>,
    pub integrand: Integrand,
    pub dual: bool,
    pub deform_eps: f64,
}

impl Evaluator {
    pub fn new(
        name: &str,
        e_cm_sq: f64,
        settings: &Settings,
        external_momenta: Vec<LorentzVector<f64>>,
        on_shell_flag: usize,
    ) -> Evaluator {
        let mut qs = vec![LorentzVector::new()];
        for (i, x) in external_momenta[1..].iter().enumerate() {
            let r = &qs[i] + x;
            qs.push(r);
        }

        let mut parameterizer = Parameterizer::new(e_cm_sq, settings.alpha, 0, 0).unwrap();
        parameterizer.set_qs(qs.clone());

        // TODO: for now the internal masses are assumed to be 0
        let mut deformer = Deformer::new(
            e_cm_sq,
            settings.mu_sq,
            settings.m1_fac,
            settings.m2_fac,
            settings.m3_fac,
            settings.m4_fac,
            settings.gamma1,
            settings.gamma2,
            settings.soft_fac,
            0,
            &vec![0.; qs.len()],
        )
        .unwrap();
        deformer.set_qs(&qs);

        let mut integrand = Integrand::new(name, 0, 0, settings.mu_sq).unwrap();
        integrand.set_externals(external_momenta);
        integrand.set_on_shell_flag(on_shell_flag);

        Evaluator {
            parameterizer,
            deformer,
            integrand,
            id: ONE_LOOP_ID,
            dual: settings.dual,
            deform_eps: settings.deform_eps,
        }
    }

    pub fn new_two_loop(
        id: usize,
        e_cm_sq: f64,
        settings: &Settings,
        external_momenta: Vec<LorentzVector<f64>>,
        on_shell_flag: usize,
    ) -> Evaluator {
        let mut parameterizer = Parameterizer::new(e_cm_sq, settings.alpha, 0, 0).unwrap();

        let mut deformer = Deformer::new(
            e_cm_sq,
            settings.mu_sq,
            settings.m1_fac,
            settings.m2_fac,
            settings.m3_fac,
            settings.m4_fac,
            settings.gamma1,
            settings.gamma2,
            settings.soft_fac,
            0,
            &[0.; 7],
        )
        .unwrap();

        let mut integrand = match id {
            DOUBLE_BOX_ID => {
                Integrand::new("box2L_direct_integration", 0, 0, settings.mu_sq).unwrap()
            }
            DOUBLE_TRIANGLE_ID => {
                Integrand::new("triangle2L_direct_integration", 0, 0, settings.mu_sq).unwrap()
            }
            TRIANGLE_BOX_ID => {
                Integrand::new("trianglebox_direct_integration", 0, 0, settings.mu_sq).unwrap()
            }
            TRIANGLE_BOX_ALTERNATIVE_ID => Integrand::new(
                "trianglebox_alternative_direct_integration",
                0,
                0,
                settings.mu_sq,
            )
            .unwrap(),
            DIAGONAL_BOX_ID => {
                Integrand::new("diagonalbox_direct_integration", 0, 0, settings.mu_sq).unwrap()
            }
            DOUBLE_BOX_SB_ID => {
                Integrand::new("box2L_direct_integration_SB", 0, 0, settings.mu_sq).unwrap()
            }
            _ => unreachable!("Unknown id"),
        };

        if id == TRIANGLE_BOX_ALTERNATIVE_ID
            || id == DOUBLE_BOX_ID
            || id == DIAGONAL_BOX_ID
            || id == DOUBLE_BOX_SB_ID
        {
            parameterizer.set_qs(vec![
                LorentzVector::new(),
                external_momenta[0].clone(),
                LorentzVector::new(),
                external_momenta[1].clone(),
                &external_momenta[1] + &external_momenta[2],
            ]);
        }

        deformer.set_external_momenta(&external_momenta);
        integrand.set_externals(external_momenta.clone());
        integrand.set_on_shell_flag(on_shell_flag);

        Evaluator {
            parameterizer,
            deformer,
            integrand,
            id,
            dual: settings.dual,
            deform_eps: settings.deform_eps,
        }
    }

    fn evaluate(&mut self, x: &[f64]) -> Result<Complex, &'static str> {
        if self.id == ONE_LOOP_ID {
            let (k_m, jac_k) = self.parameterizer.map(&LorentzVector::from_slice(x))?;

            let (d, j) = if self.dual {
                let (kk, j) = self.deformer.jacobian_using_dual(&k_m, false);
                (&kk.to_complex(false) + &k_m, j)
            } else {
                self.deformer.deform(&k_m)?
            };

            let v = self.integrand.evaluate(&[d])?;

            Ok(v * j * jac_k)
        } else {
            if self.id == TRIANGLE_BOX_ALTERNATIVE_ID
                || self.id == DOUBLE_BOX_ID
                || self.id == DIAGONAL_BOX_ID
                || self.id == DOUBLE_BOX_SB_ID
            {
                self.parameterizer.set_mode("weinzierl")?;
                self.parameterizer.set_channel(1);
            }
            let (k_m, jac_k) = self
                .parameterizer
                .map(&LorentzVector::from_slice(&x[..4]))?;
            if self.id == TRIANGLE_BOX_ALTERNATIVE_ID || self.id == DOUBLE_BOX_ID {
                // 3 here goes with integrand channel 1
                // 4 here goes with integrand channel 2
                self.parameterizer.set_channel(3);
            }
            if self.id == DIAGONAL_BOX_ID || self.id == DOUBLE_BOX_SB_ID {
                self.parameterizer.set_channel(4);
            }
            let (l_m, jac_l) = self
                .parameterizer
                .map(&LorentzVector::from_slice(&x[4..]))?;

            let (d, j) = if self.dual {
                // use the dual loop jacobian
                let (kk, ll, j) = self
                    .deformer
                    .jacobian_using_dual9_two_loops(self.id, &k_m, &l_m);
                ((kk, ll), j)
            } else {
                let d = self.deformer.deform_two_loops(self.id, &k_m, &l_m)?;
                let j =
                    self.deformer
                        .numerical_jacobian_two_loops(self.id, &k_m, &l_m, (&d.0, &d.1));
                (d, j)
            };

            if self.deform_eps > 0. {
                let deform_magnitude = (d.0.imag().euclidean_square()
                    + d.1.imag().euclidean_square())
                    / (k_m.euclidean_square() + l_m.euclidean_square());
                if deform_magnitude < self.deform_eps {
                    return Ok(Complex::default());
                }
            }

            if self.id == TRIANGLE_BOX_ALTERNATIVE_ID {
                // disable two-loop channels for now
                // self.integrand.set_channel(1);
            }

            let v = self.integrand.evaluate(&[d.0, d.1])?;

            Ok(v * j * jac_k * jac_l)
        }
    }

    pub fn evaluate_two_loop_test_function(
        &mut self,
        k: &LorentzVector<f64>,
        l: &LorentzVector<f64>,
        do_deformation: bool,
    ) -> Complex {
        let (k_m, jac_k) = self.parameterizer.map(&k).unwrap();
        let (l_m, jac_l) = self.parameterizer.map(&l).unwrap();

        let (d, j) = if do_deformation {
            let d = self.deformer.deform_two_loops(self.id, &k_m, &l_m).unwrap();
            let j = self
                .deformer
                .numerical_jacobian_two_loops(self.id, &k_m, &l_m, (&d.0, &d.1));
            (d, j)
        } else {
            (
                (k_m.to_complex(true), l_m.to_complex(true)),
                Complex::new(1., 0.),
            )
        };

        // yields pi/400 for the euclidean and -pi/400 for the non-euclidean
        let mu = Complex::new(0.0, -1.); // needs to be more than default -1e9
        let factor = mu.im.powf(2. + 6.);
        let denom = (d.0.square() - mu).powf(6.) * (d.1.square() - mu).powf(6.);

        //println!("k={}, l={}, factor={:e}, inv_denom={:e}, j={:e}, j_k={:e}, j_l={:e}", d.0, d.1, factor, 1./denom, j, jac_k, jac_l);
        //println!("{:e} {:e} {:e}", d.0.square(), d.0.square() - mu, (d.0.square() - mu).powf(8.));
        //println!("{:e} {:e} {:e}", d.1.square(), d.1.square() - mu, (d.1.square() - mu).powf(8.));

        // note that there is a precision issue for complex division!
        factor / denom * j * jac_k * jac_l
    }
}
