use deformation::deformation::Deformer;
use deformation::parameterization::Parameterizer;
use integrand::integrands::Integrand;
use vector::LorentzVector;
use {Complex, REGION_ALL, REGION_EXT, REGION_INT};

pub struct Integrator {
    do_regions: bool,
    do_multichanneling: bool,
    parameterizer: Parameterizer,
    deformer: Deformer,
    integrand: Integrand,
}

impl Integrator {
    pub fn new(
        name: &str,
        do_regions: bool,
        do_multichanneling: bool,
        e_cm_sq: f64,
        mu_sq: f64,
        external_momenta: Vec<LorentzVector<f64>>,
    ) -> Result<Integrator, &str> {
        // TODO: compute qs and set qs!
        let mut qs = vec![LorentzVector::new()];
        let l = external_momenta.len() - 1;
        for (i, x) in external_momenta[..l].iter().enumerate() {
            let r = &qs[i] + x;
            qs.push(r);
        }

        let mut parameterizer = Parameterizer::new(e_cm_sq, 0, 0)?;
        parameterizer.set_qs(qs.clone());

        // TODO: for now the internal masses are assumed to be 0
        let mut deformer = Deformer::new(e_cm_sq, mu_sq, 0, vec![0.; qs.len()])?;
        deformer.set_qs(&qs);

        let mut integrand = Integrand::new(name, 0, 0, mu_sq)?;
        integrand.set_externals(external_momenta);

        Ok(Integrator {
            do_regions,
            do_multichanneling,
            parameterizer,
            deformer,
            integrand,
        })
    }

    pub fn evaluate(&mut self, k: &LorentzVector<f64>) -> Complex {
        if self.do_regions {
            let r1 = self.evaluate_region(k, REGION_EXT, 0, false);
            // r1_neg = self.evaluate_region(k, REGION_EXT, 0, true); // external opposite direction

            let mut r2 = Complex::new(0., 0.);
            if self.do_multichanneling {
                for i in 1..4 {
                    // TODO: self.qs.len()
                    r2 += self.evaluate_region(k, REGION_INT, i, false);
                }
            } else {
                r2 = self.evaluate_region(k, REGION_INT, 0, false);
            }

            //(r1 + r1_neg) / 2.0 + r2
            r1 + r2
        } else {
            if self.do_multichanneling {
                let mut r = Complex::new(0., 0.);
                for i in 1..4 {
                    r += self.evaluate_region(k, REGION_ALL, i, false);
                }
                r
            } else {
                self.evaluate_region(k, REGION_ALL, 0, false)
            }
        }
    }

    fn evaluate_region(
        &mut self,
        k: &LorentzVector<f64>,
        region: usize,
        channel: usize,
        minus: bool,
    ) -> Complex {
        self.parameterizer.set_channel(channel);
        if region == 1 {
            // external region
            self.parameterizer.set_mode("weinzierl").unwrap();
        } else {
            if channel == 0 {
                self.parameterizer.set_mode("log").unwrap();
            } else {
                self.parameterizer.set_mode("weinzierl").unwrap();
            }
        }

        self.parameterizer.set_region(region);
        let (mut k_mapped, mut jac_k) = self.parameterizer.map(k).unwrap();

        if minus {
            k_mapped = -k_mapped;
            jac_k = -jac_k;
        }

        self.deformer.set_region(region);
        let (ks, jac_deform) = self.deformer.deform(&k_mapped).unwrap();

        self.integrand.set_region(region);
        self.integrand.set_channel(channel);
        let out = self.integrand.evaluate(&[ks]).unwrap();

        let result = out * jac_deform * jac_k;

        // for Cuba, since it evaluates the boundaries
        if result.is_finite() {
            result
        } else {
            Complex::new(0., 0.)
        }
    }
}
