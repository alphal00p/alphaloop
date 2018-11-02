use deformation::deformation::Deformer;
use deformation::deformation::{DOUBLE_BOX_ID, DOUBLE_TRIANGLE_ID, TRIANGLE_BOX_ID};
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
}

impl Topology {
    pub fn build_evaluator(&self, mu_sq: f64, alpha: Option<f64>) -> Evaluator {
        let ext = self
            .external_momenta
            .iter()
            .map(|p| LorentzVector::from_slice(&p))
            .collect();

        match self.loops {
            1 => Evaluator::new(&self.name, self.e_cm_sq, alpha, mu_sq, ext),
            2 => Evaluator::new_two_loop(self.id, self.e_cm_sq, alpha, mu_sq, ext),
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
) -> i32 {
    // master is -1
    let r = user_data.evaluator[(core + 1) as usize].evaluate(x);

    if r.is_finite() {
        if r.re.abs() > user_data.running_max {
            user_data.running_max = r.re.abs();
            //println!("max {:e}", r.re);
        }

        //println!("{:e}", r.re);
        f[0] = r.re;
    //f[1] = r.im;
    } else {
        println!("Bad point: {}", r);
        f[0] = 0.;
        //f[1] = 0.;
    }

    0
}

#[derive(Clone)]
pub struct Evaluator {
    pub id: usize,
    pub parameterizer: Parameterizer,
    pub deformer: Deformer,
    pub integrand: Integrand,
}

impl Evaluator {
    pub fn new(
        name: &str,
        e_cm_sq: f64,
        alpha: Option<f64>,
        mu_sq: f64,
        external_momenta: Vec<LorentzVector<f64>>,
    ) -> Evaluator {
        let mut qs = vec![LorentzVector::new()];
        for (i, x) in external_momenta[1..].iter().enumerate() {
            let r = &qs[i] + x;
            qs.push(r);
        }

        let mut parameterizer = Parameterizer::new(e_cm_sq, alpha, 0, 0).unwrap();
        parameterizer.set_qs(qs.clone());

        // TODO: for now the internal masses are assumed to be 0
        let mut deformer = Deformer::new(e_cm_sq, mu_sq, 0, vec![0.; qs.len()]).unwrap();
        deformer.set_qs(&qs);

        let mut integrand = Integrand::new(name, 0, 0, mu_sq).unwrap();
        integrand.set_externals(external_momenta);

        Evaluator {
            parameterizer,
            deformer,
            integrand,
            id: ONE_LOOP_ID,
        }
    }

    pub fn new_two_loop(
        id: usize,
        e_cm_sq: f64,
        alpha: Option<f64>,
        mu_sq: f64,
        external_momenta: Vec<LorentzVector<f64>>,
    ) -> Evaluator {
        let parameterizer = Parameterizer::new(e_cm_sq, alpha, 0, 0).unwrap();

        let mut deformer = Deformer::new(e_cm_sq, mu_sq, 0, vec![0.; 7]).unwrap();

        let mut integrand = match id {
            DOUBLE_BOX_ID => Integrand::new("box2L_direct_integration", 0, 0, mu_sq).unwrap(),
            DOUBLE_TRIANGLE_ID => {
                Integrand::new("triangle2L_direct_integration", 0, 0, mu_sq).unwrap()
            }
            TRIANGLE_BOX_ID => {
                Integrand::new("trianglebox_direct_integration", 0, 0, mu_sq).unwrap()
            }
            _ => unreachable!("Unknown id"),
        };

        integrand.set_externals(external_momenta.clone());
        deformer.set_external_momenta(external_momenta);

        Evaluator {
            parameterizer,
            deformer,
            integrand,
            id,
        }
    }

    fn evaluate(&mut self, x: &[f64]) -> Complex {
        if self.id == ONE_LOOP_ID {
            let (k_m, jac_k) = self
                .parameterizer
                .map(&LorentzVector::from_slice(x))
                .unwrap();

            let (d, j) = self.deformer.deform(&k_m).unwrap();
            let v = self.integrand.evaluate(&[d]).unwrap();

            v * j * jac_k
        } else {
            let (k_m, jac_k) = self
                .parameterizer
                .map(&LorentzVector::from_slice(&x[..4]))
                .unwrap();
            let (l_m, jac_l) = self
                .parameterizer
                .map(&LorentzVector::from_slice(&x[4..]))
                .unwrap();

            let d = self.deformer.deform_two_loops(self.id, &k_m, &l_m).unwrap();
            let j = self
                .deformer
                .numerical_jacobian_two_loops(self.id, &k_m, &l_m, (&d.0, &d.1));

            let v = self.integrand.evaluate(&[d.0, d.1]).unwrap();

            v * j * jac_k * jac_l
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
