extern crate cuba;
extern crate deformation;
extern crate integrand;
extern crate integrator;
extern crate num;
extern crate rand;
extern crate vector;

use deformation::deformation::{Deformer, DOUBLE_BOX_ID, DOUBLE_TRIANGLE_ID};
use deformation::parameterization::Parameterizer;
use integrand::integrands::Integrand;
use rand::prelude::*;
use vector::LorentzVector;

use cuba::{CubaIntegrator, CubaVerbosity};
use integrator::integrator::Integrator;

struct UserData {
    f: [Integrator; 4], // number of cores
}

#[inline]
fn test_integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    _nvec: usize,
    core: i32,
) -> i32 {
    let r = user_data.f[core as usize].evaluate_two_loop(
        &LorentzVector::from_slice(&x[..4]),
        &LorentzVector::from_slice(&x[4..]),
    );

    if r.is_finite() {
        f[0] = r.im;
        //f[1] = r.im;
    } else {
        f[0] = 0.;
        //f[1] = 0.;
    }

    0
}

fn test_double_triangle_integrator() -> Integrator {
    let mu_sq = -1e9;
    let e_cm_sq = 1000. * 1000.;
    let external_momenta = vec![LorentzVector::from(
        (e_cm_sq
            + 1.3758384614384497e+02_f64.powi(2)
            + 8.1217573038820291e+01_f64.powi(2)
            + 3.0672606911725950e+02_f64.powi(2))
        .sqrt(),
        1.3758384614384497e+02,
        8.1217573038820291e+01,
        -3.0672606911725950e+02,
    )];

    Integrator::new_two_loop(DOUBLE_TRIANGLE_ID, e_cm_sq, mu_sq, external_momenta)
}

fn test_double_box_integrator() -> Integrator {
    let mu_sq = -1e9;
    let e_cm_sq = 1000. * 1000.;

    // seed 2
    let external_momenta = vec![
        LorentzVector::from(
            -4.7809952420694083e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            -4.6752449673455959e+02,
        ),
        LorentzVector::from(
            -5.0850678957797919e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            4.6752449673455959e+02,
        ),
        LorentzVector::from(
            -5.0850678957797919e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            4.6752449673455959e+02,
        ),
        LorentzVector::from(
            4.5782801395958194e+02,
            1.3758384614384497e+02,
            8.1217573038820291e+01,
            -3.0672606911725950e+02,
        ),
        LorentzVector::from(
            5.2877829982533808e+02,
            -1.3758384614384497e+02,
            -8.1217573038820291e+01,
            3.0672606911725950e+02,
        ),
    ];

    Integrator::new_two_loop(DOUBLE_BOX_ID, e_cm_sq, mu_sq, external_momenta)
}

fn main() {
    //performance_test();
    //return;

    let mut ci = CubaIntegrator::new(test_integrand);
    ci.set_epsabs(0.).set_mineval(10).set_nstart(10000).set_nincrease(5000).set_maxeval(600000000);

    //let i = test_double_box_integrator();
    let i = test_double_triangle_integrator();

    let r = ci.vegas(
        8,
        1,
        CubaVerbosity::Progress,
        UserData {
            f: [i.clone(), i.clone(), i.clone(), i.clone()],
        },
    );
    println!("{:#?}", r);
}

fn performance_test() {
    // do a performance test
    let mu_sq = -1e9;
    let e_cm_sq = 1000. * 1000.;
    let external_momenta = vec![LorentzVector::from(
        (e_cm_sq
            + 1.3758384614384497e+02_f64.powi(2)
            + 8.1217573038820291e+01_f64.powi(2)
            + 3.0672606911725950e+02_f64.powi(2))
        .sqrt(),
        1.3758384614384497e+02,
        8.1217573038820291e+01,
        -3.0672606911725950e+02,
    )];

    let mut parameterizer = Parameterizer::new(e_cm_sq, 0, 0).unwrap();
    parameterizer.set_mode("log").unwrap();
    parameterizer.set_qs(external_momenta.clone()); // should be Q

    let mut deformer = Deformer::new(e_cm_sq, mu_sq, 0, vec![0.; external_momenta.len()]).unwrap();
    let mut integrand = Integrand::new("triangle2L_direct_integration", 0, 0, mu_sq).unwrap();

    integrand.set_externals(external_momenta.clone());
    deformer.set_external_momenta(external_momenta);

    let mut rng = thread_rng();

    for _ in 0..200000 {
        let k = LorentzVector::from(rng.gen(), rng.gen(), rng.gen(), rng.gen());
        let l = LorentzVector::from(rng.gen(), rng.gen(), rng.gen(), rng.gen());

        let (k_m, jac_k) = parameterizer.map(&k).unwrap();
        let (l_m, jac_l) = parameterizer.map(&l).unwrap();

        let d = deformer
            .deform_two_loops(DOUBLE_TRIANGLE_ID, &k_m, &l_m)
            .unwrap();
        let j = deformer.numerical_jacobian_two_loops(DOUBLE_TRIANGLE_ID, &k_m, &l_m, (&d.0, &d.1));
        
        integrand.evaluate(&[d.0, d.1]).unwrap();
    }
}
