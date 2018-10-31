extern crate clap;
extern crate cuba;
extern crate deformation;
extern crate integrand;
extern crate integrator;
extern crate num;
extern crate rand;
extern crate vector;

use clap::{App, Arg};
use deformation::deformation::{Deformer, DOUBLE_BOX_ID, DOUBLE_TRIANGLE_ID, TRIANGLE_BOX_ID};
use deformation::parameterization::Parameterizer;
use integrand::integrands::Integrand;
use rand::prelude::*;
use std::str::FromStr;
use vector::LorentzVector;

use cuba::{CubaIntegrator, CubaVerbosity};
use integrator::integrator::Integrator;

struct UserData {
    evaluator: Vec<Integrator>, // one evaluator per core
    running_max: f64,
}

#[inline]
fn integrand(x: &[f64], f: &mut [f64], user_data: &mut UserData, _nvec: usize, core: i32) -> i32 {
    // master is -1
    let r = user_data.evaluator[(core + 1) as usize].evaluate_two_loop(
        &LorentzVector::from_slice(&x[..4]),
        &LorentzVector::from_slice(&x[4..]),
    );

    if r.is_finite() {
        if r.re.abs() > user_data.running_max {
            user_data.running_max = r.re.abs();
            println!("max {:e}", r.re);
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

fn double_triangle_integrator() -> Integrator {
    let mu_sq = -1e9;
    let e_cm_sq = 1.0;
    let external_momenta = vec![LorentzVector::from(1.0, 0., 0., 0.)];

    println!(
        "Starting integration of double triangle with externals {:#?}",
        external_momenta
    );

    Integrator::new_two_loop(DOUBLE_TRIANGLE_ID, e_cm_sq, mu_sq, external_momenta)
}

fn triangle_box_integrator() -> Integrator {
    let mu_sq = -1e9;
    let e_cm_sq = 1.0;
    let external_momenta = vec![
        LorentzVector::from(1.0, 0., 0., 0.),
        -LorentzVector::from(19. / 32., 0., 0., 105f64.sqrt() / 32.),
        -LorentzVector::from(1. - 19. / 32., 0., 0., -105f64.sqrt() / 32.),
    ];

    println!(
        "Starting integration of triangle box with externals {:#?}",
        external_momenta
    );

    Integrator::new_two_loop(TRIANGLE_BOX_ID, e_cm_sq, mu_sq, external_momenta)
}

fn double_box_integrator() -> Integrator {
    let mu_sq = -1e9;
    let e_cm_sq = 1.0;

    // seed 2
    let external_momenta = vec![
        LorentzVector::from(
            -4.7809952420694083e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            -4.6752449673455959e+02,
        ) * 0.0001,
        LorentzVector::from(
            -5.0850678957797919e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            4.6752449673455959e+02,
        ) * 0.0001,
        LorentzVector::from(
            -5.0850678957797919e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            4.6752449673455959e+02,
        ) * 0.0001,
        LorentzVector::from(
            4.5782801395958194e+02,
            1.3758384614384497e+02,
            8.1217573038820291e+01,
            -3.0672606911725950e+02,
        ) * 0.0001,
        LorentzVector::from(
            5.2877829982533808e+02,
            -1.3758384614384497e+02,
            -8.1217573038820291e+01,
            3.0672606911725950e+02,
        ) * 0.0001,
    ];

    println!(
        "Starting integration of double box with externals {:#?}",
        external_momenta
    );

    Integrator::new_two_loop(DOUBLE_BOX_ID, e_cm_sq, mu_sq, external_momenta)
}

fn main() {
    let matches = App::new("Feynman diagram integrator")
        .version("0.1")
        .about("Numerically integrate your favourite integrals")
        .arg(
            Arg::with_name("cores")
                .short("c")
                .long("cores")
                .value_name("NUMCORES")
                .default_value("4")
                .help("Set the number of cores"),
        )
        .arg(
            Arg::with_name("config")
                .short("f")
                .long("config")
                .value_name("CONFIG_FILE")
                .help("Set the configuration file"),
        )
        .arg(
            Arg::with_name("bench")
                .long("bench")
                .help("Run a benchmark instead"),
        )
        .arg(
            Arg::with_name("topology")
                .short("t")
                .long("topology")
                .value_name("TOPOLOGY")
                .default_value("box")
                .possible_values(&["triangle", "box", "trianglebox"])
                .help("Set the topology"),
        )
        .get_matches();

    if matches.is_present("bench") {
        performance_test();
        return;
    }

    let cores = usize::from_str(matches.value_of("cores").unwrap()).unwrap();

    let mut ci = CubaIntegrator::new(integrand);
    ci.set_epsabs(0.)
        .set_mineval(10)
        .set_nstart(10000)
        .set_nincrease(50000)
        .set_maxeval(600000000)
        .set_pseudo_random(true)
        .set_cores(cores, 1000);

    let i = match matches.value_of("topology").unwrap() {
        "box" => double_box_integrator(),
        "triangle" => double_triangle_integrator(),
        "trianglebox" => triangle_box_integrator(),
        _ => unreachable!(),
    };

    let r = ci.vegas(
        8,
        1,
        CubaVerbosity::Progress,
        UserData {
            evaluator: vec![i; cores + 1],
            running_max: 0f64,
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

        let (k_m, _jac_k) = parameterizer.map(&k).unwrap();
        let (l_m, _jac_l) = parameterizer.map(&l).unwrap();

        let d = deformer
            .deform_two_loops(DOUBLE_TRIANGLE_ID, &k_m, &l_m)
            .unwrap();

        let _j =
            deformer.numerical_jacobian_two_loops(DOUBLE_TRIANGLE_ID, &k_m, &l_m, (&d.0, &d.1));

        integrand.evaluate(&[d.0, d.1]).unwrap();
    }
}
