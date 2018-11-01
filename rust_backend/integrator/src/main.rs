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

use integrator::aggregator::{Aggregator, IntegrationSettings};
use integrator::evaluator::Evaluator;

fn double_triangle_integrator() -> Evaluator {
    let mu_sq = -1e9;
    let e_cm_sq = 1.0;
    let external_momenta = vec![LorentzVector::from(1.0, 0., 0., 0.)];

    println!(
        "Starting integration of double triangle with externals {:#?}",
        external_momenta
    );

    Evaluator::new_two_loop(DOUBLE_TRIANGLE_ID, e_cm_sq, mu_sq, external_momenta)
}

fn triangle_box_integrator() -> Evaluator {
    let mu_sq = -1e9;
    let e_cm_sq = 1.0;
    let external_momenta = vec![
        LorentzVector::from(1.0, 0., 0., 0.),
        -LorentzVector::from(19. / 32., 0., 0., 105f64.sqrt() / 32.),
        -LorentzVector::from(1. - 19. / 32., 0., 0., -105f64.sqrt() / 32.),
    ];

    /*let external_momenta = vec![
        LorentzVector::from(-90., 0., 0., 0.) * 0.01111,
        LorentzVector::from(19.6586, -7.15252, -0.206016, 8.96383) * 0.01111,
        LorentzVector::from(26.874, 7.04203, -0.0501295, -12.9055) * 0.01111,
        LorentzVector::from(43.4674, 0.110491, 0.256146, 3.9417) * 0.01111,
    ];*/

    println!(
        "Starting integration of triangle box with externals {:#?}",
        external_momenta
    );

    Evaluator::new_two_loop(TRIANGLE_BOX_ID, e_cm_sq, mu_sq, external_momenta)
}

fn box_integrator(loops: usize) -> Evaluator {
    let mu_sq = -1e9;
    let e_cm_sq = 0.947492;

    // seed 2
    let external_momenta = vec![
        LorentzVector::from(
            -4.7809952420694083e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            -4.6752449673455959e+02,
        ) * 0.001,
        LorentzVector::from(
            -5.0850678957797919e+02,
            -0.0000000000000000e+00,
            -0.0000000000000000e+00,
            4.6752449673455959e+02,
        ) * 0.001,
        LorentzVector::from(
            4.5782801395958194e+02,
            1.3758384614384497e+02,
            8.1217573038820291e+01,
            -3.0672606911725950e+02,
        ) * 0.001,
        LorentzVector::from(
            5.2877829982533808e+02,
            -1.3758384614384497e+02,
            -8.1217573038820291e+01,
            3.0672606911725950e+02,
        ) * 0.001,
    ];

    if loops == 1 {
        println!(
            "Starting integration of one-loop box with externals {:#?}",
            external_momenta
        );

        Evaluator::new("box1L_direct_integration", e_cm_sq, mu_sq, external_momenta)
    } else {
        println!(
            "Starting integration of double box with externals {:#?}",
            external_momenta
        );

        Evaluator::new_two_loop(DOUBLE_BOX_ID, e_cm_sq, mu_sq, external_momenta)
    }
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
            Arg::with_name("loops")
                .short("l")
                .long("loops")
                .value_name("LOOPS")
                .default_value("2")
                .help("Set the number of loops (only used for the box)"),
        )
        .arg(
            Arg::with_name("multichanneling")
                .short("m")
                .long("multichanneling")
                .help("Use multichanneling"),
        )
        .arg(
            Arg::with_name("regions")
                .short("r")
                .long("regions")
                .help("Use regions"),
        )
        .arg(
            Arg::with_name("samples")
                .short("s")
                .long("samples")
                .value_name("SAMPLES")
                .default_value("20000000")
                .help("Number of samples per integration"),
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
    let loops = usize::from_str(matches.value_of("loops").unwrap()).unwrap();
    let samples = usize::from_str(matches.value_of("samples").unwrap()).unwrap();

    let multichanneling = matches.is_present("multichanneling");
    let regions = matches.is_present("regions");

    let settings = IntegrationSettings {
        cores,
        samples,
        param_mode: "linear".to_owned(),
    };

    let i = match matches.value_of("topology").unwrap() {
        "box" => box_integrator(loops),
        "triangle" => double_triangle_integrator(),
        "trianglebox" => triangle_box_integrator(),
        _ => unreachable!(),
    };

    let mut aggregator = Aggregator::new(loops, regions, multichanneling, i, settings);
    let r = aggregator.aggregate();
    println!("Result: {:e} +- {:e}", r.0, r.1);
}
