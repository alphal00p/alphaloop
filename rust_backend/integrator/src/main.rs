extern crate clap;
extern crate cuba;
extern crate deformation;
extern crate integrand;
extern crate integrator;
extern crate num;
extern crate rand;
extern crate vector;

use clap::{App, Arg};
use deformation::deformation::{Deformer, DOUBLE_TRIANGLE_ID};
use deformation::parameterization::Parameterizer;
use integrand::integrands::Integrand;
use rand::prelude::*;
use std::str::FromStr;
use vector::LorentzVector;

use integrator::aggregator::{Aggregator, Settings};

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

    let mut parameterizer = Parameterizer::new(e_cm_sq, None, 0, 0).unwrap();
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
                .help("Set the number of cores"),
        )
        .arg(
            Arg::with_name("samples")
                .short("s")
                .long("samples")
                .value_name("SAMPLES")
                .help("Number of samples per integration"),
        )
        .arg(
            Arg::with_name("survey_n_points")
                .short("sn")
                .long("survey_n_points")
                .value_name("SURVEY_N_POINTS")
                .help("Number of samples per survey iteration"),
        )
        .arg(
            Arg::with_name("survey_n_iterations")
                .short("si")
                .long("survey_n_iterations")
                .value_name("SURVEY_N_ITERATIONS")
                .help("Number of iterations for the survey"),
        )
        .arg(
            Arg::with_name("refine_n_points")
                .short("rn")
                .long("refine_n_points")
                .value_name("REFINE_N_POINTS")
                .help("Number of samples per refine run"),
        )
        .arg(
            Arg::with_name("refine_n_runs")
                .short("rr")
                .long("refine_n_runs")
                .value_name("REFINE_N_RUNS")
                .help("Number of independent refine runs"),
        )
        .arg(
            Arg::with_name("config")
                .short("f")
                .long("config")
                .value_name("CONFIG_FILE")
                .default_value("config.yaml")
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
                .possible_values(&[
                    "box",
                    "double-triangle",
                    "double-box",
                    "double-box-SB",
                    "triangle-box",
                    "triangle-box-alternative",
                    "cross-box"
                ])
                .help("Set the active topology"),
        )
        .get_matches();

    if matches.is_present("bench") {
        performance_test();
    }

    let mut p = Settings::from_file(matches.value_of("config").unwrap());

    if let Some(x) = matches.value_of("topology") {
        p.active_topology = x.to_owned();
    }

    if let Some(x) = matches.value_of("cores") {
        p.cores = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("samples") {
        p.max_eval = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("survey_n_points") {
        p.survey_n_points = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("survey_n_iterations") {
        p.survey_n_iterations = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("refine_n_points") {
        p.refine_n_points = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("refine_n_runs") {
        p.refine_n_runs = usize::from_str(x).unwrap();
    }

    let mut aggregator = Aggregator::new(p);
    let r = aggregator.aggregate();
    println!("Result: {:e} +- {:e}", r.0, r.1);
}
