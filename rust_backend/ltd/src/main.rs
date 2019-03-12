extern crate arrayvec;
extern crate clap;
extern crate cuba;
extern crate dual_num;
extern crate itertools;
extern crate ltd as ltdlib;
extern crate nalgebra as na;
extern crate num;
extern crate num_traits;
extern crate rand;
extern crate serde;
extern crate serde_yaml;
use serde::{Deserialize, Serialize};
extern crate vector;

use clap::{App, Arg};
use rand::prelude::*;
use std::str::FromStr;
use std::time::Instant;

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};

mod cts;
mod ltd;
mod topologies;
mod utils;

use ltdlib::{IntegratedPhase, Settings};

#[derive(Serialize, Deserialize)]
struct CubaResultDef {
    pub neval: i64,
    pub fail: i32,
    pub result: Vec<f64>,
    pub error: Vec<f64>,
    pub prob: Vec<f64>,
}

impl CubaResultDef {
    fn new(o: &CubaResult) -> CubaResultDef {
        CubaResultDef {
            neval: o.neval,
            fail: o.fail,
            result: o.result.clone(),
            error: o.error.clone(),
            prob: o.prob.clone(),
        }
    }
}

#[derive(Debug)]
struct UserData {
    topo: Vec<topologies::Topology>,
    sample_count: usize,
    integrated_phase: IntegratedPhase,
}

#[inline(always)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    _nvec: usize,
    core: i32,
) -> Result<(), &'static str> {
    let res = user_data.topo[(core + 1) as usize].evaluate(x);
    user_data.sample_count += 1;

    if user_data.sample_count % 100000 == 0 {
        println!("Sample: {:?} {:e}", x, res);
    }

    if res.re.is_finite() {
        match user_data.integrated_phase {
            IntegratedPhase::Real => {
                f[0] = res.re;
            }
            IntegratedPhase::Imag => {
                f[0] = res.im;
            }
            IntegratedPhase::Both => {
                f[0] = res.re;
                f[1] = res.im;
            }
        }
    } else {
        println!("Bad point: {:?}", x);
        f[0] = 0.;
    }

    Ok(())
}

fn bench(topo: &topologies::Topology, max_eval: usize) {
    let mut x = vec![0.; 3 * topo.n_loops];
    let mut rng = rand::thread_rng();

    let now = Instant::now();
    for _ in 0..max_eval {
        for xi in x.iter_mut() {
            *xi = rng.gen();
        }

        let _r = topo.evaluate(&x);
    }

    println!("{:#?}", now.elapsed());
}

fn main() {
    let matches = App::new("Feynman diagram integrator")
        .version("0.1")
        .about("Numerically integrate your favourite integrals with LTD")
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
            Arg::with_name("topologies")
                .short("l")
                .long("topologies")
                .value_name("TOPOLOGY_FILE")
                .default_value("../LTD/topologies.yaml")
                .help("Set the topology file"),
        )
        .arg(
            Arg::with_name("config")
                .short("f")
                .long("config")
                .value_name("CONFIG_FILE")
                .default_value("../LTD/hyperparameters.yaml")
                .help("Set the configuration file"),
        )
        .arg(
            Arg::with_name("deformation")
                .short("d")
                .long("deformation")
                .value_name("DEFORMATION")
                .default_value("none")
                .help("Set the deformation"),
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
                .help("Set the active topology"),
        )
        .get_matches();

    let mut settings = Settings::from_file(matches.value_of("config").unwrap());

    let mut cores = 1;
    if let Some(x) = matches.value_of("cores") {
        cores = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("samples") {
        settings.integrator.n_max = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("topology") {
        settings.general.topology = x.to_owned();
    }

    let topology_file = matches.value_of("topologies").unwrap();

    if let Some(x) = matches.value_of("deformation_strategy") {
        settings.general.deformation_strategy = x.to_owned();
    }
    let mut ci = CubaIntegrator::new(integrand);

    ci.set_mineval(10)
        .set_nstart(settings.integrator.n_start as i64)
        .set_nincrease(settings.integrator.n_increase as i64)
        .set_maxeval(settings.integrator.n_max as i64)
        .set_epsabs(0.)
        .set_epsrel(1e-15)
        .set_seed(1)
        .set_cores(cores, 1000);

    // load the example file
    let topologies = topologies::Topology::from_file(topology_file, &settings);
    let topo = topologies
        .get(&settings.general.topology)
        .expect("Unknown topology");

    if matches.is_present("bench") {
        bench(&topo, settings.integrator.n_max);
        return;
    }

    println!(
        "Integrating {} with {} samples and deformation {}",
        settings.general.topology, settings.integrator.n_max, settings.general.deformation_strategy
    );

    let cuba_result = match settings.integrator.integrator.as_ref() {
        "vegas" => ci.vegas(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            CubaVerbosity::Progress,
            0,
            UserData {
                topo: vec![topo.clone(); cores + 1],
                sample_count: 0,
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
        "cuhre" => ci.cuhre(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            CubaVerbosity::Progress,
            UserData {
                topo: vec![topo.clone(); cores + 1],
                sample_count: 0,
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
        x => panic!("Unknown integrator {}", x),
    };
    println!("{:#?}", cuba_result);

    let f = OpenOptions::new()
        .create(true)
        .write(true)
        .open(settings.general.topology.clone() + "_res.dat")
        .expect("Unable to create result file");
    let mut result_file = BufWriter::new(f);

    // write the result to a file
    writeln!(
        &mut result_file,
        "{}",
        serde_yaml::to_string(&CubaResultDef::new(&cuba_result)).unwrap()
    )
    .unwrap();
    writeln!(&mut result_file, "...").unwrap(); // write end-marker, for easy streaming
}
