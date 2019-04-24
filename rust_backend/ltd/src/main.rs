extern crate arrayvec;
extern crate clap;
extern crate cuba;
extern crate dual_num;
extern crate itertools;
extern crate ltd;
extern crate nalgebra as na;
extern crate num;
extern crate num_traits;
extern crate rand;
extern crate serde;
extern crate serde_yaml;
use serde::{Deserialize, Serialize};
extern crate vector;

use clap::{App, Arg, SubCommand};
use num_traits::ToPrimitive;
use rand::prelude::*;
use std::str::FromStr;
use std::time::Instant;

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};

use ltd::integrand::Integrand;
use ltd::topologies::{LTDCache, Topology};
use ltd::{float, IntegratedPhase, Integrator, Settings};

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

struct UserData {
    integrand: Vec<Integrand>,
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
    let res = user_data.integrand[(core + 1) as usize].evaluate(x);

    if res.is_finite() {
        match user_data.integrated_phase {
            IntegratedPhase::Real => {
                f[0] = res.re.to_f64().unwrap();
            }
            IntegratedPhase::Imag => {
                f[0] = res.im.to_f64().unwrap();
            }
            IntegratedPhase::Both => {
                f[0] = res.re.to_f64().unwrap();
                f[1] = res.im.to_f64().unwrap();
            }
        }
    } else {
        f[0] = 0.;
    }

    Ok(())
}

fn bench(topo: &Topology, settings: &Settings) {
    let mut x = vec![0.; 3 * topo.n_loops];
    let mut rng = rand::thread_rng();

    let mut cache = LTDCache::<float>::new(&topo);

    let now = Instant::now();
    for _ in 0..settings.integrator.n_max {
        for xi in x.iter_mut() {
            *xi = rng.gen();
        }

        let _r = topo.evaluate(&x, &mut cache);
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
            Arg::with_name("topology")
                .short("t")
                .long("topology")
                .value_name("TOPOLOGY")
                .help("Set the active topology"),
        )
        .subcommand(SubCommand::with_name("bench").about("Run a benchmark"))
        .subcommand(
            SubCommand::with_name("inspect")
                .about("Inspect a single input point")
                .arg(Arg::with_name("point").required(true).min_values(3))
                .arg(
                    Arg::with_name("use_f128")
                        .short("f128")
                        .long("use_f128")
                        .help("Use f128 evaluation"),
                )
                .arg(
                    Arg::with_name("full_integrand")
                        .long("full_integrand")
                        .help("Evaluate the integrand and possibly its rotated vesion"),
                ),
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
        settings.general.deformation_strategy = x.into();
    }
    let mut ci = CubaIntegrator::new(integrand);

    ci.set_mineval(10)
        .set_nstart(settings.integrator.n_start as i64)
        .set_nincrease(settings.integrator.n_increase as i64)
        .set_maxeval(settings.integrator.n_max as i64)
        .set_epsabs(0.)
        .set_epsrel(1e-15)
        .set_seed(settings.integrator.seed)
        .set_cores(cores, 1000);

    // load the example file
    let mut topologies = Topology::from_file(topology_file, &settings);
    let topo = topologies
        .get_mut(&settings.general.topology)
        .expect("Unknown topology");
    topo.process();

    if let Some(_) = matches.subcommand_matches("bench") {
        bench(&topo, &settings);
        return;
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        let pt: Vec<_> = matches
            .values_of("point")
            .unwrap()
            .map(|x| f64::from_str(x).unwrap())
            .collect();
        if pt.len() != 3 * topo.n_loops {
            panic!(
                "Dimension of the input point is incorrect. It should be {} but is {}.",
                topo.n_loops * 3,
                pt.len()
            );
        }

        if matches.is_present("full_integrand") {
            settings.general.screen_log_core = Some(1);
            settings.general.log_points_to_screen = true;
            let mut i = Integrand::new(&topo, settings.clone(), 1);
            i.evaluate(&pt);
            return;
        }

        // TODO: prevent code repetition
        if matches.is_present("use_f128") {
            let mut cache = LTDCache::<f128::f128>::new(&topo);
            let (x, k_def, jac_para, jac_def, result) = topo.clone().evaluate(&pt, &mut cache);
            match topo.n_loops {
                1 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result, x, k_def[0], jac_para, jac_def
                    );
                }
                2 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result, x, k_def[0], k_def[1], jac_para, jac_def
                        );
                }
                3 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    );
                }
                4 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}n  | n={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], k_def[3], jac_para, jac_def
                    );
                }
                _ => {}
            }
        } else {
            let mut cache = LTDCache::<float>::new(&topo);
            let (x, k_def, jac_para, jac_def, result) = topo.clone().evaluate(&pt, &mut cache);
            match topo.n_loops {
                1 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result, x, k_def[0], jac_para, jac_def
                    );
                }
                2 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result, x, k_def[0], k_def[1], jac_para, jac_def
                        );
                }
                3 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    );
                }
                4 => {
                    println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}n  | n={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], k_def[3], jac_para, jac_def
                    );
                }
                _ => {}
            }
        }

        return;
    }

    println!(
        "Integrating {} with {} samples and deformation '{}'",
        settings.general.topology, settings.integrator.n_max, settings.general.deformation_strategy
    );

    match topo.analytical_result_real {
        Some(_) => println!(
            "Analytic result: {:e}",
            num::Complex::<f64>::new(
                topo.analytical_result_real.unwrap(),
                topo.analytical_result_imag.unwrap()
            )
        ),
        _ => println!("Analytic result not available."),
    }

    let cuba_result = match settings.integrator.integrator {
        Integrator::Vegas => ci.vegas(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            CubaVerbosity::Progress,
            0,
            UserData {
                integrand: (0..=cores)
                    .map(|i| Integrand::new(topo, settings.clone(), i))
                    .collect(),
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
        Integrator::Suave => ci.suave(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_new,
            settings.integrator.n_min,
            settings.integrator.flatness,
            CubaVerbosity::Progress,
            UserData {
                integrand: (0..=cores)
                    .map(|i| Integrand::new(topo, settings.clone(), i))
                    .collect(),
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
        Integrator::Cuhre => ci.cuhre(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            CubaVerbosity::Progress,
            UserData {
                integrand: (0..=cores)
                    .map(|i| Integrand::new(topo, settings.clone(), i))
                    .collect(),
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
    };
    println!("{:#?}", cuba_result);
    match topo.analytical_result_real {
        Some(_) => println!(
            "Analytic result: {:e}",
            num::Complex::<f64>::new(
                topo.analytical_result_real.unwrap(),
                topo.analytical_result_imag.unwrap()
            )
        ),
        _ => println!("Analytic result not available."),
    }
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
