extern crate arrayvec;
extern crate clap;
extern crate cuba;
extern crate dual_num;
#[macro_use]
extern crate itertools;
extern crate f128;
extern crate ltd;
extern crate nalgebra as na;
extern crate num;
extern crate num_traits;
extern crate rand;
extern crate serde;
extern crate serde_yaml;
use serde::{Deserialize, Serialize};

use arrayvec::ArrayVec;
use clap::{App, Arg, ArgMatches, SubCommand};
use ltd::topologies::{Cut, CutList};
use ltd::LorentzVector;
use num_traits::real::Real;
use num_traits::ToPrimitive;
use num_traits::{NumCast, Zero};
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

        let _r = topo.evaluate(&x, &mut cache, &None);
    }

    println!("{:#?}", now.elapsed());
}

fn surface_prober<'a>(topo: &Topology, _settings: &Settings, matches: &ArgMatches<'a>) {
    let mut loop_momenta = vec![LorentzVector::<f128::f128>::default(); topo.n_loops];
    let e_cm = topo.e_cm_squared.sqrt();

    let mut k_def: ArrayVec<[LorentzVector<num::Complex<f128::f128>>; ltd::MAX_LOOP]> = (0..topo
        .n_loops)
        .map(|_| LorentzVector::default())
        .collect();
    let mut cache = LTDCache::new(topo);

    let mut rng = rand::thread_rng();

    let mut positive_lm = vec![LorentzVector::<f128::f128>::default(); topo.n_loops];
    let mut negative_lm = vec![LorentzVector::<f128::f128>::default(); topo.n_loops];

    let ids: Vec<_> = match matches.values_of("ids") {
        Some(x) => x.map(|x| usize::from_str(x).unwrap()).collect(),
        None => vec![],
    };

    let samples = usize::from_str(matches.value_of("samples").unwrap()).unwrap();

    for (surf_index, surf) in topo.surfaces.iter().enumerate() {
        if !ids.is_empty() && !ids.contains(&surf_index) {
            continue;
        }

        println!(
            "-> group={}, ellipsoid={}, prop={:?} cut={}, mom_map={:?}, signs={:?}, marker={}, shift={}",
            surf.group, surf.ellipsoid, (surf.onshell_ll_index, surf.onshell_prop_index), CutList(&surf.cut), surf.sig_ll_in_cb,
            surf.signs, surf.delta_sign, surf.shift
        );

        for _ in 0..samples {
            let mut did_break = false;
            let mut bisection_flag = 0;

            // try to bisect
            for _ in 0..100_000_000 {
                // bisect if we have found a positive and negative bound,
                // otherwise randomly sample to find the bounds
                if bisection_flag == 3 {
                    for (lm, pl, nl) in izip!(
                        loop_momenta.iter_mut(),
                        positive_lm.iter(),
                        negative_lm.iter()
                    ) {
                        *lm = (pl + nl) * Into::<f128::f128>::into(0.5);
                    }
                } else {
                    for lm in &mut loop_momenta {
                        for index in 1..4 {
                            lm[index] = ((rng.gen::<f64>() * 2.0 - 1.0) * e_cm * 10.).into();
                        }
                    }
                }

                let mut res = f128::f128::zero();

                let mut cut_index = 0;
                for (cut, ll) in izip!(surf.cut.iter(), topo.loop_lines.iter()) {
                    if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = *cut {
                        if surf.signs[cut_index] == 0 {
                            cut_index += 1;
                            continue;
                        }

                        // construct the cut energy
                        let mut mom = LorentzVector::<f128::f128>::default();
                        for (&cut_sign, lm) in ll.signature.iter().zip(loop_momenta.iter()) {
                            mom += lm * cut_sign.into();
                        }
                        // compute the postive cut energy
                        let q: LorentzVector<f128::f128> = ll.propagators[i].q.cast();
                        let energy = ((mom + q).spatial_squared()
                            + <f128::f128 as NumCast>::from(ll.propagators[i].m_squared).unwrap())
                        .sqrt();

                        res += energy * Into::<f128::f128>::into(surf.signs[cut_index]);

                        cut_index += 1;
                    }
                }

                // now for the surface term
                let mut mom = LorentzVector::<f128::f128>::default();
                let onshell_ll = &topo.loop_lines[surf.onshell_ll_index];
                let onshell_prop = &onshell_ll.propagators[surf.onshell_prop_index];
                for (&surf_sign, lm) in onshell_ll.signature.iter().zip(loop_momenta.iter()) {
                    mom += lm * surf_sign.into();
                }
                let energy = ((mom + onshell_prop.q.cast()).spatial_squared()
                    + <f128::f128 as NumCast>::from(onshell_prop.m_squared).unwrap())
                .sqrt();

                res += energy * Into::<f128::f128>::into(surf.delta_sign);
                res += <f128::f128 as NumCast>::from(surf.shift.t).unwrap();

                // update the bounds
                if res < f128::f128::zero() {
                    for (nl, ll) in negative_lm.iter_mut().zip(loop_momenta.iter()) {
                        *nl = ll.clone();
                    }
                    bisection_flag |= 2;
                }
                if res > f128::f128::zero() {
                    for (pl, ll) in positive_lm.iter_mut().zip(loop_momenta.iter()) {
                        *pl = ll.clone();
                    }
                    bisection_flag |= 1;
                }

                if res.abs() < 1e-14.into() {
                    // set the loop momenta
                    for (k, lm) in k_def.iter_mut().zip(&loop_momenta) {
                        *k = lm.to_complex(true);
                        k.t = num::Complex::default();
                    }

                    // do a full evaluation
                    topo.compute_complex_cut_energies(&k_def, &mut cache);

                    let mut result = num::Complex::<f128::f128>::default();
                    for (cuts, mat) in topo.ltd_cut_options.iter().zip(topo.cb_to_lmb_mat.iter()) {
                        // for each cut coming from the same cut structure
                        for cut in cuts.iter() {
                            result += topo.evaluate_cut(&mut k_def, cut, mat, &mut cache).unwrap();
                        }
                    }

                    println!("Result: {}", result);
                    if result.norm() > Into::<f128::f128>::into(1e8) {
                        println!("^ large result for {:?}", loop_momenta);
                    }
                    did_break = true;
                    break;
                }
            }

            if !did_break {
                println!(
                    "Could not find postive or negative value for surface with cut={} and os={:?}",
                    CutList(&surf.cut),
                    (surf.onshell_ll_index, surf.onshell_prop_index)
                );
            }
        }
    }
}

fn inspect<'a>(topo: &Topology, settings: &mut Settings, matches: &ArgMatches<'a>) {
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
        let (x, k_def, jac_para, jac_def, result) = topo.clone().evaluate(&pt, &mut cache, &None);
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
        let (x, k_def, jac_para, jac_def, result) = topo.clone().evaluate(&pt, &mut cache, &None);
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
            SubCommand::with_name("probe")
                .about("Sample points on hyperboloids and ellipsoids")
                .arg(
                    Arg::with_name("ids")
                        .long("ids")
                        .min_values(1)
                        .help("Only sample these surface ids"),
                )
                .arg(
                    Arg::with_name("samples")
                        .short("s")
                        .default_value("100")
                        .help("Number of samples per surface"),
                ),
        )
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

    if let Some(matches) = matches.subcommand_matches("probe") {
        surface_prober(&topo, &settings, matches);
        return;
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        inspect(&topo, &mut settings, matches);
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
