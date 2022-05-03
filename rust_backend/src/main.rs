#[cfg(feature = "use_mpi")]
use ltd::integrand::OwnedIntegrandSample;
#[cfg(feature = "use_mpi")]
use mpi::point_to_point::{Destination, Source};
#[cfg(feature = "use_mpi")]
use mpi::request::WaitGuard;
#[cfg(feature = "use_mpi")]
use mpi::topology::Communicator;

use clap::{App, Arg, ArgMatches, SubCommand};
use color_eyre::Report;
use lorentz_vector::LorentzVector;
use num::Complex;
use num_traits::ToPrimitive;
use rand::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::str::FromStr;
use std::time::Instant;

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};

use ltd::integrand::IntegrandImplementation;
use ltd::integrand::{Integrand, IntegrandSample};
use ltd::squared_topologies::{SquaredTopology, SquaredTopologySet};
use ltd::topologies::Topology;
use ltd::{float, IntegrandType, IntegratedPhase, Integrator, Settings};

use ltd::dashboard::{Dashboard, StatusUpdate, StatusUpdateSender};

use havana::{AverageAndErrorAccumulator, Sample};

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

pub enum Integrands {
    CrossSection(Integrand<SquaredTopologySet>),
}

enum Diagram {
    CrossSection(SquaredTopologySet),
}

struct UserData<'a> {
    n_loops: usize,
    integrand: Vec<Integrands>,
    internal_parallelization: bool,
    integrated_phase: IntegratedPhase,
    #[cfg(feature = "use_mpi")]
    world: &'a mpi::topology::SystemCommunicator,
    #[cfg(not(feature = "use_mpi"))]
    phantom_data: std::marker::PhantomData<&'a usize>,
}

#[cfg(feature = "use_mpi")]
pub fn evaluate_mpi_worker(mut integrand: Integrands, world: &mpi::topology::SystemCommunicator) {
    eprintln!("Slave started: {} rank", world.rank());

    let phase = match &integrand {
        Integrands::CrossSection(t) => t.settings.integrator.integrated_phase,
        Integrands::Topology(t) => t.settings.integrator.integrated_phase,
    };

    let mut f = vec![0.; 100];
    let mut iter = 0;
    loop {
        // TODO: check status tag to break the loop
        let (msg, _status) = world.any_process().receive_vec();

        let msg: Vec<OwnedIntegrandSample> = bincode::deserialize(&msg).unwrap();

        if phase != IntegratedPhase::Both {
            f.resize(msg.len(), 0.);
        } else {
            f.resize(msg.len() * 2, 0.);
        };

        // compute all the points
        for (i, owned_sample) in msg.iter().enumerate() {
            let (s, w) = match owned_sample {
                OwnedIntegrandSample::Flat(w, x) => (IntegrandSample::Flat(*w, x), *w),
                OwnedIntegrandSample::Nested(x) => (IntegrandSample::Nested(x), x.get_weight()),
            };

            // TODO: this iter number is not correct, but we are also not using it
            let res = match &mut integrand {
                Integrands::Topology(t) => t.evaluate(s, w, iter),
                Integrands::CrossSection(t) => t.evaluate(s, w, iter),
            };

            if res.is_finite() {
                match phase {
                    IntegratedPhase::Real => {
                        f[i] = res.re.to_f64().unwrap();
                    }
                    IntegratedPhase::Imag => {
                        f[i] = res.im.to_f64().unwrap();
                    }
                    IntegratedPhase::Both => {
                        f[i * 2] = res.re.to_f64().unwrap();
                        f[i * 2 + 1] = res.im.to_f64().unwrap();
                    }
                }
            } else {
                if phase != IntegratedPhase::Both {
                    f[i] = 0.
                } else {
                    f[i * 2] = 0.;
                    f[i * 2 + 1] = 0.;
                }
            }
        }

        iter += 1;

        mpi::request::scope(|scope| {
            let _sreq =
                WaitGuard::from(world.process_at_rank(0).immediate_send(scope, f.as_slice()));
        });
    }
}

fn havana_integrate<'a, F>(settings: &Settings, user_data_generator: F) -> CubaResult
where
    F: Fn() -> UserData<'a>,
{
    let mut num_points = 0;

    let mut samples = vec![Sample::new(); settings.integrator.n_start];
    let mut f = vec![0.; settings.integrator.n_start];
    let mut integral = AverageAndErrorAccumulator::new();

    let mut rng = rand::thread_rng();

    let mut user_data = user_data_generator();

    let mut grid = match &user_data.integrand[0] {
        Integrands::CrossSection(t) => t.topologies[0].create_grid(),
    };

    #[cfg(feature = "use_mpi")]
    let mut samples_per_worker = Vec::with_capacity(samples.len());

    let mut iter = 1;
    while num_points < settings.integrator.n_max {
        let cur_points = settings.integrator.n_start + settings.integrator.n_increase * (iter - 1);
        samples.resize(cur_points, Sample::new());
        f.resize(cur_points, 0.);

        for sample in &mut samples[..cur_points] {
            grid.sample(&mut rng, sample);
        }

        let cores = user_data.integrand.len();
        // the number of points per core for all cores but the last, which may have fewer
        let nvec_per_core = (cur_points - 1) / cores + 1;

        #[cfg(not(feature = "use_mpi"))]
        user_data.integrand[..cores]
            .into_par_iter()
            .zip(f.par_chunks_mut(nvec_per_core))
            .zip(samples.par_chunks(nvec_per_core))
            .for_each(|((integrand_f, ff), xi)| {
                for (ffi, s) in ff.iter_mut().zip(xi.iter()) {
                    let fc = match integrand_f {
                        Integrands::CrossSection(t) => {
                            t.evaluate(IntegrandSample::Nested(s), s.get_weight(), iter)
                        }
                    };

                    let f = match settings.integrator.integrated_phase {
                        IntegratedPhase::Real => fc.re,
                        IntegratedPhase::Imag => fc.im,
                        IntegratedPhase::Both => unimplemented!(),
                    };

                    *ffi = f;
                }
            });

        #[cfg(feature = "use_mpi")]
        {
            if user_data.world.size() == 0 {
                panic!("No workers registered");
            }

            let workers = (user_data.world.size() - 1) as usize;
            let n_samples_per_worker = (samples[..cur_points].len() - 1) / workers + 1; // the last worker may have less

            samples_per_worker.clear();
            for x in samples[..cur_points].chunks(n_samples_per_worker) {
                let v: Vec<_> = x
                    .iter()
                    .map(|s| OwnedIntegrandSample::Nested(s.clone()))
                    .collect();
                samples_per_worker.push(bincode::serialize(&v).unwrap());
            }

            mpi::request::scope(|scope| {
                for (i, x) in samples_per_worker.iter().enumerate() {
                    let _sreq = WaitGuard::from(
                        user_data
                            .world
                            .process_at_rank(i as i32 + 1)
                            .immediate_send(scope, x.as_slice()),
                    );
                }
            });

            for _ in 0..workers {
                let (msg, status) = user_data.world.any_process().receive_vec::<f64>();
                let worker_id = status.source_rank() as usize - 1; // TODO: make sure it's positive

                if user_data.integrated_phase == IntegratedPhase::Both {
                    f[worker_id * n_samples_per_worker * 2
                        ..worker_id * n_samples_per_worker * 2 + msg.len()]
                        .copy_from_slice(&msg);
                } else {
                    f[worker_id * n_samples_per_worker
                        ..worker_id * n_samples_per_worker + msg.len()]
                        .copy_from_slice(&msg);
                }
            }
        }

        for (s, f) in samples[..cur_points].iter().zip(&f[..cur_points]) {
            grid.add_training_sample(s, *f).unwrap();
            integral.add_sample(*f * s.get_weight(), Some(s));
        }

        grid.update(
            settings.integrator.learning_rate,
            settings.integrator.n_bins,
            settings.integrator.train_on_avg,
        );
        integral.update_iter();

        // now merge all statistics and observables into the first
        let (first, others) = user_data.integrand[..cores].split_at_mut(1);
        for i in others {
            match (&mut first[0], i) {
                (Integrands::CrossSection(i1), Integrands::CrossSection(i2)) => {
                    i1.merge_statistics(i2)
                }
            }
        }

        #[cfg(not(feature = "use_mpi"))]
        match &mut user_data.integrand[0] {
            Integrands::CrossSection(i) => i.broadcast_statistics(),
        }

        #[cfg(feature = "use_mpi")]
        println!(
            "Iteration {}:\n{} {:.2} χ²",
            integral.cur_iter,
            ltd::utils::format_uncertainty(integral.avg, integral.err),
            integral.chi_sq / integral.cur_iter as f64,
        );

        if let havana::Grid::DiscreteGrid(g) = &grid {
            g.discrete_dimensions[0].plot("grid_disc.svg").unwrap();
        }

        iter += 1;
        num_points += cur_points;
    }

    // TODO: support multiple dimensions in the ouput
    CubaResult {
        neval: settings.integrator.n_max as i64,
        fail: 0,
        result: vec![integral.avg],
        error: vec![integral.err],
        prob: vec![integral.chi_sq],
    }
}

/// Integrate with Vegas, optionally using survey and refining rounds
fn vegas_integrate<'a, F>(
    name: &str,
    n_loops: usize,
    settings: &Settings,
    mut ci: CubaIntegrator,
    user_data_generator: F,
) -> CubaResult
where
    F: Fn() -> UserData<'a>,
{
    let state_filename = if let Some(ref prefix) = settings.integrator.state_filename_prefix {
        prefix.clone() + &name.clone() + "_state.dat"
    } else {
        name.to_owned() + "_state.dat"
    };
    let survey_filename = if let Some(ref prefix) = settings.integrator.state_filename_prefix {
        prefix.clone() + &name.clone() + "_survey.dat"
    } else {
        name.to_owned() + "_survey.dat"
    };

    ci.set_use_only_last_sample(false)
        .set_keep_state_file(false)
        .set_reset_vegas_integrator(false);

    if settings.integrator.load_from_state_file {
        ci.set_save_state_file(state_filename.clone());
    }

    if settings.integrator.refine_n_runs > 0 {
        // Assign cuba flags according to this chosen survey+refine strategy
        ci.set_save_state_file(state_filename.clone())
            .set_keep_state_file(true)
            .set_reset_vegas_integrator(true);

        // First perform the survey, making sure that there is no previously existing file
        let _ = std::fs::remove_file(&state_filename);
        let _ = std::fs::remove_file(&survey_filename);

        // Keep track of the total number of failed points
        let mut total_fails = 0;

        // Now we can run the survey, using the specified number of iterations and sample sampel points
        println!(
            ">>> Now running Vegas survey with {} iterations of {} points:",
            settings.integrator.survey_n_iterations, settings.integrator.survey_n_points
        );

        ci.set_nstart(settings.integrator.survey_n_points as i64)
            .set_nincrease(0 as i64)
            .set_maxeval(
                (settings.integrator.survey_n_iterations * settings.integrator.survey_n_points)
                    as i64,
            )
            .set_reset_vegas_integrator(settings.integrator.reset_vegas_integrator)
            .set_use_only_last_sample(settings.integrator.use_only_last_sample)
            .set_keep_state_file(settings.integrator.keep_state_file);
        let survey_result = ci.vegas(
            3 * n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_vec,
            if settings.integrator.dashboard {
                CubaVerbosity::Silent
            } else {
                CubaVerbosity::Progress
            },
            1, // Save grid in slot 1
            vegas_integrand,
            user_data_generator(),
        );

        total_fails += survey_result.fail;
        println!(">>> Survey result : {:#?}", survey_result);

        // Now move the saved state
        let _ = std::fs::rename(&state_filename, &survey_filename);
        println!("Survey grid files saved in file {}", survey_filename);

        let mut vegas_central_values: Vec<Vec<f64>> =
            Vec::with_capacity(settings.integrator.refine_n_runs);
        let mut vegas_errors: Vec<Vec<f64>> = Vec::with_capacity(settings.integrator.refine_n_runs);

        // We can now start the self.settings.refine_n_runs independent runs
        ci.set_nstart(settings.integrator.refine_n_points as i64)
            .set_nincrease(0_i64)
            .set_maxeval(settings.integrator.refine_n_points as i64)
            .set_reset_vegas_integrator(true)
            .set_use_only_last_sample(true)
            .set_keep_state_file(false);

        for i_run in 0..settings.integrator.refine_n_runs {
            // Udate the seed
            ci.set_seed((settings.integrator.seed + (i_run as i32) + 1) as i32);
            // Reinitialise the saved state to the survey grids
            let _ = std::fs::copy(&survey_filename, &state_filename);
            println!(
                ">>> Now running Vegas refine run #{} with {} points:",
                i_run, settings.integrator.refine_n_points
            );
            let refine_result = ci.vegas(
                3 * n_loops,
                if settings.integrator.integrated_phase == IntegratedPhase::Both {
                    2
                } else {
                    1
                },
                settings.integrator.n_vec,
                if settings.integrator.dashboard {
                    CubaVerbosity::Silent
                } else {
                    CubaVerbosity::Progress
                },
                0,
                vegas_integrand,
                user_data_generator(),
            );

            total_fails += refine_result.fail;
            // Make sure to remove any saved state left over
            let _ = std::fs::remove_file(&state_filename);
            vegas_central_values.push(refine_result.result.clone());
            vegas_errors.push(refine_result.error.clone());
            println!(">>> Refine result #{}: {:#?}", i_run + 1, refine_result);
        }

        // Now combine the result
        let mut combined_central: Vec<f64> = vec![0.; survey_result.result.len()];
        let mut combined_error: Vec<f64> = vec![0.; survey_result.result.len()];
        for (central, error) in vegas_central_values.iter().zip(vegas_errors.iter()) {
            for i_component in 0..central.len() {
                combined_central[i_component] += central[i_component] / error[i_component].powi(2);
                combined_error[i_component] += 1. / error[i_component].powi(2);
            }
        }
        for i_component in 0..combined_central.len() {
            combined_central[i_component] /= combined_error[i_component];
            combined_error[i_component] = 1. / combined_error[i_component].sqrt();
        }
        // Now return the corresponding CubaResult
        // TODO: For now only the first integrand component is handled
        // The corresponding probability is also not aggregated
        CubaResult {
            neval: (settings.integrator.survey_n_points * settings.integrator.survey_n_iterations
                + settings.integrator.refine_n_points * settings.integrator.refine_n_runs)
                as i64,
            fail: total_fails,
            result: combined_central.clone(),
            error: combined_error.clone(),
            prob: vec![-1.; combined_error.len()],
        }
    } else {
        ci.vegas(
            3 * n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_vec,
            if settings.integrator.dashboard {
                CubaVerbosity::Silent
            } else {
                CubaVerbosity::Progress
            },
            0,
            vegas_integrand,
            user_data_generator(),
        )
    }
}

#[inline(always)]
#[allow(unused_variables)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    nvec: usize,
    core: i32,
    weight: &[f64],
    iter: usize,
) -> Result<(), &'static str> {
    let integrator_settings = match &user_data.integrand[0] {
        Integrands::CrossSection(t) => &t.settings.integrator,
    };

    let n_start = integrator_settings.n_start;
    let n_increase = integrator_settings.n_increase;

    #[cfg(not(feature = "use_mpi"))]
    if user_data.internal_parallelization {
        let cores = user_data.integrand.len();
        // the number of points per core for all cores but the last, which may have fewer
        let nvec_per_core = (nvec - 1) / cores + 1;
        let loops = user_data.n_loops;
        let phase = user_data.integrated_phase;
        let f_len = if IntegratedPhase::Both == user_data.integrated_phase {
            2
        } else {
            1
        };

        user_data.integrand[..cores]
            .into_par_iter()
            .zip(f.par_chunks_mut(f_len * nvec_per_core))
            .zip(x.par_chunks(3 * loops * nvec_per_core))
            .enumerate()
            .for_each(|(i, ((integrand_f, ff), xi))| {
                for (ii, (y, fff)) in xi.chunks(3 * loops).zip(ff.chunks_mut(f_len)).enumerate() {
                    let w = if weight.len() > 0 {
                        // NOTE: only correct for Vegas
                        weight[i * nvec_per_core + ii] * (n_start + (iter - 1) * n_increase) as f64
                    } else {
                        1.
                    };

                    let res = match integrand_f {
                        Integrands::CrossSection(c) => {
                            c.evaluate(IntegrandSample::Flat(w, y), w, iter)
                        }
                    };

                    if res.is_finite() {
                        match phase {
                            IntegratedPhase::Real => {
                                fff[0] = res.re.to_f64().unwrap();
                            }
                            IntegratedPhase::Imag => {
                                fff[0] = res.im.to_f64().unwrap();
                            }
                            IntegratedPhase::Both => {
                                fff[0] = res.re.to_f64().unwrap();
                                fff[1] = res.im.to_f64().unwrap();
                            }
                        }
                    } else {
                        if f_len == 1 {
                            fff[0] = 0.
                        } else {
                            fff[0] = 0.;
                            fff[1] = 0.;
                        }
                    }
                }
            });

        // now merge all statistics and observables into the first
        let (first, others) = user_data.integrand[..cores].split_at_mut(1);
        for i in others {
            match (&mut first[0], i) {
                (Integrands::CrossSection(i1), Integrands::CrossSection(i2)) => {
                    i1.merge_statistics(i2)
                }
            }
        }

        match &mut user_data.integrand[0] {
            Integrands::CrossSection(i) => i.broadcast_statistics(),
        }
    } else {
        for (i, y) in x.chunks(3 * user_data.n_loops).enumerate() {
            let w = if weight.len() > 0 {
                // NOTE: only correct for Vegas
                weight[i] * (n_start + (iter - 1) * n_increase) as f64
            } else {
                1.
            };

            let res = match &mut user_data.integrand[(core + 1) as usize] {
                Integrands::CrossSection(c) => c.evaluate(IntegrandSample::Flat(w, y), w, iter),
            };

            if res.is_finite() {
                match user_data.integrated_phase {
                    IntegratedPhase::Real => {
                        f[i] = res.re.to_f64().unwrap();
                    }
                    IntegratedPhase::Imag => {
                        f[i] = res.im.to_f64().unwrap();
                    }
                    IntegratedPhase::Both => {
                        f[i * 2] = res.re.to_f64().unwrap();
                        f[i * 2 + 1] = res.im.to_f64().unwrap();
                    }
                }
            } else {
                if user_data.integrated_phase != IntegratedPhase::Both {
                    f[i] = 0.
                } else {
                    f[i * 2] = 0.;
                    f[i * 2 + 1] = 0.;
                }
            }
        }

        match &mut user_data.integrand[0] {
            Integrands::CrossSection(i) => i.broadcast_statistics(),
        }
    }

    #[cfg(feature = "use_mpi")]
    {
        if user_data.world.size() == 0 {
            panic!("No workers registered");
        }

        let workers = (user_data.world.size() - 1) as usize;
        let samples_per_worker = (weight.len() - 1) / workers + 1; // the last worker may have less

        let mut samples = Vec::with_capacity(weight.len());
        for (w, xx) in weight.iter().zip(x.chunks(3 * user_data.n_loops)) {
            let w = w
                * (integrator_settings.n_start + (iter - 1) * integrator_settings.n_increase)
                    as f64;

            samples.push(OwnedIntegrandSample::Flat(w, xx.to_vec()))
        }

        let mut d = Vec::with_capacity(samples.len());
        for x in samples.chunks(samples_per_worker) {
            d.push(bincode::serialize(&x).unwrap());
        }

        mpi::request::scope(|scope| {
            for (i, x) in d.iter().enumerate() {
                let _sreq = WaitGuard::from(
                    user_data
                        .world
                        .process_at_rank(i as i32 + 1)
                        .immediate_send(scope, x.as_slice()),
                );
            }
        });

        for _ in 0..workers {
            let (msg, status) = user_data.world.any_process().receive_vec::<f64>();
            let worker_id = status.source_rank() as usize - 1;

            if user_data.integrated_phase == IntegratedPhase::Both {
                f[worker_id * samples_per_worker * 2
                    ..worker_id * samples_per_worker * 2 + msg.len()]
                    .copy_from_slice(&msg);
            } else {
                f[worker_id * samples_per_worker..worker_id * samples_per_worker + msg.len()]
                    .copy_from_slice(&msg);
            }
        }
    }

    Ok(())
}

#[inline(always)]
fn vegas_integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    nvec: usize,
    core: i32,
    weight: &[f64],
    iter: usize,
) -> Result<(), &'static str> {
    integrand(x, f, user_data, nvec, core, weight, iter)
}

#[inline(always)]
fn cuhre_integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    nvec: usize,
    core: i32,
) -> Result<(), &'static str> {
    integrand(x, f, user_data, nvec, core, &[], 1)
}

#[inline(always)]
fn suave_integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    nvec: usize,
    core: i32,
    weight: &[f64],
    iter: usize,
) -> Result<(), &'static str> {
    integrand(x, f, user_data, nvec, core, weight, iter)
}

#[inline(always)]
fn divonne_integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    nvec: usize,
    core: i32,
    _phase: usize,
) -> Result<(), &'static str> {
    integrand(x, f, user_data, nvec, core, &[], 1)
}

fn bench(
    diagram: &Diagram,
    status_update_sender: StatusUpdateSender,
    settings: &Settings,
    num_samples: usize,
) {
    let mut integrand = match diagram {
        Diagram::CrossSection(sqt) => Integrands::CrossSection(Integrand::new(
            sqt.get_maximum_loop_count(),
            sqt.clone(),
            settings.clone(),
            true,
            status_update_sender,
            1,
            None,
        )),
    };

    let n_loops = match &diagram {
        Diagram::CrossSection(sqt) => sqt.get_maximum_loop_count(),
    };

    let mut x = vec![0.; 3 * n_loops];
    let mut rng: StdRng = SeedableRng::seed_from_u64(100);

    let now = Instant::now();
    for _ in 0..num_samples {
        for xi in x.iter_mut() {
            *xi = rng.gen();
        }

        match &mut integrand {
            Integrands::CrossSection(sqt) => sqt.evaluate(IntegrandSample::Flat(1., &x), 1., 1),
        };
    }

    println!("{:#?}", now.elapsed());
}

fn inspect<'a>(
    diagram: &mut Diagram,
    status_update_sender: StatusUpdateSender,
    settings: &mut Settings,
    matches: &ArgMatches<'a>,
) {
    let (n_loops, e_cm_squared) = match &diagram {
        Diagram::CrossSection(sqt) => (sqt.get_maximum_loop_count(), sqt.e_cm_squared),
    };

    let mut pt: Vec<_> = matches
        .values_of("point")
        .unwrap()
        .map(|x| f64::from_str(x.trim_end_matches(',')).unwrap())
        .collect();
    if pt.len() != 3 * n_loops {
        panic!(
            "Dimension of the input point is incorrect. It should be {} but is {}.",
            n_loops * 3,
            pt.len()
        );
    }

    if matches.is_present("momentum_space") {
        // map the point back from momentum-space to the unit hypercube
        for (i, x) in pt.chunks_exact_mut(3).enumerate() {
            let r = Topology::inv_parametrize::<f128::f128>(
                &LorentzVector::from_args(0., x[0], x[1], x[2]).cast(),
                e_cm_squared,
                i,
                &settings,
            );
            x[0] = f128::f128::to_f64(&r.0[0]).unwrap();
            x[1] = f128::f128::to_f64(&r.0[1]).unwrap();
            x[2] = f128::f128::to_f64(&r.0[2]).unwrap();
        }
    }

    if matches.is_present("full_integrand") {
        let result = match diagram {
            Diagram::CrossSection(sqt) => Integrand::new(
                sqt.get_maximum_loop_count(),
                sqt.clone(),
                settings.clone(),
                true,
                status_update_sender,
                1,
                None,
            )
            .evaluate(IntegrandSample::Flat(1., &pt), 1., 1),
        };
        println!("result={:e}\n  | x={:?}\n", result, pt);
        return;
    }

    // TODO: prevent code repetition
    if matches.is_present("use_f128") {
        let result = match diagram {
            Diagram::CrossSection(sqt) => {
                if matches.is_present("use_ltd") {
                    for t in &mut sqt.topologies {
                        t.settings.cross_section.integrand_type = IntegrandType::LTD;
                    }
                }

                let mut cache = sqt.create_caches();
                sqt.clone()
                    .evaluate::<f128::f128>(IntegrandSample::Flat(1., &pt), &mut cache, None)
            }
        };

        println!("result={:e}\n  | x={:?}\n", result, pt);
    } else {
        let result = match diagram {
            Diagram::CrossSection(sqt) => {
                if matches.is_present("use_ltd") {
                    for t in &mut sqt.topologies {
                        t.settings.cross_section.integrand_type = IntegrandType::LTD;
                    }
                }

                let mut cache = sqt.create_caches();
                sqt.clone()
                    .evaluate::<float>(IntegrandSample::Flat(1., &pt), &mut cache, None)
            }
        };
        println!("result={:e}\n  | x={:?}\n", result, pt);
    }
}

fn main() -> Result<(), Report> {
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
            Arg::with_name("debug")
                .long("debug")
                .value_name("LEVEL")
                .help("Set the debug level. Higher means more verbose."),
        )
        .arg(
            Arg::with_name("samples")
                .short("s")
                .long("samples")
                .value_name("SAMPLES")
                .help("Number of samples per integration"),
        )
        .arg(
            Arg::with_name("n_start")
                .long("n_start")
                .value_name("N_START")
                .help("Number of starting samples for Vegas"),
        )
        .arg(
            Arg::with_name("n_increase")
                .long("n_increase")
                .value_name("N_INCREASE")
                .help("Number of increase samples for Vegas"),
        )
        .arg(
            Arg::with_name("integrator")
                .long("integrator")
                .value_name("INTEGRATOR")
                .help("Select the integrator (Vegas, Cuhre, Suave, Divonne)"),
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
            Arg::with_name("seed")
                .long("seed")
                .value_name("SEED")
                .help("Specify the integration seed"),
        )
        .arg(
            Arg::with_name("target")
                .long("target")
                .multiple(true)
                .allow_hyphen_values(true)
                .value_name("TARGET")
                .help("Specify the integration target a <real> <imag>"),
        )
        .arg(
            Arg::with_name("topology")
                .short("t")
                .long("topology")
                .value_name("TOPOLOGY")
                .help("Set the active topology"),
        )
        .arg(
            Arg::with_name("amplitude")
                .short("a")
                .long("amplitude")
                .value_name("AMPLITUDE")
                .help("Set the active amplitude"),
        )
        .arg(
            Arg::with_name("state_filename_prefix")
                .long("state_filename_prefix")
                .value_name("STATE_FILENAME")
                .help("Set the prefix to apply to vegas grid file"),
        )
        .arg(
            Arg::with_name("log_file_prefix")
                .long("log_file_prefix")
                .value_name("LOG_FILE_PREFIX")
                .help("Set the prefix to apply to the integration statistics log"),
        )
        .arg(
            Arg::with_name("res_file_prefix")
                .long("res_file_prefix")
                .value_name("RES_FILE_PREFIX")
                .help("Set the prefix to apply to the result file"),
        )
        .arg(
            Arg::with_name("cross_section")
                .long("cross_section")
                .value_name("SQUARED_TOPOLOGY_FILE")
                .help("Set the squared topology file"),
        )
        .arg(
            Arg::with_name("cross_section_set")
                .long("cross_section_set")
                .value_name("SQUARED_TOPOLOGY_SET_FILE")
                .help("Set the squared topology set file"),
        )
        .subcommand(
            SubCommand::with_name("integrated_ct")
                .about("Gives the integrated CT for the selected amplitude"),
        )
        .subcommand(
            SubCommand::with_name("bench").about("Run a benchmark").arg(
                Arg::with_name("samples")
                    .required(true)
                    .long("samples")
                    .short("s")
                    .value_name("SAMPLES")
                    .help("Number of samples for benchmark"),
            ),
        )
        .subcommand(
            SubCommand::with_name("inspect")
                .about("Inspect a single input point")
                .arg(
                    Arg::with_name("point")
                        .short("p")
                        .required(true)
                        .min_values(3)
                        .allow_hyphen_values(true),
                )
                .arg(
                    Arg::with_name("use_f128")
                        .short("f128")
                        .long("use_f128")
                        .help("Use f128 evaluation"),
                )
                .arg(
                    Arg::with_name("use_ltd")
                        .short("ltd")
                        .long("use_ltd")
                        .help("Use LTD instead of cLTD for the evaluation"),
                )
                .arg(
                    Arg::with_name("momentum_space")
                        .short("m")
                        .long("momentum_space")
                        .help("Set if the point is specified in momentum space"),
                )
                .arg(
                    Arg::with_name("full_integrand")
                        .long("full_integrand")
                        .help("Evaluate the integrand and possibly its rotated vesion"),
                ),
        )
        .get_matches();

    let mut settings = Settings::from_file(matches.value_of("config").unwrap())?;

    let mut cores = 1;
    if let Some(x) = matches.value_of("cores") {
        cores = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("debug") {
        settings.general.debug = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("seed") {
        settings.integrator.seed = i32::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("samples") {
        settings.integrator.n_max = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("n_start") {
        settings.integrator.n_start = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("n_increase") {
        settings.integrator.n_increase = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("integrator") {
        settings.integrator.integrator = serde_yaml::from_str::<Integrator>(x).unwrap();
    }

    if let Some(x) = matches.value_of("topology") {
        settings.general.topology = x.to_owned();
    }

    if let Some(x) = matches.value_of("state_filename_prefix") {
        settings.integrator.state_filename_prefix = Some(x.to_owned());
    }

    if let Some(x) = matches.value_of("res_file_prefix") {
        settings.general.res_file_prefix = x.to_owned();
    }

    if let Some(x) = matches.value_of("deformation_strategy") {
        settings.general.deformation_strategy = x.into();
    }

    let mut target = None;
    if let Some(t) = matches.values_of("target") {
        let tt: Vec<_> = t
            .map(|x| f64::from_str(x.trim_end_matches(',')).unwrap())
            .collect();
        if tt.len() != 2 {
            panic!("Expected two numbers for target");
        }
        target = Some(Complex::new(tt[0], tt[1]));
    }

    if !settings.observables.is_empty()
        && ((cores > 1 && !settings.integrator.internal_parallelization)
            || (settings.integrator.integrator != Integrator::Vegas
                && settings.integrator.integrator != Integrator::Suave
                && settings.integrator.integrator != Integrator::Havana))
    {
        println!("Removing observable functions because we are not running in single core or because a not supported integrator is selected.");
        settings.observables.clear();
    }

    if settings.integrator.dashboard && !settings.integrator.internal_parallelization && cores > 0
        || (settings.integrator.integrator != Integrator::Vegas
            && settings.integrator.integrator != Integrator::Suave
            && settings.integrator.integrator != Integrator::Havana)
    {
        println!("Cannot use dashboard with Cuba parallelization and cores != 0 or an integrator other than Vegas or Suave");
        settings.integrator.dashboard = false;
    }

    if matches.is_present("bench") || matches.is_present("inspect") {
        settings.integrator.dashboard = false;
    }

    #[cfg(feature = "use_mpi")]
    {
        settings.integrator.dashboard = false;
        settings.integrator.internal_parallelization = false;
    }

    let mut dashboard = Dashboard::new(
        settings.integrator.dashboard,
        settings.integrator.show_plot,
        settings.integrator.quiet_mode,
    );

    let mut diagram = if let Some(cs_opt) = matches.value_of("cross_section") {
        Diagram::CrossSection(SquaredTopologySet::from_one(SquaredTopology::from_file(
            cs_opt, &settings,
        )?))
    } else if let Some(css_opt) = matches.value_of("cross_section_set") {
        Diagram::CrossSection(SquaredTopologySet::from_file(css_opt, &settings)?)
    } else {
        panic!("Specify cross_section or cross_section_set");
    };

    if let Some(matches) = matches.subcommand_matches("bench") {
        bench(
            &diagram,
            dashboard.status_update_sender,
            &settings,
            matches.value_of("samples").unwrap().parse().unwrap(),
        );
        return Ok(());
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        inspect(
            &mut diagram,
            dashboard.status_update_sender,
            &mut settings,
            matches,
        );
        return Ok(());
    }

    match &mut diagram {
        Diagram::CrossSection(sqt) => {
            dashboard
                .status_update_sender
                .send(StatusUpdate::Message(format!(
                    "Integrating {} with {} samples and deformation '{}'",
                    sqt.name, settings.integrator.n_max, settings.general.deformation_strategy
                )))
                .unwrap();
            sqt.print_info(&mut dashboard.status_update_sender);
        }
    }

    let n_loops = match &diagram {
        Diagram::CrossSection(sqt) => sqt.get_maximum_loop_count(),
    };

    #[cfg(feature = "use_mpi")]
    let (_universe, world) = {
        let universe = mpi::initialize().unwrap();
        let world = universe.world();
        let rank = world.rank();

        // if we are not the root, we listen for jobs
        if rank != 0 {
            evaluate_mpi_worker(
                match &diagram {
                    Diagram::CrossSection(sqt) => Integrands::CrossSection(Integrand::new(
                        sqt.get_maximum_loop_count(),
                        sqt.clone(),
                        settings.clone(),
                        true,
                        dashboard.status_update_sender.clone(),
                        1,
                        target,
                    )),
                },
                &world,
            );

            return Ok(());
        } else {
            cores = 0;
        }
        (universe, world)
    };

    let mut ci = CubaIntegrator::new();

    if settings.integrator.internal_parallelization {
        cores = 1.max(cores);

        if settings.integrator.n_vec < cores {
            dashboard
                .status_update_sender
                .send(StatusUpdate::Message(format!(
                "n_vec is less than number of cores: setting it equal to the number of cores. For performance reasons, it should be set a 100 times higher."
            ))).unwrap();
            settings.integrator.n_vec = cores;
        }
    }

    ci.set_mineval(10)
        .set_nstart(settings.integrator.n_start as i64)
        .set_nincrease(settings.integrator.n_increase as i64)
        .set_maxeval(settings.integrator.n_max as i64)
        .set_epsrel(settings.integrator.eps_rel)
        .set_epsabs(settings.integrator.eps_abs)
        .set_border(settings.integrator.border)
        .set_maxpass(settings.integrator.maxpass as i32)
        .set_maxchisq(settings.integrator.maxchisq)
        .set_mindeviation(settings.integrator.mindeviation)
        .set_seed(settings.integrator.seed)
        .set_batch((settings.integrator.n_vec as i64).max(1000))
        .set_cores(
            if settings.integrator.internal_parallelization {
                0
            } else {
                cores
            },
            settings.integrator.n_vec,
        );

    let name = &match &diagram {
        Diagram::CrossSection(sqt) => &sqt.name,
    }
    .clone();

    let num_integrands = if settings.integrator.internal_parallelization {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cores)
            .build_global()
            .unwrap();
        cores
    } else {
        // cuba needs one more core, as the master could theoretically also do some work
        cores + 1
    };

    let user_data_generator = || UserData {
        n_loops,
        integrand: (0..num_integrands)
            .map(|i| match &diagram {
                Diagram::CrossSection(sqt) => Integrands::CrossSection(Integrand::new(
                    sqt.get_maximum_loop_count(),
                    sqt.clone(),
                    settings.clone(),
                    true,
                    dashboard.status_update_sender.clone(),
                    i,
                    target,
                )),
            })
            .collect(),
        internal_parallelization: settings.integrator.internal_parallelization,
        integrated_phase: settings.integrator.integrated_phase,
        #[cfg(feature = "use_mpi")]
        world: &world,
        #[cfg(not(feature = "use_mpi"))]
        phantom_data: std::marker::PhantomData,
    };

    let cuba_result = match settings.integrator.integrator {
        Integrator::Havana => havana_integrate(&settings, user_data_generator),
        Integrator::Vegas => vegas_integrate(name, n_loops, &settings, ci, user_data_generator),
        Integrator::Suave => ci.suave(
            3 * n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_vec,
            settings.integrator.n_new,
            settings.integrator.n_min,
            settings.integrator.flatness,
            if settings.integrator.dashboard {
                CubaVerbosity::Silent
            } else {
                CubaVerbosity::Progress
            },
            suave_integrand,
            user_data_generator(),
        ),
        Integrator::Divonne => ci.divonne(
            3 * n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_vec,
            &[],
            CubaVerbosity::Progress,
            divonne_integrand,
            user_data_generator(),
        ),
        Integrator::Cuhre => ci.cuhre(
            3 * n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_vec,
            CubaVerbosity::Progress,
            cuhre_integrand,
            user_data_generator(),
        ),
    };
    println!("{:#?}", cuba_result);

    let f = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(format!(
            "{}{}.dat",
            settings.general.res_file_prefix,
            settings.general.topology.clone() + "_res"
        ))
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

    Ok(())
}
