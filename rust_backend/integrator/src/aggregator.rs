use std::collections::HashMap;
use std::fs::File;

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};
use evaluator::{integrand, Evaluator, Topology, UserData};
use {REGION_ALL, REGION_EXT, REGION_INT};

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct Settings {
    pub active_topology: String,
    pub on_shell_flag: usize,
    pub mu_sq: f64,
    pub dual: bool,
    pub multi_channeling: bool,
    pub alpha: Option<f64>,
    pub regions: bool,
    pub min_eval: usize,
    pub nstart: usize,
    pub nincrease: usize,
    pub max_eval: usize,
    pub cores: usize,
    pub seed: i32,
    pub param_mode: String,
    pub topologies: HashMap<String, Topology>,
    pub survey_n_iterations: usize,
    pub survey_n_points: usize,
    pub refine_n_runs: usize,
    pub refine_n_points: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            active_topology: "box".to_owned(),
            on_shell_flag: 0,
            mu_sq: -1e9,
            dual: true,
            multi_channeling: false,
            regions: false,
            alpha: None,
            min_eval: 10,
            nstart: 1000,
            nincrease: 500,
            max_eval: 1000000,
            cores: 4,
            seed: 0,
            param_mode: "log".to_owned(),
            topologies: HashMap::new(),
            survey_n_iterations: 5,
            survey_n_points: 10000000,
            refine_n_points: 10000000,
            refine_n_runs: 20,
        }
    }
}

impl Settings {
    pub fn from_file(filename: &str) -> Settings {
        let f = File::open(filename).unwrap();
        serde_yaml::from_reader(f).unwrap()
    }
}

pub struct Aggregator {
    n_loops: usize,
    do_regions: bool,
    do_multichanneling: bool,
    evaluator: Evaluator,
    settings: Settings,
}

impl Aggregator {
    pub fn new(settings: Settings) -> Aggregator {
        Aggregator {
            n_loops: settings.topologies[&settings.active_topology].loops,
            do_regions: settings.regions,
            do_multichanneling: settings.multi_channeling,
            evaluator: settings.topologies[&settings.active_topology].build_evaluator(
                settings.mu_sq,
                settings.dual,
                settings.alpha,
            ),
            settings,
        }
    }

    fn integrate(&mut self, region: usize, channel: usize) -> CubaResult {
        println!(
            "Integration of {} with region {} and channel {}",
            self.settings.active_topology, region, channel
        );

        // create an evaluator
        let mut eval = self.evaluator.clone();

        eval.parameterizer.set_region(region);
        eval.parameterizer.set_channel(channel);
        eval.deformer.set_region(region);
        eval.integrand.set_region(region);
        eval.integrand.set_channel(channel);

        if region == REGION_EXT {
            eval.parameterizer.set_mode("weinzierl").unwrap();
        } else {
            if channel == 0 {
                eval.parameterizer
                    .set_mode(&self.settings.param_mode)
                    .unwrap();
            } else {
                eval.parameterizer.set_mode("weinzierl").unwrap();
            }
        }

        let mut ci = CubaIntegrator::new(integrand);
        ci.set_epsabs(0.)
            .set_mineval(self.settings.min_eval as i64)
            .set_nstart(self.settings.nstart as i64)
            .set_nincrease(self.settings.nincrease as i64)
            .set_maxeval(self.settings.max_eval as i64)
            .set_seed(self.settings.seed as i32)
            .set_cores(self.settings.cores, 1000)
            .set_use_only_last_sample(false)
            .set_save_state_file("vegas_saved_state.dat".to_string())
            .set_keep_state_file(false)
            .set_reset_vegas_integrator(false);

        let final_result = {
            if self.settings.refine_n_runs > 0 {
                // Assign cuba flags according to this chosen survey+refine strategy
                ci.set_save_state_file("vegas_saved_state.dat".to_string())
                    .set_keep_state_file(true)
                    .set_reset_vegas_integrator(true);

                // First perform the survey, making sure that there is no previously existing file
                // "vegas_saved_state" or "vegas_survey.dat"
                let _ = std::fs::remove_file("vegas_saved_state.dat");
                let _ = std::fs::remove_file("vegas_survey.dat");

                // Keep track of the total number of failed points
                let mut total_fails = 0;

                // Now we can run the survey, using the specified number of iterations and sample sampel points
                println!(
                    ">>> Now running Vegas survey with {} iterations of {} points:",
                    self.settings.survey_n_iterations, self.settings.survey_n_points
                );

                ci.set_nstart(self.settings.survey_n_points as i64)
                    .set_nincrease(0 as i64)
                    .set_maxeval(
                        (self.settings.survey_n_iterations * self.settings.survey_n_points) as i64,
                    );
                let survey_result = ci.vegas(
                    4 * self.n_loops,
                    1,
                    CubaVerbosity::Progress,
                    1, // Save grid in slot 1
                    UserData {
                        evaluator: vec![eval; self.settings.cores + 1],
                        running_max: 0f64,
                    },
                );
                total_fails += survey_result.fail;
                println!(">>> Survey result : {:#?}", survey_result);

                // Now move the saved state created to file 'vegas_survey.dat'
                let _ = std::fs::rename("vegas_saved_state.dat", "vegas_survey.dat");
                println!("Survey grid files saved in file 'vegas_survey.dat'.");

                let mut vegas_central_values: Vec<Vec<f64>> =
                    Vec::with_capacity(self.settings.refine_n_runs);
                let mut vegas_errors: Vec<Vec<f64>> =
                    Vec::with_capacity(self.settings.refine_n_runs);

                // We can now start the self.settings.refine_n_runs independent runs
                ci.set_nstart(self.settings.refine_n_points as i64)
                    .set_nincrease(0_i64)
                    .set_maxeval(self.settings.refine_n_points as i64)
                    .set_reset_vegas_integrator(true)
                    .set_use_only_last_sample(true)
                    .set_keep_state_file(false);

                for i_run in 0..self.settings.refine_n_runs {
                    // Udate the seed
                    ci.set_seed((self.settings.seed + (i_run as i32) + 1) as i32);
                    // Reinitialise the saved state to the survey grids
                    let _ = std::fs::copy("vegas_survey.dat", "vegas_saved_state.dat");
                    // create a new evaluator for each refine
                    let mut eval = self.evaluator.clone();
                    println!(
                        ">>> Now running Vegas refine run #{} with {} points:",
                        i_run, self.settings.refine_n_points
                    );
                    let refine_result = ci.vegas(
                        4 * self.n_loops,
                        1,
                        CubaVerbosity::Progress,
                        0, // Do not save grids
                        UserData {
                            evaluator: vec![eval; self.settings.cores + 1],
                            running_max: 0f64,
                        },
                    );
                    total_fails += refine_result.fail;
                    // Make sure to remove any saved state left over
                    let _ = std::fs::remove_file("vegas_saved_state.dat");
                    vegas_central_values.push(refine_result.result.clone());
                    vegas_errors.push(refine_result.error.clone());
                    println!(">>> Refine result #{}: {:#?}", i_run + 1, refine_result);
                }

                // Now combine the result
                let mut combined_central: Vec<f64> = vec![0.; survey_result.result.len()];
                let mut combined_error: Vec<f64> = vec![0.; survey_result.result.len()];
                for (central, error) in vegas_central_values.iter().zip(vegas_errors.iter()) {
                    for i_component in 0..central.len() {
                        combined_central[i_component] +=
                            central[i_component] / error[i_component].powi(2);
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
                    neval: (self.settings.survey_n_points * self.settings.survey_n_iterations
                        + self.settings.refine_n_points * self.settings.refine_n_runs)
                        as i64,
                    fail: total_fails,
                    result: combined_central.clone(),
                    error: combined_error.clone(),
                    prob: vec![-1.; combined_error.len()],
                }
            } else {
                ci.vegas(
                    4 * self.n_loops,
                    1,
                    CubaVerbosity::Progress,
                    0, // Not saving the grid
                    UserData {
                        evaluator: vec![eval; self.settings.cores + 1],
                        running_max: 0f64,
                    },
                )
            }
        };
        println!(">>> Final integration result : {:#?}", final_result);
        final_result
    }

    /// Aggregate the results from several integrations over regions and channels.
    pub fn aggregate(&mut self) -> (f64, f64) {
        let mut average = 0.;
        let mut error_sq = 0.;

        if self.do_regions {
            // do the external region
            let ext = self.integrate(REGION_EXT, 0);
            average = ext.result[0];
            error_sq = ext.error[0] * ext.error[0];

            if self.do_multichanneling {
                // TODO: the topology should define how many channels we have
                // ask from the evaluator
                for i in 1..4 {
                    let r = self.integrate(REGION_INT, i);
                    average += r.result[0];
                    error_sq += r.error[0] * r.error[0];
                }
            } else {
                let int = self.integrate(REGION_INT, 0);
                average += int.result[0];
                error_sq += int.error[0] * int.error[0];
            }
        } else {
            if self.do_multichanneling {
                for i in 1..4 {
                    let r = self.integrate(REGION_ALL, i);
                    average += r.result[0];
                    error_sq += r.error[0] * r.error[0];
                }
            } else {
                let r = self.integrate(REGION_ALL, 0);

                average = r.result[0];
                error_sq = r.error[0] * r.error[0];
            }
        }

        (average, error_sq.sqrt())
    }
}
