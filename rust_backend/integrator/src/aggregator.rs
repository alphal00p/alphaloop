use std::collections::HashMap;
use std::fs::File;

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};
use evaluator::{integrand, Evaluator, Topology, UserData};
use {REGION_ALL, REGION_EXT, REGION_INT};

#[derive(Debug, Deserialize)]
#[serde(default)]
pub struct Settings {
    pub active_topology: String,
    pub mu_sq: f64,
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
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            active_topology: "box".to_owned(),
            mu_sq: -1e9,
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
            evaluator: settings.topologies[&settings.active_topology]
                .build_evaluator(settings.mu_sq, settings.alpha),
            settings,
        }
    }

    fn integrate(&mut self, region: usize, channel: usize) -> CubaResult {
        println!("Integration of {} with region {} and channel {}", self.settings.active_topology, region, channel);

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
            .set_mineval(self.settings.min_eval as i32)
            .set_nstart(self.settings.nstart as i32)
            .set_nincrease(self.settings.nincrease as i32)
            .set_maxeval(self.settings.max_eval as i32)
            .set_seed(self.settings.seed as i32)
            .set_cores(self.settings.cores, 1000);

        let r = ci.vegas(
            4 * self.n_loops,
            1,
            CubaVerbosity::Progress,
            0, // don't store the grid for now
            UserData {
                evaluator: vec![eval; self.settings.cores + 1],
                running_max: 0f64,
            },
        );
        println!("{:#?}", r);
        r
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
