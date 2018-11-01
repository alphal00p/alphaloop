use {REGION_ALL, REGION_EXT, REGION_INT};

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};

use evaluator::{Evaluator, UserData, integrand};

#[derive(Debug, Clone)]
pub struct IntegrationSettings {
    pub cores: usize,
    pub samples: usize,
    pub param_mode: String,
}

#[derive(Clone)]
pub struct Aggregator {
    n_loops: usize,
    do_regions: bool,
    do_multichanneling: bool,
    evaluator: Evaluator,
    settings: IntegrationSettings,
}

impl Aggregator {
    pub fn new(
        n_loops: usize,
        do_regions: bool,
        do_multichanneling: bool,
        evaluator: Evaluator,
        settings: IntegrationSettings,
    ) -> Aggregator {
        Aggregator {
            n_loops,
            do_regions,
            do_multichanneling,
            evaluator,
            settings,
        }
    }

    fn integrate(&mut self, region: usize, channel: usize) -> CubaResult {
        println!("Integration region {} and channel {}", region, channel);

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
            .set_mineval(10)
            .set_nstart(10000)
            .set_nincrease(5000)
            .set_maxeval(self.settings.samples as i32)
            .set_pseudo_random(true)
            .set_cores(self.settings.cores, 1000);

        let r = ci.vegas(
            4 * self.n_loops,
            1,
            CubaVerbosity::Progress,
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
