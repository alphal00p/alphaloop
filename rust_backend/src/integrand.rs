use crate::dashboard::{StatusUpdate, StatusUpdateSender};
use crate::observables::EventManager;
use crate::squared_topologies::MAX_SG_LOOP;
use crate::{IntegratedPhase, Settings};
use f128::f128;
use havana::{Grid, Sample};
use num::Complex;
use num_traits::{Float, FromPrimitive, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use std::time::Instant;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum OwnedIntegrandSample {
    Flat(f64, Vec<f64>),
    Nested(Sample),
}

#[derive(Debug, Copy, Clone, Serialize)]
pub enum IntegrandSample<'a> {
    Flat(f64, &'a [f64]),
    #[serde(borrow)]
    Nested(&'a Sample),
}

impl<'a> IntegrandSample<'a> {
    pub fn to_flat(&self) -> &'a [f64] {
        match self {
            IntegrandSample::Flat(_, x) => *x,
            IntegrandSample::Nested(x) => match x {
                Sample::ContinuousGrid(_w, v) => &v,
                Sample::DiscreteGrid(_, _, s) => {
                    if let Some(cs) = s {
                        IntegrandSample::Nested(cs.as_ref()).to_flat()
                    } else {
                        unreachable!("No continuous grid found in sample")
                    }
                }
                _ => unimplemented!(),
            },
        }
    }
}

pub trait IntegrandImplementation: Clone {
    type Cache: Default;

    fn create_stability_check(&self, num_checks: usize) -> Vec<Self>;

    fn get_target(&self) -> Option<Complex<f64>>;

    fn set_partial_fractioning(&mut self, enable: bool);

    fn create_grid(&self) -> Grid;

    fn evaluate_f64<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut Self::Cache,
        events: Option<&mut EventManager>,
    ) -> Complex<f64>;

    fn evaluate_f128<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut Self::Cache,
        events: Option<&mut EventManager>,
    ) -> Complex<f128>;

    fn set_precision(&mut self, prec: usize);

    fn create_cache(&self) -> Self::Cache;
}

#[derive(Debug, Clone)]
pub struct IntegrandStatistics {
    pub phase: IntegratedPhase,
    pub target: Option<Complex<f64>>,
    pub running_max_re: (f64, f64, usize),
    pub running_max_im: (f64, f64, usize),
    pub total_samples: usize,
    pub nan_point_count: usize,
    pub unstable_point_count: Vec<usize>,
    pub regular_point_count: usize,
    pub total_sample_time: f64,
    pub n_loops: usize,
    pub running_max_coordinate_re: [f64; 3 * MAX_SG_LOOP],
    pub running_max_coordinate_im: [f64; 3 * MAX_SG_LOOP],
    pub running_max_stability: f64,
    pub integrand_evaluation_timing: u128,
    pub integrand_evaluation_timing_count: usize,
    pub integrand_evaluation_timing_record: bool,
}

impl IntegrandStatistics {
    pub fn new(
        n_loops: usize,
        phase: IntegratedPhase,
        num_stability_levels: usize,
        target: Option<Complex<f64>>,
    ) -> IntegrandStatistics {
        IntegrandStatistics {
            n_loops,
            target,
            phase,
            running_max_re: (0., 0., 0),
            running_max_im: (0., 0., 0),
            total_samples: 0,
            regular_point_count: 0,
            unstable_point_count: vec![0; num_stability_levels],
            nan_point_count: 0,
            total_sample_time: 0.,
            running_max_coordinate_re: [0.; 3 * MAX_SG_LOOP],
            running_max_coordinate_im: [0.; 3 * MAX_SG_LOOP],
            running_max_stability: 0.,
            integrand_evaluation_timing: 0,
            integrand_evaluation_timing_count: 0,
            integrand_evaluation_timing_record: true,
        }
    }

    pub fn merge(&mut self, other: &mut IntegrandStatistics) {
        self.total_samples += other.total_samples;
        self.regular_point_count += other.regular_point_count;

        for (u, ou) in self
            .unstable_point_count
            .iter_mut()
            .zip(&mut other.unstable_point_count)
        {
            *u += *ou;
            *ou = 0;
        }

        self.nan_point_count += other.nan_point_count;
        self.total_sample_time += other.total_sample_time;
        self.integrand_evaluation_timing += other.integrand_evaluation_timing;
        self.integrand_evaluation_timing_count += other.integrand_evaluation_timing_count;

        if self.running_max_re.0 < other.running_max_re.0 {
            self.running_max_re = other.running_max_re;
            self.running_max_coordinate_re[..3 * self.n_loops]
                .copy_from_slice(&other.running_max_coordinate_re[..3 * self.n_loops]);
            if self.phase != IntegratedPhase::Imag {
                self.running_max_stability = other.running_max_stability;
            }
        }
        if self.running_max_im.0 < other.running_max_im.0 {
            self.running_max_im = other.running_max_im;
            self.running_max_coordinate_im[..3 * self.n_loops]
                .copy_from_slice(&other.running_max_coordinate_im[..3 * self.n_loops]);
            if self.phase != IntegratedPhase::Real {
                self.running_max_stability = other.running_max_stability;
            }
        }

        other.total_samples = 0;
        other.regular_point_count = 0;
        other.nan_point_count = 0;
        other.total_sample_time = 0.;
        other.integrand_evaluation_timing = 0;
        other.integrand_evaluation_timing_count = 0;
    }
}

/// A structure that integrates and keeps statistics
pub struct Integrand<I: IntegrandImplementation> {
    pub n_loops: usize,
    pub settings: Settings,
    pub topologies: Vec<I>,
    pub cache: I::Cache,
    pub integrand_statistics: IntegrandStatistics,
    pub cur_iter: usize,
    pub id: usize,
    pub event_manager: EventManager,
    pub status_update_sender: StatusUpdateSender,
}

macro_rules! check_stability_precision {
    ($($name:ident, $eval_fn:ident, $ty:ty),*) => {
        $(
    fn $name(
        &mut self,
        weight: f64,
        x: IntegrandSample<'_>,
        prec: usize,
        relative_precision: f64,
        num_samples: usize,
        escalate_for_large_weight_threshold: f64,
        minimal_precision_to_skip_further_checks: f64,
        cache: &mut I::Cache,
        event_manager: &mut EventManager,
        timing: bool,
    ) -> ($ty, $ty, Complex<$ty>, Complex<$ty>, bool) {
        event_manager.integrand_evaluation_timing = 0;
        event_manager.time_integrand_evaluation = timing;

        for t in &mut self.topologies {
            t.set_precision(prec);
        }

        let result = self.topologies[0].$eval_fn(x, cache, Some(event_manager));

        if timing {
            self.integrand_statistics.integrand_evaluation_timing +=
                event_manager.integrand_evaluation_timing;
            self.integrand_statistics.integrand_evaluation_timing_count += 1;
            event_manager.integrand_evaluation_timing = 0;

            // only record the first 60 timings since it's slow to time
            // note: this only works when nvec > 60
            self.integrand_statistics.integrand_evaluation_timing_record = self.integrand_statistics.integrand_evaluation_timing_count < 60;
        }

        // even though the point may be stable, we may want to escalate if it's large
        let escalate = if escalate_for_large_weight_threshold > 0. {
            match self.settings.integrator.integrated_phase {
                IntegratedPhase::Real => {
                    result.re.abs().to_f64().unwrap() * weight
                        > self.integrand_statistics.running_max_re.0
                            * escalate_for_large_weight_threshold
                }
                IntegratedPhase::Both => {
                    result.re.abs().to_f64().unwrap() * weight
                        > self.integrand_statistics.running_max_re.0
                            * escalate_for_large_weight_threshold
                        || result.im.abs().to_f64().unwrap() * weight
                            > self.integrand_statistics.running_max_im.0
                                * escalate_for_large_weight_threshold
                }
                IntegratedPhase::Imag => {
                    result.im.abs().to_f64().unwrap() * weight
                        > self.integrand_statistics.running_max_im.0
                            * escalate_for_large_weight_threshold
                }
            }
        } else {
            false
        };

        let mut ret = (
            Into::<$ty>::into(relative_precision),
            Into::<$ty>::into(self.settings.general.absolute_precision),
            result,
            result,
            escalate
        );

        if self.topologies.len() == 1
            || num_samples == 1
            || escalate
        {
            return ret;
        }


        // collect smallest result and biggest result
        let mut min = result;
        let mut max = result;

        // do not track events of rotated topologies
        let track_events = event_manager.track_events;
        event_manager.track_events = false;

        for (rot_index, rot_topo) in self.topologies[1..num_samples].iter_mut().enumerate() {
            if self.settings.general.debug > 2 {
                println!("Evaluating integrand with rotated topologies");
            }

            let result_rot = rot_topo.$eval_fn(x, cache, Some(event_manager));

            // compute the number of similar digits
            match self.settings.integrator.integrated_phase {
                IntegratedPhase::Real => {
                    if result_rot.re < min.re {
                        min.re = result_rot.re
                    }
                    if result_rot.re > max.re {
                        max.re = result_rot.re
                    }
                }
                IntegratedPhase::Imag => {
                    if result_rot.im < min.im {
                        min.im = result_rot.im
                    }
                    if result_rot.im > max.im {
                        max.im = result_rot.im
                    }
                }
                IntegratedPhase::Both => {
                    if result_rot.re < min.re {
                        min.re = result_rot.re
                    }
                    if result_rot.re > max.re {
                        max.re = result_rot.re
                    }
                    if result_rot.im < min.im {
                        min.im = result_rot.im
                    }
                    if result_rot.im > max.im {
                        max.im = result_rot.im
                    }
                }
            };

            // determine if the difference is largest in the re or im
            let (num, num_rot) = if max.re - min.re > max.im - min.im {
                (max.re, min.re)
            } else {
                (max.im, min.im)
            };

            let d = if num.is_zero() && num_rot.is_zero() {
                <$ty>::from_usize(100).unwrap()
            } else {
                if num == -num_rot {
                    -<$ty>::from_usize(100).unwrap()
                } else {
                    -Float::abs((num - num_rot) / (num + num_rot)).log10()
                }
            };

            let diff = Float::abs(num - num_rot);

            ret = (d, Float::abs(num - num_rot), min, max, false);

            // if we see the point is unstable, stop further samples
            if !result_rot.is_finite()
                || !min.is_finite()
                || !max.is_finite()
                || d < Into::<$ty>::into(relative_precision)
                || diff > Into::<$ty>::into(self.settings.general.absolute_precision) {
                break;
            }

            // if the first point is very stable, we stop further samples
            if rot_index == 0 && result_rot.is_finite() && d > Into::<$ty>::into(minimal_precision_to_skip_further_checks) {
                break;
            }
        }

        event_manager.track_events = track_events;
        ret
    }

    )*

    }
}

impl<I: IntegrandImplementation> Integrand<I> {
    pub fn new(
        n_loops: usize,
        topology: I,
        settings: Settings,
        track_events: bool,
        status_update_sender: StatusUpdateSender,
        id: usize,
        target: Option<Complex<f64>>,
    ) -> Integrand<I> {
        // create extra topologies with rotated kinematics to check the uncertainty
        let mut topologies = vec![topology.clone()];

        topologies.extend(
            topology.create_stability_check(
                settings
                    .general
                    .stability_checks
                    .iter()
                    .map(|sc| sc.n_samples)
                    .max()
                    .unwrap_or(0),
            ),
        );

        Integrand {
            n_loops,
            topologies,
            cache: topology.create_cache(),
            integrand_statistics: IntegrandStatistics::new(
                n_loops,
                settings.integrator.integrated_phase,
                settings.general.stability_checks.len(),
                target.or(topology.get_target()),
            ),
            settings: settings.clone(),
            id,
            cur_iter: 1, // always start at iter 1
            event_manager: EventManager::new(
                track_events,
                settings.clone(),
                status_update_sender.clone(),
            ),
            status_update_sender,
        }
    }

    check_stability_precision!(check_stability_float, evaluate_f64, f64);
    check_stability_precision!(check_stability_quad, evaluate_f128, f128);

    pub fn merge_statistics(&mut self, other: &mut Integrand<I>) {
        self.event_manager.merge_samples(&mut other.event_manager);

        self.integrand_statistics
            .merge(&mut other.integrand_statistics);
    }

    pub fn broadcast_statistics(&mut self) {
        if self.id == 0 {
            self.status_update_sender
                .send(StatusUpdate::Statistics(self.integrand_statistics.clone()))
                .unwrap();

            self.event_manager.update_live_result();
        }
    }

    /// Evalute a point generated from the Monte Carlo generator with `weight` and current iteration number `iter`.
    pub fn evaluate(&mut self, x: IntegrandSample<'_>, weight: f64, iter: usize) -> Complex<f64> {
        let start_time = Instant::now(); // time the evaluation

        if self.cur_iter != iter {
            // the first integrand accumulates all the results from the others
            if self.id == 0 {
                self.event_manager.update_result();

                self.status_update_sender
                    .send(StatusUpdate::Statistics(self.integrand_statistics.clone()))
                    .unwrap();
            }

            self.cur_iter += 1;
        }

        // Set a global seed equal in all topologies if we need it for rng operations
        if self.settings.general.multi_channeling {
            if let Some(i) = self.settings.general.multi_channeling_channel {
                if i < 0 {
                    unimplemented!("Not implemented for now");
                    /*let global_seed: [u8; 32] = rand::thread_rng().gen();
                    for topo in self.topologies.iter_mut() {
                        topo.global_seed = Some(global_seed);
                    }*/
                }
            }
        }

        let mut event_manager = std::mem::replace(&mut self.event_manager, EventManager::default());
        let mut cache = std::mem::replace(&mut self.cache, I::Cache::default());

        // go through the stability pipeline
        let mut result = Complex::zero();
        let mut stability_level = 0;
        let mut stable = false;
        let mut stable_digits = 0.;
        let mut stability_checks =
            std::mem::replace(&mut self.settings.general.stability_checks, vec![]);
        for (level, stability_check) in stability_checks.iter().enumerate() {
            // clear events when there is instability
            event_manager.clear(false);

            if x.to_flat().chunks(3).any(|v| {
                v[0] < stability_check.accepted_radius_range_in_x_space.0
                    || v[0] > stability_check.accepted_radius_range_in_x_space.1
            }) {
                if self.settings.general.debug > 1 {
                    println!("Skipping stability check at level {} since radius is outside accepted range", level);
                }

                self.integrand_statistics.unstable_point_count[level] += 1;
                continue;
            }

            for t in &mut self.topologies {
                t.set_partial_fractioning(stability_check.use_pf);
            }

            let abs_diff;
            let escalate;
            if stability_check.prec == 16 {
                let (d, diff, min_rot, max_rot, esc) = self.check_stability_float(
                    weight,
                    x,
                    16,
                    stability_check.relative_precision,
                    stability_check.n_samples,
                    if level + 1 == stability_checks.len() {
                        -1.
                    } else {
                        stability_check.escalate_for_large_weight_threshold
                    },
                    stability_check.minimal_precision_to_skip_further_checks,
                    &mut cache,
                    &mut event_manager,
                    level == 0 && self.integrand_statistics.integrand_evaluation_timing_record,
                );

                stable_digits = d;
                abs_diff = diff;
                escalate = esc;
                result = (min_rot + max_rot) / 2.;
            } else {
                let (d, diff, min_rot, max_rot, esc) = self.check_stability_quad(
                    weight,
                    x,
                    stability_check.prec,
                    stability_check.relative_precision,
                    stability_check.n_samples,
                    if level + 1 == stability_checks.len() {
                        -1.
                    } else {
                        stability_check.escalate_for_large_weight_threshold
                    },
                    stability_check.minimal_precision_to_skip_further_checks,
                    &mut cache,
                    &mut event_manager,
                    level == 0,
                );

                let sum = min_rot + max_rot;
                stable_digits = d.to_f64().unwrap();
                abs_diff = diff.to_f64().unwrap();
                escalate = esc;
                result = Complex::new(sum.re.to_f64().unwrap() / 2., sum.im.to_f64().unwrap() / 2.);
            }

            stable = result.is_finite()
                && stable_digits >= stability_check.relative_precision
                && abs_diff <= self.settings.general.absolute_precision
                && !escalate;

            stability_level = level;

            if stable {
                break;
            }

            self.integrand_statistics.unstable_point_count[level] += 1;
        }

        std::mem::swap(
            &mut self.settings.general.stability_checks,
            &mut stability_checks,
        );
        std::mem::swap(&mut self.cache, &mut cache);

        self.integrand_statistics.total_samples += 1;

        // check if we have deemed the point unstable
        if !stable {
            if !result.is_finite() {
                self.status_update_sender
                    .send(StatusUpdate::Message(format!("NaN point {:?}", x)))
                    .unwrap();

                self.integrand_statistics.nan_point_count += 1;
                event_manager.clear(true); // throw away all events and treat them as rejected
                result = Complex::default();
            } else {
                if stable_digits < self.settings.general.minimal_precision_for_returning_result {
                    self.status_update_sender
                        .send(StatusUpdate::Message(format!(
                            "Unstable point {:?}: {:.5e}",
                            x,
                            result * weight
                        )))
                        .unwrap();

                    event_manager.clear(true); // throw away all events and treat them as rejected
                    result = Complex::default();
                }
            }
        } else {
            self.integrand_statistics.regular_point_count += 1;
        }

        if self.integrand_statistics.running_max_re.0 < result.re.abs() * weight {
            self.integrand_statistics.running_max_re =
                (result.re.abs() * weight, weight, stability_level);
            // TODO: replace by full sample point logging
            let len = x.to_flat().len(); // could be shorter than self.n_loops
            self.integrand_statistics.running_max_coordinate_re[..len].copy_from_slice(x.to_flat());

            if self.settings.integrator.integrated_phase != IntegratedPhase::Imag {
                self.integrand_statistics.running_max_stability = stable_digits;
            }
        }
        if self.integrand_statistics.running_max_im.0 < result.im.abs() * weight {
            self.integrand_statistics.running_max_im =
                (result.im.abs() * weight, weight, stability_level);
            let len = x.to_flat().len();
            self.integrand_statistics.running_max_coordinate_im[..len].copy_from_slice(x.to_flat());

            if self.settings.integrator.integrated_phase != IntegratedPhase::Real {
                self.integrand_statistics.running_max_stability = stable_digits;
            }
        }

        event_manager.process_events(result, weight);
        std::mem::swap(&mut self.event_manager, &mut event_manager);

        self.integrand_statistics.total_sample_time +=
            Instant::now().duration_since(start_time).as_secs_f64() * 1e6;

        result
    }
}
