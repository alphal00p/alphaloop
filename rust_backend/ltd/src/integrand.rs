use arrayvec::ArrayVec;
use color_eyre::{Help, Report};
use dashboard::{StatusUpdate, StatusUpdateSender};
use eyre::WrapErr;
use f128::f128;
use float;
use num::Complex;
use num_traits::{Float, FromPrimitive, NumCast, ToPrimitive, Zero};
use observables::EventManager;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;
use vector::LorentzVector;
use {FloatLike, IntegratedPhase, Settings, MAX_LOOP};

pub trait IntegrandImplementation: Clone {
    type Cache: Default;

    fn create_stability_check(&self, num_checks: usize) -> Vec<Self>;

    fn evaluate_float<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut Self::Cache,
        events: Option<&mut EventManager>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]>,
        float,
        Complex<float>,
        Complex<float>,
    );

    fn evaluate_f128<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut Self::Cache,
        events: Option<&mut EventManager>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<f128>>; MAX_LOOP]>,
        f128,
        Complex<f128>,
        Complex<f128>,
    );

    fn create_cache(&self) -> Self::Cache;
}

#[derive(Debug, Copy, Clone)]
pub struct IntegrandStatistics {
    pub phase: IntegratedPhase,
    pub running_max: Complex<float>,
    pub total_samples: usize,
    pub nan_point_count: usize,
    pub unstable_point_count: usize,
    pub unstable_f128_point_count: usize,
    pub regular_point_count: usize,
    pub total_sample_time: f64,
    pub n_loops: usize,
    pub running_max_coordinate_re: [f64; 3 * MAX_LOOP],
    pub running_max_coordinate_im: [f64; 3 * MAX_LOOP],
    pub running_max_stability: (f64, f64),
}

impl IntegrandStatistics {
    pub fn new(n_loops: usize, phase: IntegratedPhase) -> IntegrandStatistics {
        IntegrandStatistics {
            n_loops,
            phase,
            running_max: Complex::default(),
            total_samples: 0,
            regular_point_count: 0,
            unstable_point_count: 0,
            unstable_f128_point_count: 0,
            nan_point_count: 0,
            total_sample_time: 0.,
            running_max_coordinate_re: [0.; 3 * MAX_LOOP],
            running_max_coordinate_im: [0.; 3 * MAX_LOOP],
            running_max_stability: (0., 0.),
        }
    }

    pub fn merge(&mut self, other: &mut IntegrandStatistics) {
        self.total_samples += other.total_samples;
        self.regular_point_count += other.regular_point_count;
        self.unstable_point_count += other.unstable_point_count;
        self.unstable_f128_point_count += other.unstable_f128_point_count;
        self.nan_point_count += other.nan_point_count;
        self.total_sample_time += other.total_sample_time;

        if self.running_max.re < other.running_max.re {
            self.running_max.re = other.running_max.re;
            if self.phase != IntegratedPhase::Imag {
                self.running_max_stability = other.running_max_stability;
            }
        }
        if self.running_max.im < other.running_max.im {
            self.running_max.im = other.running_max.im;
            if self.phase != IntegratedPhase::Real {
                self.running_max_stability = other.running_max_stability;
            }
        }

        other.total_samples = 0;
        other.regular_point_count = 0;
        other.unstable_point_count = 0;
        other.unstable_f128_point_count = 0;
        other.nan_point_count = 0;
        other.total_sample_time = 0.;
    }
}

/// A structure that integrates and keeps statistics
pub struct Integrand<I: IntegrandImplementation> {
    pub n_loops: usize,
    pub settings: Settings,
    pub topologies: Vec<I>,
    pub cache: I::Cache,
    pub log: BufWriter<File>,
    pub quadruple_upgrade_log: BufWriter<File>,
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
        x: &[f64],
        result: Complex<$ty>,
        num_samples: usize,
        cache: &mut I::Cache,
        event_manager: &mut EventManager,
    ) -> ($ty, $ty, Complex<$ty>, Complex<$ty>) {
        let mut ret = (
            <$ty>::from_f64(self.settings.general.relative_precision).unwrap(),
            <$ty>::from_f64(self.settings.general.absolute_precision).unwrap(),
            Complex::default(),
            Complex::default(),
        );

        if !self.settings.general.numerical_instability_check
            || self.topologies.len() == 1
            || num_samples == 1
        {
            return ret;
        }

        // collect smallest result and biggest result
        let mut min = result;
        let mut max = result;

        // do not track events of rotated topologies
        let track_events = event_manager.track_events;
        event_manager.track_events = false;

        // TODO: we are also findings the centers for the rotated topologies
        // inherit them from the first topology
        for (rot_index, rot_topo) in self.topologies[1..num_samples].iter_mut().enumerate() {
            if self.settings.general.debug > 2 {
                println!("Evaluating integrand with rotated topologies");
            }

            let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, result_rot) =
                rot_topo.$eval_fn(x, cache, Some(event_manager));

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

            ret = (d, Float::abs(num - num_rot), min, max);

            // if we see the point is unstable, stop further samples
            if !result_rot.is_finite()
                || !min.is_finite()
                || !max.is_finite()
                || d < NumCast::from(self.settings.general.relative_precision).unwrap()
                || diff > NumCast::from(self.settings.general.absolute_precision).unwrap() {
                break;
            }

            // if the first point is very stable, we stop further samples
            if rot_index == 0 && result_rot.is_finite() && d > NumCast::from(self.settings.general.minimal_precision_to_skip_further_checks).unwrap() {
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
        mut settings: Settings,
        track_events: bool,
        status_update_sender: StatusUpdateSender,
        id: usize,
    ) -> Integrand<I> {
        // create extra topologies with rotated kinematics to check the uncertainty
        let mut rng = rand::thread_rng();
        let mut topologies = vec![topology.clone()];

        if settings.general.force_f128 {
            settings.general.num_f64_samples = 1;
        }

        topologies.extend(
            topology.create_stability_check(
                settings
                    .general
                    .num_f64_samples
                    .max(settings.general.num_f128_samples),
            ),
        );

        let log_filename = format!("{}{}.log", settings.general.log_file_prefix, id);
        let quad_log_filename = format!("{}_f128_{}.log", settings.general.log_file_prefix, id);

        Integrand {
            n_loops,
            topologies,
            cache: topology.create_cache(),
            integrand_statistics: IntegrandStatistics::new(
                n_loops,
                settings.integrator.integrated_phase,
            ),
            log: BufWriter::new(
                File::create(&log_filename)
                    .wrap_err_with(|| format!("Could not create log file {}", log_filename))
                    .suggestion("Check if this topology is in the specified topology file.")
                    .unwrap(),
            ),
            quadruple_upgrade_log: BufWriter::new(
                File::create(&quad_log_filename)
                    .wrap_err_with(|| format!("Could not create log file {}", quad_log_filename))
                    .suggestion("Check if this topology is in the specified topology file.")
                    .unwrap(),
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

    fn print_info<T: FloatLike>(
        &mut self,
        n_loops: usize,
        new_max: bool,
        unstable: bool,
        x: &[f64],
        k_def: &ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        jac_para: T,
        jac_def: Complex<T>,
        result: Complex<T>,
        rot_result: Complex<T>,
        stable_digits: T,
    ) {
        if new_max
            || unstable
            || !result.re.is_finite()
            || !result.im.is_finite()
            || self.settings.general.debug > 0
        {
            let sample_or_max = if new_max { "MAX" } else { "Sample" };
            let log_to_screen = self.settings.general.log_points_to_screen
                && (self.settings.general.screen_log_core == None
                    || self.settings.general.screen_log_core == Some(self.id));

            match n_loops {
                1 => {
                    if log_to_screen {
                        eprintln!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], jac_para, jac_def
                    );
                    }
                    writeln!(self.log, "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], jac_para, jac_def).unwrap();
                }
                2 => {
                    if log_to_screen {
                        eprintln!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], jac_para, jac_def
                    );
                    }
                    writeln!(self.log, "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], jac_para, jac_def).unwrap();
                }
                3 => {
                    if log_to_screen {
                        eprintln!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    );
                    }
                    writeln!(self.log,
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    ).unwrap();
                }
                4 => {
                    if log_to_screen {
                        eprintln!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | n={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], k_def[2], k_def[3], jac_para, jac_def
                    );
                    }
                    writeln!(self.log,
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | n={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], k_def[2], k_def[3], jac_para, jac_def
                    ).unwrap();
                }
                _ => {}
            }
        }
    }

    fn print_statistics(&mut self) {
        let s = &self.integrand_statistics;

        if self.settings.general.log_stats_to_screen
            && (self.settings.general.screen_log_core == None
                || self.settings.general.screen_log_core == Some(self.id))
        {
            eprintln!(
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | unstable f128 points={} ({:.2}%)\n  | nan points={} ({:.2}%)\n",
            s.running_max, s.total_samples, s.regular_point_count, s.regular_point_count as f64 / s.total_samples as f64 * 100.,
            s.unstable_point_count, s.unstable_point_count as f64 / s.total_samples as f64 * 100.,
            s.unstable_f128_point_count, s.unstable_f128_point_count as f64 / s.total_samples as f64 * 100.,
            s.nan_point_count, s.nan_point_count as f64 / s.total_samples as f64 * 100.);
        }

        writeln!(self.log,
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | unstable f128 points={} ({:.2}%)\n  | nan points={} ({:.2}%)\n",
            s.running_max, s.total_samples, s.regular_point_count, s.regular_point_count as f64 / s.total_samples as f64 * 100.,
            s.unstable_point_count, s.unstable_point_count as f64 / s.total_samples as f64 * 100.,
            s.unstable_f128_point_count, s.unstable_f128_point_count as f64 / s.total_samples as f64 * 100.,
            s.nan_point_count, s.nan_point_count as f64 / s.total_samples as f64 * 100.).unwrap();
    }

    check_stability_precision!(check_stability_float, evaluate_float, float);
    check_stability_precision!(check_stability_quad, evaluate_f128, f128);

    pub fn merge_statistics(&mut self, other: &mut Integrand<I>) {
        self.event_manager.merge_samples(&mut other.event_manager);

        self.integrand_statistics
            .merge(&mut other.integrand_statistics);
    }

    pub fn broadcast_statistics(&mut self) {
        if self.id == 0 {
            self.status_update_sender
                .send(StatusUpdate::Statistics(self.integrand_statistics))
                .unwrap();

            self.event_manager.update_live_result();
        }
    }

    /// Evalute a point generated from the Monte Carlo generator with `weight` and current iteration number `iter`.
    pub fn evaluate(&mut self, x: &[f64], mut weight: f64, iter: usize) -> Complex<float> {
        let start_time = Instant::now(); // time the evaluation

        // NOTE: only correct for Vegas
        weight *= (self.settings.integrator.n_start
            + (self.cur_iter - 1) * self.settings.integrator.n_increase) as f64;

        if self.cur_iter != iter {
            // the first integrand accumulates all the results from the others
            if self.id == 0 {
                self.event_manager.update_result();

                self.status_update_sender
                    .send(StatusUpdate::Statistics(self.integrand_statistics))
                    .unwrap();
            }

            self.cur_iter += 1;
        }

        if self.settings.general.integration_statistics
            && self.integrand_statistics.total_samples > 0
            && self.integrand_statistics.total_samples % self.settings.general.statistics_interval
                == 0
        {
            self.print_statistics();
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

        // print warnings for unstable points and try to make sure the screen output does not
        // get corrupted
        if !self.settings.integrator.internal_parallelization
            && self.integrand_statistics.unstable_point_count as f64
                / self.integrand_statistics.total_samples as f64
                * 100.
                > self.settings.general.unstable_point_warning_percentage
            && self.integrand_statistics.total_samples % self.settings.general.statistics_interval
                == (self.id + 1) * 20
        {
            self.status_update_sender
                .send(StatusUpdate::Message(format!(
                "WARNING on core {}: {:.2}% of points are unstable, {:.2}% are not saved by f128",
                self.id,
                self.integrand_statistics.unstable_point_count as f64 / self.integrand_statistics.total_samples as f64 * 100.,
                self.integrand_statistics.unstable_f128_point_count as f64 / self.integrand_statistics.total_samples as f64 * 100.
            )))
                .unwrap();
        }

        let mut event_manager = std::mem::replace(&mut self.event_manager, EventManager::default());

        let mut cache = std::mem::replace(&mut self.cache, I::Cache::default());
        let (x, k_def, jac_para, jac_def, mut result) =
            self.topologies[0].evaluate_float(x, &mut cache, Some(&mut event_manager));
        let (d, diff, min_rot, max_rot) = self.check_stability_float(
            x,
            result,
            self.settings.general.num_f64_samples,
            &mut cache,
            &mut event_manager,
        );
        std::mem::swap(&mut self.cache, &mut cache);
        self.integrand_statistics.total_samples += 1;

        let do_f128 = match self.settings.integrator.integrated_phase {
            IntegratedPhase::Real => {
                result.re.abs() * weight
                    > self.integrand_statistics.running_max.re
                        * self.settings.general.force_f128_for_large_weight_threshold
            }
            IntegratedPhase::Both => {
                result.re.abs() * weight
                    > self.integrand_statistics.running_max.re
                        * self.settings.general.force_f128_for_large_weight_threshold
                    || result.im.abs() * weight
                        > self.integrand_statistics.running_max.im
                            * self.settings.general.force_f128_for_large_weight_threshold
            }
            IntegratedPhase::Imag => {
                result.im.abs() * weight
                    > self.integrand_statistics.running_max.im
                        * self.settings.general.force_f128_for_large_weight_threshold
            }
        };

        let mut stability = (d, 0.);
        if self.settings.general.force_f128
            || do_f128
            || !result.is_finite()
            || !min_rot.is_finite()
            || !max_rot.is_finite()
            || d < NumCast::from(self.settings.general.relative_precision).unwrap()
            || diff > NumCast::from(self.settings.general.absolute_precision).unwrap()
        {
            // clear events when there is instability
            event_manager.clear(false);

            if self.settings.general.integration_statistics {
                let loops = self.n_loops;
                self.print_info(
                    loops, false, true, x, &k_def, jac_para, jac_def, min_rot, max_rot, d,
                );
            }

            if self.settings.general.log_quad_upgrade {
                match self.n_loops {
                    1 => writeln!(
                        self.quadruple_upgrade_log,
                        "{:?} {} {} {}",
                        x, k_def[0].x.re, k_def[0].y.re, k_def[0].z.re
                    )
                    .unwrap(),
                    2 => writeln!(
                        self.quadruple_upgrade_log,
                        "{} {} {} {} {} {}",
                        k_def[0].x.re,
                        k_def[0].y.re,
                        k_def[0].z.re,
                        k_def[1].x.re,
                        k_def[1].y.re,
                        k_def[1].z.re
                    )
                    .unwrap(),
                    3 => writeln!(
                        self.quadruple_upgrade_log,
                        "{} {} {} {} {} {} {} {} {}",
                        k_def[0].x.re,
                        k_def[0].y.re,
                        k_def[0].z.re,
                        k_def[1].x.re,
                        k_def[1].y.re,
                        k_def[1].z.re,
                        k_def[2].x.re,
                        k_def[2].y.re,
                        k_def[2].z.re
                    )
                    .unwrap(),
                    4 => writeln!(
                        self.quadruple_upgrade_log,
                        "{} {} {} {} {} {} {} {} {} {} {} {}",
                        k_def[0].x.re,
                        k_def[0].y.re,
                        k_def[0].z.re,
                        k_def[1].x.re,
                        k_def[1].y.re,
                        k_def[1].z.re,
                        k_def[2].x.re,
                        k_def[2].y.re,
                        k_def[2].z.re,
                        k_def[3].x.re,
                        k_def[3].y.re,
                        k_def[3].z.re
                    )
                    .unwrap(),
                    _ => {}
                }
                self.quadruple_upgrade_log.flush().unwrap();
            }

            // compute the point again with f128 to see if it is stable then
            std::mem::swap(&mut self.cache, &mut cache);
            let (_, k_def_f128, jac_para_f128, jac_def_f128, result_f128) =
                self.topologies[0].evaluate_f128(x, &mut cache, Some(&mut event_manager));
            // NOTE: for this check we use a f64 rotation matrix at the moment!
            let (d_f128, diff_f128, min_rot_f128, max_rot_f128) = self.check_stability_quad(
                x,
                result_f128,
                self.settings.general.num_f128_samples,
                &mut cache,
                &mut event_manager,
            );
            std::mem::swap(&mut self.cache, &mut cache);
            self.integrand_statistics.unstable_point_count += 1;
            stability = (stability.0, d_f128.to_f64().unwrap());

            if !result_f128.is_finite()
                || !min_rot_f128.is_finite()
                || !max_rot_f128.is_finite()
                || d_f128 < NumCast::from(self.settings.general.relative_precision).unwrap()
                || diff_f128 > NumCast::from(self.settings.general.absolute_precision).unwrap()
            {
                if self.settings.general.integration_statistics {
                    let loops = self.n_loops;
                    self.print_info(
                        loops,
                        false,
                        true,
                        x,
                        &k_def_f128,
                        jac_para_f128,
                        jac_def_f128,
                        min_rot_f128,
                        max_rot_f128,
                        d_f128,
                    );
                }

                if !result_f128.is_finite() {
                    self.status_update_sender
                        .send(StatusUpdate::Message(format!("NaN point {:?}", x)))
                        .unwrap();

                    self.integrand_statistics.nan_point_count += 1;
                    event_manager.clear(true); // throw away all events and treat them as rejected
                    result = Complex::default();
                } else {
                    self.integrand_statistics.unstable_f128_point_count += 1;
                    //println!("f128 fail : x={:?}, min result={}, max result={}", x, min_rot_f128, max_rot_f128);

                    if d_f128
                        > NumCast::from(
                            self.settings.general.minimal_precision_for_returning_result,
                        )
                        .unwrap()
                    {
                        result = Complex::new(
                            <float as NumCast>::from(result_f128.re).unwrap(),
                            <float as NumCast>::from(result_f128.im).unwrap(),
                        );
                    } else {
                        self.status_update_sender
                            .send(StatusUpdate::Message(format!(
                                "Unstable quad point {:?}: {:.5e}",
                                x,
                                result_f128 * f128::from_f64(weight).unwrap()
                            )))
                            .unwrap();

                        event_manager.clear(true); // throw away all events and treat them as rejected
                        result = Complex::default();
                    }
                }
            } else {
                // we have saved the integration!
                // TODO: also modify the other parameters for the print_info below?
                result = Complex::new(
                    <float as NumCast>::from(result_f128.re).unwrap(),
                    <float as NumCast>::from(result_f128.im).unwrap(),
                );
            }
        } else {
            self.integrand_statistics.regular_point_count += 1;
        }

        let mut new_max = false;
        if self.integrand_statistics.running_max.re < result.re.abs() * weight {
            new_max = true;
            self.integrand_statistics.running_max.re = result.re.abs() * weight;
            self.integrand_statistics.running_max_coordinate_re[..3 * self.n_loops]
                .copy_from_slice(x);

            if self.settings.integrator.integrated_phase != IntegratedPhase::Imag {
                self.integrand_statistics.running_max_stability = stability;
            }
        }
        if self.integrand_statistics.running_max.im < result.im.abs() * weight {
            new_max = true;
            self.integrand_statistics.running_max.im = result.im.abs() * weight;
            self.integrand_statistics.running_max_coordinate_im[..3 * self.n_loops]
                .copy_from_slice(x);

            if self.settings.integrator.integrated_phase != IntegratedPhase::Real {
                self.integrand_statistics.running_max_stability = stability;
            }
        }

        if self.settings.general.integration_statistics {
            let loops = self.n_loops;
            self.print_info(
                loops, new_max, false, x, &k_def, jac_para, jac_def, min_rot, max_rot, d,
            );
        }

        event_manager.process_events(result, weight);
        std::mem::swap(&mut self.event_manager, &mut event_manager);

        self.integrand_statistics.total_sample_time +=
            Instant::now().duration_since(start_time).as_secs_f64() * 1e6;

        result
    }
}
