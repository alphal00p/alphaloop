use arrayvec::ArrayVec;
use f128::f128;
use float;
use num::Complex;
use num_traits::ops::inv::Inv;
use num_traits::{Float, FloatConst, FromPrimitive, NumCast, One, ToPrimitive, Zero};
use observables;
use observables::{Event, EventFilter, Observable};
use rand::Rng;
use std::fs::File;
use std::io::{BufWriter, Write};
use topologies::{CachePrecisionSelector, LTDCacheAllPrecisions, Topology};
use vector::LorentzVector;
use {FloatLike, IntegratedPhase, ObservableMode, Settings, MAX_LOOP};

pub trait IntegrandImplementation: Clone {
    type Cache: Default;

    fn rotate(&self, angle: float, axis: (float, float, float)) -> Self;
    fn evaluate_float<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut Self::Cache,
        events: Option<&mut Vec<Event>>,
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
        events: Option<&mut Vec<Event>>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<f128>>; MAX_LOOP]>,
        f128,
        Complex<f128>,
        Complex<f128>,
    );

    fn create_cache(&self) -> Self::Cache;
}

/// A structure that integrates and keeps statistics
pub struct Integrand<I: IntegrandImplementation> {
    pub n_loops: usize,
    pub settings: Settings,
    pub topologies: Vec<I>,
    pub cache: I::Cache,
    pub running_max: Complex<float>,
    pub running_min_lambda: float,
    pub total_samples: usize,
    pub nan_point_count: usize,
    pub unstable_point_count: usize,
    pub unstable_f128_point_count: usize,
    pub regular_point_count: usize,
    pub log: BufWriter<File>,
    pub quadruple_upgrade_log: BufWriter<File>,
    pub cur_iter: usize,
    pub id: usize,
    pub event_filter: Box<dyn EventFilter + Send>,
    pub observables: Vec<Box<dyn Observable + Send>>,
    pub events: Vec<Event>, // event buffer,
    pub track_events: bool,
}

impl IntegrandImplementation for Topology {
    type Cache = LTDCacheAllPrecisions;

    /// Create a rotated version of this topology. The axis needs to be normalized.
    fn rotate(&self, angle: float, axis: (float, float, float)) -> Topology {
        let cos_t = angle.cos();
        let sin_t = angle.sin();
        let cos_t_bar = float::one() - angle.cos();

        let rot_matrix: [[float; 3]; 3] = [
            [
                cos_t + axis.0 * axis.0 * cos_t_bar,
                axis.0 * axis.1 * cos_t_bar - axis.2 * sin_t,
                axis.0 * axis.2 * cos_t_bar + axis.1 * sin_t,
            ],
            [
                axis.0 * axis.1 * cos_t_bar + axis.2 * sin_t,
                cos_t + axis.1 * axis.1 * cos_t_bar,
                axis.1 * axis.2 * cos_t_bar - axis.0 * sin_t,
            ],
            [
                axis.0 * axis.2 * cos_t_bar - axis.1 * sin_t,
                axis.1 * axis.2 * cos_t_bar + axis.0 * sin_t,
                cos_t + axis.2 * axis.2 * cos_t_bar,
            ],
        ];

        let mut rotated_topology = self.clone();
        rotated_topology.name = rotated_topology.name + "_rot";
        rotated_topology.rotation_matrix = rot_matrix.clone();

        for e in &mut rotated_topology.external_kinematics {
            let old_x = float::from_f64(e.x).unwrap();
            let old_y = float::from_f64(e.y).unwrap();
            let old_z = float::from_f64(e.z).unwrap();
            e.x = (rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z)
                .to_f64()
                .unwrap();
            e.y = (rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z)
                .to_f64()
                .unwrap();
            e.z = (rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z)
                .to_f64()
                .unwrap();
        }

        for ll in &mut rotated_topology.loop_lines {
            for p in &mut ll.propagators {
                let old_x = float::from_f64(p.q.x).unwrap();
                let old_y = float::from_f64(p.q.y).unwrap();
                let old_z = float::from_f64(p.q.z).unwrap();
                p.q.x = (rot_matrix[0][0] * old_x
                    + rot_matrix[0][1] * old_y
                    + rot_matrix[0][2] * old_z)
                    .to_f64()
                    .unwrap();
                p.q.y = (rot_matrix[1][0] * old_x
                    + rot_matrix[1][1] * old_y
                    + rot_matrix[1][2] * old_z)
                    .to_f64()
                    .unwrap();
                p.q.z = (rot_matrix[2][0] * old_x
                    + rot_matrix[2][1] * old_y
                    + rot_matrix[2][2] * old_z)
                    .to_f64()
                    .unwrap();
            }
        }

        for surf in &mut rotated_topology.surfaces {
            let old_x = surf.shift.x;
            let old_y = surf.shift.y;
            let old_z = surf.shift.z;
            surf.shift.x =
                rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z;
            surf.shift.y =
                rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z;
            surf.shift.z =
                rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z;
        }

        // now rotate the fixed deformation vectors
        for d_lim in &mut rotated_topology.fixed_deformation {
            for d in &mut d_lim.deformation_per_overlap {
                for source in &mut d.deformation_sources {
                    let old_x = source.x;
                    let old_y = source.y;
                    let old_z = source.z;
                    source.x = rot_matrix[0][0] * old_x
                        + rot_matrix[0][1] * old_y
                        + rot_matrix[0][2] * old_z;
                    source.y = rot_matrix[1][0] * old_x
                        + rot_matrix[1][1] * old_y
                        + rot_matrix[1][2] * old_z;
                    source.z = rot_matrix[2][0] * old_x
                        + rot_matrix[2][1] * old_y
                        + rot_matrix[2][2] * old_z;
                }
            }
        }
        // now rotate the numerators
        rotated_topology.numerator = rotated_topology.numerator.rotate(rot_matrix);
        if rotated_topology.settings.general.use_amplitude {
            for diag in rotated_topology.amplitude.diagrams.iter_mut() {
                diag.numerator = diag.numerator.rotate(rot_matrix);
            }
        }

        rotated_topology
    }

    #[inline]
    fn evaluate_float<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut LTDCacheAllPrecisions,
        _events: Option<&mut Vec<Event>>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]>,
        float,
        Complex<float>,
        Complex<float>,
    ) {
        self.evaluate(x, cache.get())
    }

    #[inline]
    fn evaluate_f128<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut LTDCacheAllPrecisions,
        _events: Option<&mut Vec<Event>>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<f128>>; MAX_LOOP]>,
        f128,
        Complex<f128>,
        Complex<f128>,
    ) {
        self.evaluate(x, cache.get())
    }

    #[inline]
    fn create_cache(&self) -> LTDCacheAllPrecisions {
        LTDCacheAllPrecisions::new(self)
    }
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
    ) -> ($ty, $ty, Complex<$ty>, Complex<$ty>) {
        if !self.settings.general.numerical_instability_check
            || self.topologies.len() == 1
            || num_samples == 1
        {
            return (
                <$ty>::from_f64(self.settings.general.relative_precision).unwrap(),
                <$ty>::from_f64(self.settings.general.absolute_precision).unwrap(),
                Complex::default(),
                Complex::default(),
            );
        }

        // collect smallest result and biggest result
        let mut min = result;
        let mut max = result;

        // TODO: we are also findings the centers for the rotated topologies
        // inherit them from the first topology
        for rot_topo in &mut self.topologies[1..num_samples] {
            if self.settings.general.debug > 2 {
                println!("Evaluating integrand with rotated topologies");
            }

            let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, result_rot) =
                rot_topo.$eval_fn(x, cache, None);

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
        }

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

        (d, Float::abs(num - num_rot), min, max)
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
        id: usize,
    ) -> Integrand<I> {
        // create extra topologies with rotated kinematics to check the uncertainty
        let mut rng = rand::thread_rng();
        let mut topologies = vec![topology.clone()];
        for _ in 0..settings
            .general
            .num_f64_samples
            .max(settings.general.num_f128_samples)
        {
            let angle =
                float::from_f64(rng.gen::<f64>() * 2.).unwrap() * <float as FloatConst>::PI();
            let mut rv = (
                float::from_f64(rng.gen()).unwrap(),
                float::from_f64(rng.gen()).unwrap(),
                float::from_f64(rng.gen()).unwrap(),
            ); // rotation axis
            let inv_norm = (rv.0 * rv.0 + rv.1 * rv.1 + rv.2 * rv.2).sqrt().inv();
            rv = (rv.0 * inv_norm, rv.1 * inv_norm, rv.2 * inv_norm);

            topologies.push(topology.rotate(angle, rv));
        }

        let mut observables: Vec<Box<dyn Observable + Send>> = vec![];
        for o in &settings.observables.active_observables {
            match o {
                ObservableMode::Jet1PT => {
                    observables.push(Box::new(observables::Jet1PTObservable::new(
                        settings.observables.Jet1PT.x_min,
                        settings.observables.Jet1PT.x_max,
                        settings.observables.Jet1PT.n_bins,
                        settings.observables.Jet1PT.dR,
                        settings.observables.Jet1PT.write_to_file,
                        settings.observables.Jet1PT.filename.clone(),
                    )));
                }
                ObservableMode::CrossSection => {
                    observables.push(Box::new(observables::CrossSectionObservable::default()));
                }
            }
        }

        Integrand {
            n_loops,
            topologies,
            cache: topology.create_cache(),
            running_max: Complex::default(),
            total_samples: 0,
            regular_point_count: 0,
            unstable_point_count: 0,
            unstable_f128_point_count: 0,
            nan_point_count: 0,
            running_min_lambda: float::infinity(),
            log: BufWriter::new(
                File::create(format!("{}{}.log", settings.general.log_file_prefix, id))
                    .expect("Could not create log file"),
            ),
            quadruple_upgrade_log: BufWriter::new(
                File::create(format!(
                    "{}_f128_{}.log",
                    settings.general.log_file_prefix, id
                ))
                .expect("Could not create log file"),
            ),
            settings,
            id,
            cur_iter: 1, // always start at iter 1
            event_filter: Box::new(observables::NoEventFilter::default()),
            observables,
            events: Vec::with_capacity(100),
            track_events,
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
        if new_max || unstable || !result.is_finite() || self.settings.general.debug > 0 {
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
        if self.settings.general.log_stats_to_screen
            && (self.settings.general.screen_log_core == None
                || self.settings.general.screen_log_core == Some(self.id))
        {
            eprintln!(
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | unstable f128 points={} ({:.2}%)\n  | nan points={} ({:.2}%)\n  | running min lambda={}",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_f128_point_count, self.unstable_f128_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100., self.running_min_lambda);
        }

        writeln!(self.log,
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | unstable f128 points={} ({:.2}%)\n  | nan points={} ({:.2}%)\n  | running min lambda={}",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_f128_point_count, self.unstable_f128_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100.,
            self.running_min_lambda).unwrap();
    }

    check_stability_precision!(check_stability_float, evaluate_float, float);
    check_stability_precision!(check_stability_quad, evaluate_f128, f128);

    /// Evalute a point generated from the Monte Carlo generator with `weight` and current iteration number `iter`.
    pub fn evaluate(&mut self, x: &[f64], weight: f64, iter: usize) -> Complex<float> {
        if self.cur_iter != iter {
            if self.track_events {
                for o in &mut self.observables {
                    o.update_result();
                }
            }

            self.cur_iter += 1;
        }

        if self.settings.general.integration_statistics
            && self.total_samples > 0
            && self.total_samples % self.settings.general.statistics_interval == 0
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
        if self.unstable_point_count as f64 / self.total_samples as f64 * 100.
            > self.settings.general.unstable_point_warning_percentage
            && self.total_samples % self.settings.general.statistics_interval == (self.id + 1) * 20
        {
            eprintln!(
                "WARNING on core {}: {:.2}% of points are unstable, {:.2}% are not saved by f128",
                self.id,
                self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
                self.unstable_f128_point_count as f64 / self.total_samples as f64 * 100.
            );
        }

        let mut events = std::mem::replace(&mut self.events, vec![]);

        let mut cache = std::mem::replace(&mut self.cache, I::Cache::default());
        let (x, k_def, jac_para, jac_def, mut result) =
            self.topologies[0].evaluate_float(x, &mut cache, Some(&mut events));
        let (d, diff, min_rot, max_rot) = self.check_stability_float(
            x,
            result,
            self.settings.general.num_f64_samples,
            &mut cache,
        );
        std::mem::swap(&mut self.cache, &mut cache);

        // FIXME
        /*if self.cache.overall_lambda < self.running_min_lambda {
            self.running_min_lambda = self.cache.overall_lambda;
        }*/

        self.total_samples += 1;

        if self.settings.general.force_f128
            || !result.is_finite()
            || !min_rot.is_finite()
            || !max_rot.is_finite()
            || d < NumCast::from(self.settings.general.relative_precision).unwrap()
            || diff > NumCast::from(self.settings.general.absolute_precision).unwrap()
        {
            // clear events when there is instability
            events.clear();

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
                self.topologies[0].evaluate_f128(x, &mut cache, Some(&mut events));
            // NOTE: for this check we use a f64 rotation matrix at the moment!
            let (d_f128, diff_f128, min_rot_f128, max_rot_f128) = self.check_stability_quad(
                x,
                result_f128,
                self.settings.general.num_f128_samples,
                &mut cache,
            );
            std::mem::swap(&mut self.cache, &mut cache);
            /*
            // check if the f128 computation is consistent with the f64 one
            // by checking if the f128 result is more than 2 stable digits removed from the f64 one
            let (num_f64, num_f128) =
                if self.settings.integrator.integrated_phase == IntegratedPhase::Imag {
                    (
                        <float as NumCast>::from(result.im).unwrap(),
                        <float as NumCast>::from(result_f128.im).unwrap(),
                    )
                } else {
                    (
                        <float as NumCast>::from(result.re).unwrap(),
                        <float as NumCast>::from(result_f128.re).unwrap(),
                    )
                };

            let d_comp = if num_f64.is_zero() && num_f128.is_zero() {
                float::from_f64(100.).unwrap()
            } else {
                if num_f64 == -num_f128 {
                    float::from_f64(-100.).unwrap()
                } else {
                    -((num_f64 - num_f128) / (num_f64 + num_f128)).abs().log10()
                }
            };

            if d_comp
                + float::from_f64(self.settings.general.num_digits_different_for_inconsistency)
                    .unwrap()
                < d
                && d > float::from_f64(3.).unwrap()
                && num_f64 > float::from_f64(1e-20).unwrap()
            {
                eprintln!(
                    "{} between f64 and f128!\n  | f64  ={:e}, jac={:e}\n  | f64' ={:e}, jac={:e}\n  | f128 ={:e}, jac={:e}\n  | f128'={:e}, jac={:e}\n  | x={:?}",
                    "Inconsistency".red(),
                    result, jac_def * jac_para, result_rot, 0., result_f128, jac_def_f128 * jac_para_f128, result_rot_f128, 0., x,
                );

                // now try more points for f64
                for (i, topo) in self.topologies[2..].iter().enumerate() {
                    let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, result_rot) =
                        topo.evaluate(x, &mut self.cache_float, &self.numerator);
                    eprintln!("  | rot{} ={:e}", i + 2, result_rot);
                }
            }*/

            self.unstable_point_count += 1;

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
                    self.nan_point_count += 1;
                    return Complex::default();
                } else {
                    self.unstable_f128_point_count += 1;
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
                        return Complex::default();
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
            self.regular_point_count += 1;
        }

        let mut new_max = false;
        match self.settings.integrator.integrated_phase {
            IntegratedPhase::Real => {
                if result.re.abs() > self.running_max.re.abs() {
                    self.running_max = result;
                    new_max = true;
                }
            }
            IntegratedPhase::Imag => {
                if result.im.abs() > self.running_max.im.abs() {
                    self.running_max = result;
                    new_max = true;
                }
            }
            IntegratedPhase::Both => {
                if result.norm_sqr() > self.running_max.norm_sqr() {
                    self.running_max = result;
                    new_max = true;
                }
            }
        }

        if self.settings.general.integration_statistics {
            let loops = self.n_loops;
            self.print_info(
                loops, new_max, false, x, &k_def, jac_para, jac_def, min_rot, max_rot, d,
            );
        }
        if self.track_events {
            // filter the events, potentially changing the integrand result
            result = self.event_filter.process_event_group(&mut events, weight);

            // give the events to an observable function
            for o in &mut self.observables {
                o.process_event_group(&mut events, weight);
            }
            events.clear();
        }

        std::mem::swap(&mut self.events, &mut events);

        result
    }
}
