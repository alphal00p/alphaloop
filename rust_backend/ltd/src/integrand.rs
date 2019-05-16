use arrayvec::ArrayVec;
use colored::Colorize;
use f128::f128;
use float;
use num::Complex;
use num_traits::{Float, FloatConst, FromPrimitive, Inv, NumCast, One, ToPrimitive, Zero};
use rand::Rng;
use std::fs::File;
use std::io::{BufWriter, Write};
use topologies::{LTDCache, Topology};
use vector::LorentzVector;

use {FloatLike, IntegratedPhase, PythonNumerator, Settings, MAX_LOOP};

/// A structure that integrates and keeps statistics
pub struct Integrand {
    pub settings: Settings,
    pub topologies: Vec<Topology>,
    pub cache_float: LTDCache<float>,
    pub cache_f128: LTDCache<f128>,
    pub running_max: Complex<float>,
    pub total_samples: usize,
    pub nan_point_count: usize,
    pub unstable_point_count: usize,
    pub unstable_f128_point_count: usize,
    pub regular_point_count: usize,
    pub log: BufWriter<File>,
    pub id: usize,
    pub python_numerator: Option<PythonNumerator>,
}

impl Topology {
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

        rotated_topology
    }
}

impl Integrand {
    pub fn new(topology: &Topology, settings: Settings, id: usize) -> Integrand {
        // create extra topologies with rotated kinematics to check the uncertainty
        let mut rng = rand::thread_rng();
        let mut topologies = vec![topology.clone()];
        for _ in 0..5 {
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

        let python_numerator = settings
            .general
            .python_numerator
            .as_ref()
            .map(|module| PythonNumerator::new(module, topology.n_loops));

        Integrand {
            topologies,
            cache_float: LTDCache::<float>::new(&topology),
            cache_f128: LTDCache::<f128>::new(&topology),
            running_max: Complex::default(),
            total_samples: 0,
            regular_point_count: 0,
            unstable_point_count: 0,
            unstable_f128_point_count: 0,
            nan_point_count: 0,
            log: BufWriter::new(
                File::create(format!("{}{}.log", settings.general.log_file_prefix, id))
                    .expect("Could not create log file"),
            ),
            settings,
            id,
            python_numerator,
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
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | unstable f128 points={} ({:.2}%)\n  | nan points={} ({:.2}%)",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_f128_point_count, self.unstable_f128_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100.);
        }

        writeln!(self.log,
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | unstable f128 points={} ({:.2}%)\n  | nan points={} ({:.2}%)",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_f128_point_count, self.unstable_f128_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100.).unwrap();
    }

    fn check_stability<T: FloatLike>(
        &self,
        x: &[f64],
        result: Complex<T>,
        cache: &mut LTDCache<T>,
    ) -> (T, T, Complex<T>) {
        if self.settings.general.numerical_instability_check {
            let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, result_rot) =
                self.topologies[1].evaluate(x, cache, &self.python_numerator);

            // compute the number of similar digits
            let (num, num_rot) =
                if self.settings.integrator.integrated_phase == IntegratedPhase::Imag {
                    (result.im, result_rot.im)
                } else {
                    (result.re, result_rot.re)
                };

            let d = if num.is_zero() && num_rot.is_zero() {
                T::from_usize(100).unwrap()
            } else {
                if num == -num_rot {
                    -T::from_usize(100).unwrap()
                } else {
                    -Float::abs((num - num_rot) / (num + num_rot)).log10()
                }
            };
            (d, Float::abs(num - num_rot), result_rot)
        } else {
            (
                T::from_f64(self.settings.general.relative_precision).unwrap(),
                T::from_f64(self.settings.general.absolute_precision).unwrap(),
                Complex::default(),
            )
        }
    }

    pub fn evaluate(&mut self, x: &[f64]) -> Complex<float> {
        if self.settings.general.integration_statistics
            && self.total_samples > 0
            && self.total_samples % self.settings.general.statistics_interval == 0
        {
            self.print_statistics();
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

        let mut cache_float = std::mem::replace(&mut self.cache_float, LTDCache::default());
        let (x, k_def, jac_para, jac_def, mut result) =
            self.topologies[0].evaluate(x, &mut cache_float, &self.python_numerator);
        let (d, diff, result_rot) = self.check_stability(x, result, &mut cache_float);
        std::mem::swap(&mut self.cache_float, &mut cache_float);

        self.total_samples += 1;

        if !result.is_finite()
            || !result_rot.is_finite()
            || d < NumCast::from(self.settings.general.relative_precision).unwrap()
            || diff > NumCast::from(self.settings.general.absolute_precision).unwrap()
        {
            if self.settings.general.integration_statistics {
                let loops = self.topologies[0].n_loops;
                self.print_info(
                    loops, false, true, x, &k_def, jac_para, jac_def, result, result_rot, d,
                );
            }

            // compute the point again with f128 to see if it is stable then
            let mut cache_f128 = std::mem::replace(&mut self.cache_f128, LTDCache::default());
            let (_, k_def_f128, jac_para_f128, jac_def_f128, result_f128) =
                self.topologies[0].evaluate(x, &mut cache_f128, &self.python_numerator);
            // NOTE: for this check we use a f64 rotation matrix at the moment!
            let (d_f128, diff_f128, result_rot_f128) =
                self.check_stability(x, result_f128, &mut cache_f128);
            std::mem::swap(&mut self.cache_f128, &mut cache_f128);

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

            if d_comp + float::from_f64(2.).unwrap() < d
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
                        topo.evaluate(x, &mut self.cache_float, &self.python_numerator);
                    eprintln!("  | rot{} ={:e}", i + 2, result_rot);
                }
            }

            self.unstable_point_count += 1;

            if !result_f128.is_finite()
                || !result_rot_f128.is_finite()
                || d_f128 < NumCast::from(self.settings.general.relative_precision).unwrap()
                || diff_f128 > NumCast::from(self.settings.general.absolute_precision).unwrap()
            {
                if self.settings.general.integration_statistics {
                    let loops = self.topologies[0].n_loops;
                    self.print_info(
                        loops,
                        false,
                        true,
                        x,
                        &k_def_f128,
                        jac_para_f128,
                        jac_def_f128,
                        result_f128,
                        result_rot_f128,
                        d_f128,
                    );
                }

                if !result_f128.is_finite() {
                    self.nan_point_count += 1;
                    return Complex::default();
                } else {
                    self.unstable_f128_point_count += 1;

                    if self.settings.general.return_unstable_point {
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
            let loops = self.topologies[0].n_loops;
            self.print_info(
                loops, new_max, false, x, &k_def, jac_para, jac_def, result, result_rot, d,
            );
        }

        result
    }
}
