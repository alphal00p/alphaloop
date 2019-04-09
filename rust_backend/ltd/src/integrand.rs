use arrayvec::ArrayVec;
use num_traits::NumCast;
use num_traits::{Float, FromPrimitive, Inv, One, ToPrimitive, Zero};
use std::fs::File;
use std::io::{BufWriter, Write};
use topologies::Topology;
use vector::LorentzVector;
use {float, Complex};

const MAX_LOOP: usize = 3;

use IntegratedPhase;
use Settings;

/// A structure that integrates and keeps statistics
pub struct Integrand {
    pub settings: Settings,
    pub topologies: Vec<Topology>,
    pub running_max: Complex,
    pub total_samples: usize,
    pub nan_point_count: usize,
    pub unstable_point_count: usize,
    pub regular_point_count: usize,
    pub log: BufWriter<File>,
}

impl Integrand {
    pub fn new(topology: &Topology, settings: Settings, id: usize) -> Integrand {
        // create an extra topology with rotated kinematics to check the uncertainty
        let angle = float::from_f64(0.33).unwrap(); // rotation angle
        let rv = (
            float::from_f64(2.).unwrap().sqrt().inv(),
            float::from_f64(3.).unwrap().sqrt().inv(),
            float::from_f64(6.).unwrap().sqrt().inv(),
        ); // rotation axis

        let cos_t = angle.cos();
        let sin_t = angle.sin();
        let cos_t_bar = float::one() - angle.cos();

        let rot_matrix: [[float; 3]; 3] = [
            [
                cos_t + rv.0 * rv.0 * cos_t_bar,
                rv.0 * rv.1 * cos_t_bar - rv.2 * sin_t,
                rv.0 * rv.2 * cos_t_bar + rv.1 * sin_t,
            ],
            [
                rv.0 * rv.1 * cos_t_bar + rv.2 * sin_t,
                cos_t + rv.1 * rv.1 * cos_t_bar,
                rv.1 * rv.2 * cos_t_bar - rv.0 * sin_t,
            ],
            [
                rv.0 * rv.2 * cos_t_bar - rv.1 * sin_t,
                rv.1 * rv.2 * cos_t_bar + rv.0 * sin_t,
                cos_t + rv.2 * rv.2 * cos_t_bar,
            ],
        ];

        let mut rotated_topology = topology.clone();
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

        Integrand {
            topologies: vec![topology.clone(), rotated_topology],
            running_max: Complex::default(),
            total_samples: 0,
            regular_point_count: 0,
            unstable_point_count: 0,
            nan_point_count: 0,
            log: BufWriter::new(
                File::create(format!("{}{}.log", settings.general.log_file_prefix, id))
                    .expect("Could not create log file"),
            ),
            settings,
        }
    }

    fn print_info(
        &mut self,
        n_loops: usize,
        new_max: bool,
        unstable: bool,
        x: &[f64],
        k_def: ArrayVec<[LorentzVector<Complex>; MAX_LOOP]>,
        jac_para: float,
        jac_def: Complex,
        result: Complex,
        rot_result: Complex,
        stable_digits: float,
    ) {
        if new_max || unstable || !result.is_finite() || self.settings.general.debug > 0 {
            let sample_or_max = if new_max { "MAX" } else { "Sample" };
            match n_loops {
                1 => {
                    if !unstable && self.settings.general.log_to_screen {
                        eprintln!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], jac_para, jac_def
                    );
                    }
                    writeln!(self.log, "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], jac_para, jac_def).unwrap();
                }
                2 => {
                    if !unstable && self.settings.general.log_to_screen {
                        eprintln!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], jac_para, jac_def
                    );
                    }
                    writeln!(self.log, "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], jac_para, jac_def).unwrap();
                }
                3 => {
                    if !unstable && self.settings.general.log_to_screen {
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
                _ => {}
            }
        }
    }

    fn print_statistics(&mut self) {
        if self.settings.general.log_to_screen {
            eprintln!(
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | nan points={} ({:.2}%)",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100.);
        }

        writeln!(self.log,
            "Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | nan points={} ({:.2}%)",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100.).unwrap();
    }

    pub fn evaluate(&mut self, x: &[f64]) -> Complex {
        if self.settings.general.integration_statistics
            && self.total_samples > 0
            && self.total_samples % self.settings.general.statistics_interval == 0
        {
            self.print_statistics();
        }

        let (x, k_def, jac_para, jac_def, result) = self.topologies[0].evaluate(x);
        self.total_samples += 1;

        let (d, result_rot) = if self.settings.general.numerical_instability_check {
            let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, result_rot) =
                self.topologies[1].evaluate(x);

            // compute the number of similar digits
            let (num, num_rot) =
                if self.settings.integrator.integrated_phase == IntegratedPhase::Imag {
                    (result.im, result_rot.im)
                } else {
                    (result.re, result_rot.re)
                };

            let d = if num.is_zero() && num_rot.is_zero() {
                float::from_usize(100).unwrap()
            } else {
                if num == -num_rot {
                    -float::from_usize(100).unwrap()
                } else {
                    -((num - num_rot) / (num + num_rot)).abs().log10()
                }
            };
            (d, result_rot)
        } else {
            (
                float::from_f64(self.settings.general.relative_precision).unwrap(),
                Complex::default(),
            )
        };

        if !result.is_finite()
            || !result_rot.is_finite()
            || d < NumCast::from(self.settings.general.relative_precision).unwrap()
        {
            if self.settings.general.integration_statistics {
                let loops = self.topologies[0].n_loops;
                self.print_info(
                    loops, false, true, x, k_def, jac_para, jac_def, result, result_rot, d,
                );
            }

            if !result.is_finite() {
                self.nan_point_count += 1;
            } else {
                self.unstable_point_count += 1;
            }

            // if we have large numerical instability, we return 0
            return Complex::default();
        }

        self.regular_point_count += 1;

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
                loops, new_max, false, x, k_def, jac_para, jac_def, result, result_rot, d,
            );
        }

        result
    }
}
