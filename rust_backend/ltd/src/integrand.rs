use arrayvec::ArrayVec;
use topologies::Topology;
use vector::LorentzVector;
use Complex;

const MAX_LOOP: usize = 3;

use IntegratedPhase;
use Settings;

/// A structure that integrates and keeps statistics
#[derive(Debug, Clone)]
pub struct Integrand {
    pub settings: Settings,
    pub topologies: Vec<Topology>,
    pub running_max: Complex,
    pub total_samples: usize,
    pub nan_point_count: usize,
    pub unstable_point_count: usize,
    pub regular_point_count: usize,
}

impl Integrand {
    pub fn new(topology: &Topology, settings: Settings) -> Integrand {
        // create an extra topology with rotated kinematics to check the uncertainty
        let rot_matrix = [
            [0.985334, -0.141122, 0.0959282],
            [0.170454, 0.840137, -0.514893],
            [-0.00793036, 0.523693, 0.85187],
        ];
        let mut rotated_topology = topology.clone();
        rotated_topology.name = rotated_topology.name + "_rot";
        rotated_topology.rotation_matrix = rot_matrix.clone();

        for x in &mut rotated_topology.external_kinematics {
            let old_x = x.x;
            let old_y = x.y;
            let old_z = x.z;
            x.x = rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z;
            x.y = rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z;
            x.z = rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z;
        }

        for ll in &mut rotated_topology.loop_lines {
            for p in &mut ll.propagators {
                let old_x = p.q.x;
                let old_y = p.q.y;
                let old_z = p.q.z;
                p.q.x =
                    rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z;
                p.q.y =
                    rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z;
                p.q.z =
                    rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z;
            }
        }

        Integrand {
            topologies: vec![topology.clone(), rotated_topology],
            running_max: Complex::new(0., 0.),
            total_samples: 0,
            regular_point_count: 0,
            unstable_point_count: 0,
            nan_point_count: 0,
            settings,
        }
    }

    fn print_info(
        &self,
        n_loops: usize,
        new_max: bool,
        x: &[f64],
        k_def: ArrayVec<[LorentzVector<Complex>; MAX_LOOP]>,
        jac_para: f64,
        jac_def: Complex,
        result: Complex,
        rot_result: Complex,
        stable_digits: f64,
    ) {
        if new_max || !result.is_finite() || self.settings.general.debug > 0 {
            let sample_or_max = if new_max { "MAX" } else { "Sample" };
            match n_loops {
                1 => {
                    println!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], jac_para, jac_def
                    );
                }
                2 => {
                    println!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], jac_para, jac_def
                    );
                }
                3 => {
                    println!(
                        "{}\n  | result={:e}, rot={:e}, stable digits={}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        sample_or_max, result, rot_result, stable_digits, x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    );
                }
                _ => {}
            }
        }
    }

    fn print_statistics(&self) {
        println!("Statistics\n  | running max={:e}\n  | total samples={}\n  | regular points={} ({:.2}%)\n  | unstable points={} ({:.2}%)\n  | nan points={} ({:.2}%)",
            self.running_max, self.total_samples, self.regular_point_count, self.regular_point_count as f64 / self.total_samples as f64 * 100.,
            self.unstable_point_count, self.unstable_point_count as f64 / self.total_samples as f64 * 100.,
            self.nan_point_count, self.nan_point_count as f64 / self.total_samples as f64 * 100.);
    }

    pub fn evaluate(&mut self, x: &[f64]) -> Complex {
        if self.settings.general.integration_statistics
            && self.total_samples > 0
            && self.total_samples % self.settings.general.statistics_interval == 0
        {
            self.print_statistics();
        }

        let (x, k_def, jac_para, jac_def, result) = self.topologies[0].evaluate(x);
        let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, result_rot) =
            self.topologies[1].evaluate(x);
        self.total_samples += 1;

        // compute the number of similar digits
        // for now, only the real part
        let d = -((result.re - result_rot.re) / (result.re + result_rot.re))
            .abs()
            .log10();

        // FIXME: only checking the real for now
        if !result.re.is_finite()
            || !result_rot.re.is_finite()
            || d < self.settings.general.relative_precision
        {
            if self.settings.general.integration_statistics {
                self.print_info(
                    self.topologies[0].n_loops,
                    false,
                    x,
                    k_def,
                    jac_para,
                    jac_def,
                    result,
                    result_rot,
                    d,
                );
            }

            if !result.is_finite() {
                self.nan_point_count += 1;
            } else {
                self.unstable_point_count += 1;
            }

            // if we have large numerical instability, we return 0
            return Complex::new(0., 0.);
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
            self.print_info(
                self.topologies[0].n_loops,
                new_max,
                x,
                k_def,
                jac_para,
                jac_def,
                result,
                result_rot,
                d,
            );
        }

        result
    }
}
