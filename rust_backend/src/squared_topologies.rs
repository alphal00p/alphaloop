use crate::dashboard::{StatusUpdate, StatusUpdateSender};
use crate::dualklt3::Dualklt3;
use crate::dualkt2::Dualkt2;
use crate::dualkt3::Dualkt3;
use crate::dualt2::Dualt2;
use crate::dualt3::Dualt3;
use crate::integrand::{IntegrandImplementation, IntegrandSample};
use crate::observables::EventManager;
use crate::topologies::{FixedDeformationLimit, LTDCache, SOCPProblem, Topology};
use crate::utils;
use crate::IntegratedPhase;
use crate::{DeformationStrategy, FloatLike, IntegrandType, NormalisingFunction, Settings};
use arrayvec::ArrayVec;
use color_eyre::{Help, Report};
use dlopen::raw::Library;
use eyre::WrapErr;
use f128::f128;
use havana::{ContinuousGrid, DiscreteGrid, Grid};
use hyperdual::Hyperdual;
use itertools::Itertools;
use libc::{c_double, c_int};
use lorentz_vector::{Field, LorentzVector, RealNumberLike};
use nalgebra::Scalar;
use num::Complex;
use num_traits::{Float, FloatConst, FromPrimitive, Inv, One, ToPrimitive, Zero};
use rand::Rng;
use serde::Deserialize;
use smallvec::SmallVec;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::time::Instant;
use utils::Signum;

#[cfg(not(feature = "higher_loops"))]
pub const MAX_SG_LOOP: usize = 5;
#[cfg(feature = "higher_loops")]
pub const MAX_SG_LOOP: usize = 10;

pub const MAX_DUAL_SIZE: usize = 12;

mod form_integrand {
    use super::FORMIntegrandCallSignature;

    pub trait GetIntegrand {
        fn get_integrand_ltd(
            api_container: &FORMIntegrandCallSignature,
            p: &[Self],
            params: &[Self],
            conf: usize,
            res: &mut [Self],
        ) where
            Self: std::marker::Sized;

        fn get_integrand_pf(
            api_container: &FORMIntegrandCallSignature,
            p: &[Self],
            params: &[Self],
            conf: usize,
            res: &mut [Self],
        ) where
            Self: std::marker::Sized;

        fn get_integrand_ltd_mpfr(
            api_container: &FORMIntegrandCallSignature,
            p: &[Self],
            params: &[Self],
            conf: usize,
            prec: usize,
            res: &mut [Self],
        ) where
            Self: std::marker::Sized;

        fn get_integrand_pf_mpfr(
            api_container: &FORMIntegrandCallSignature,
            p: &[Self],
            params: &[Self],
            conf: usize,
            prec: usize,
            res: &mut [Self],
        ) where
            Self: std::marker::Sized;
    }

    impl GetIntegrand for f64 {
        fn get_integrand_ltd(
            api_container: &FORMIntegrandCallSignature,
            p: &[f64],
            params: &[f64],
            conf: usize,
            res: &mut [f64],
        ) {
            unsafe {
                if let Some(eval) = api_container.evaluate_ltd {
                    eval(
                        &p[0] as *const f64,
                        &params[0] as *const f64,
                        conf as i32,
                        &mut res[0] as *mut f64,
                    );
                }
            }
        }

        fn get_integrand_pf(
            api_container: &FORMIntegrandCallSignature,
            p: &[f64],
            params: &[f64],
            conf: usize,
            res: &mut [f64],
        ) {
            unsafe {
                if let Some(eval) = api_container.evaluate_pf {
                    eval(
                        &p[0] as *const f64,
                        &params[0] as *const f64,
                        conf as i32,
                        &mut res[0] as *mut f64,
                    );
                }
            }
        }

        fn get_integrand_ltd_mpfr(
            _api_container: &FORMIntegrandCallSignature,
            _p: &[f64],
            _params: &[f64],
            _conf: usize,
            _prec: usize,
            _res: &mut [f64],
        ) {
            unimplemented!("MPFR is only supported in f128 mode")
        }

        fn get_integrand_pf_mpfr(
            _api_container: &FORMIntegrandCallSignature,
            _p: &[f64],
            _params: &[f64],
            _conf: usize,
            _prec: usize,
            _res: &mut [f64],
        ) {
            unimplemented!("MPFR is only supported in f128 mode")
        }
    }

    impl GetIntegrand for f128::f128 {
        fn get_integrand_ltd(
            api_container: &FORMIntegrandCallSignature,
            p: &[f128::f128],
            params: &[f128::f128],
            conf: usize,
            res: &mut [f128::f128],
        ) {
            unsafe {
                if let Some(eval) = api_container.evaluate_ltd_f128 {
                    eval(
                        &p[0] as *const f128::f128,
                        &params[0] as *const f128::f128,
                        conf as i32,
                        &mut res[0] as *mut f128::f128,
                    );
                }
            }
        }

        fn get_integrand_pf(
            api_container: &FORMIntegrandCallSignature,
            p: &[f128::f128],
            params: &[f128::f128],
            conf: usize,
            res: &mut [f128::f128],
        ) {
            unsafe {
                if let Some(eval) = api_container.evaluate_pf_f128 {
                    eval(
                        &p[0] as *const f128::f128,
                        &params[0] as *const f128::f128,
                        conf as libc::c_int,
                        &mut res[0] as *mut f128::f128,
                    );
                }
            }
        }

        fn get_integrand_ltd_mpfr(
            api_container: &FORMIntegrandCallSignature,
            p: &[f128::f128],
            params: &[f128::f128],
            conf: usize,
            prec: usize,
            res: &mut [f128::f128],
        ) {
            unsafe {
                if let Some(eval) = api_container.evaluate_ltd_mpfr {
                    eval(
                        &p[0] as *const f128::f128,
                        &params[0] as *const f128::f128,
                        conf as i32,
                        prec as i32,
                        &mut res[0] as *mut f128::f128,
                    );
                }
            }
        }

        fn get_integrand_pf_mpfr(
            api_container: &FORMIntegrandCallSignature,
            p: &[f128::f128],
            params: &[f128::f128],
            conf: usize,
            prec: usize,
            res: &mut [f128::f128],
        ) {
            unsafe {
                if let Some(eval) = api_container.evaluate_pf_mpfr {
                    eval(
                        &p[0] as *const f128::f128,
                        &params[0] as *const f128::f128,
                        conf as i32,
                        prec as i32,
                        &mut res[0] as *mut f128::f128,
                    );
                }
            }
        }
    }
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub struct CutkoskyCut {
    pub name: String,
    pub sign: i8,
    pub power: usize,
    pub particle_id: isize,
    pub signature: (Vec<i8>, Vec<i8>),
    pub m_squared: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCutDiagramInfo {
    pub graph: Topology,
    pub conjugate_deformation: bool,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCutDiagramSet {
    pub id: usize,
    pub diagram_info: Vec<CutkoskyCutDiagramInfo>,
    pub cb_to_lmb: Option<Vec<i8>>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCuts {
    pub cuts: Vec<CutkoskyCut>,
    pub diagram_sets: Vec<CutkoskyCutDiagramSet>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MGNumeratorCallSignature {
    pub left_diagram_id: usize,
    pub proc_id: usize,
    pub right_diagram_id: usize,
}

#[derive(Debug, Clone, Deserialize)]
pub struct FORMNumeratorCallSignature {
    pub id: usize,
}

#[derive(Deserialize)]
pub struct FORMIntegrandCallSignature {
    #[serde(default)]
    pub base_path: String,
    pub id: usize,
    #[serde(skip_deserializing)]
    pub form_integrand: Option<Library>,
    #[serde(skip_deserializing)]
    pub evaluate_pf:
        Option<unsafe extern "C" fn(*const c_double, *const c_double, c_int, *mut c_double)>,
    #[serde(skip_deserializing)]
    pub evaluate_ltd:
        Option<unsafe extern "C" fn(*const c_double, *const c_double, c_int, *mut c_double)>,
    #[serde(skip_deserializing)]
    pub evaluate_pf_f128: Option<unsafe extern "C" fn(*const f128, *const f128, c_int, *mut f128)>,
    #[serde(skip_deserializing)]
    pub evaluate_ltd_f128: Option<unsafe extern "C" fn(*const f128, *const f128, c_int, *mut f128)>,
    #[serde(skip_deserializing)]
    pub evaluate_ltd_mpfr:
        Option<unsafe extern "C" fn(*const f128, *const f128, c_int, c_int, *mut f128)>,
    #[serde(skip_deserializing)]
    pub evaluate_pf_mpfr:
        Option<unsafe extern "C" fn(*const f128, *const f128, c_int, c_int, *mut f128)>,
    #[serde(default)]
    pub prec: usize, // will be set by the precision checker
    #[serde(default)]
    pub extra_calls: Vec<(usize, usize)>,
    #[serde(skip_deserializing)]
    pub extra_form_integrands: Vec<Library>,
}

impl FORMIntegrandCallSignature {
    pub fn initialise(&mut self, base_path: String) {
        self.base_path = base_path;

        let a = format!("{}/lib/libFORM_sg_{}.so", self.base_path, self.id);
        self.form_integrand =
            Some(Library::open(a).expect("Could not open library or load symbols"));

        self.evaluate_pf = Some(
            unsafe {
                self.form_integrand
                    .as_ref()
                    .unwrap()
                    .symbol(&format!("evaluate_PF_{}", self.id))
            }
            .unwrap(),
        );

        self.evaluate_pf_f128 = unsafe {
            self.form_integrand
                .as_ref()
                .unwrap()
                .symbol(&format!("evaluate_PF_{}_f128", self.id))
        }
        .ok();

        self.evaluate_ltd_mpfr = unsafe {
            self.form_integrand
                .as_ref()
                .unwrap()
                .symbol(&format!("evaluate_LTD_{}_mpfr", self.id))
        }
        .ok();

        self.evaluate_pf_mpfr = unsafe {
            self.form_integrand
                .as_ref()
                .unwrap()
                .symbol(&format!("evaluate_PF_{}_mpfr", self.id))
        }
        .ok();

        self.evaluate_ltd = unsafe {
            self.form_integrand
                .as_ref()
                .unwrap()
                .symbol(&format!("evaluate_LTD_{}", self.id))
        }
        .ok();

        self.evaluate_ltd_f128 = unsafe {
            self.form_integrand
                .as_ref()
                .unwrap()
                .symbol(&format!("evaluate_LTD_{}_f128", self.id))
        }
        .ok();

        self.extra_form_integrands = self
            .extra_calls
            .iter()
            .map(|(diag, _)| {
                Library::open(&format!("{}/lib/libFORM_sg_{}.so", self.base_path, diag))
                    .expect("Could not open library or load symbols")
            })
            .collect();
    }
}

impl Clone for FORMIntegrandCallSignature {
    fn clone(&self) -> Self {
        FORMIntegrandCallSignature {
            base_path: self.base_path.clone(),
            id: self.id,
            form_integrand: self.form_integrand.as_ref().map(|_| {
                Library::open(&format!("{}/lib/libFORM_sg_{}.so", self.base_path, self.id))
                    .expect("Could not open library or load symbols")
            }),
            evaluate_pf: self.evaluate_pf.clone(), // FIXME: is this safe?
            evaluate_ltd: self.evaluate_ltd.clone(),
            evaluate_pf_f128: self.evaluate_pf_f128.clone(),
            evaluate_ltd_f128: self.evaluate_ltd_f128.clone(),
            evaluate_ltd_mpfr: self.evaluate_ltd_mpfr.clone(),
            evaluate_pf_mpfr: self.evaluate_pf_mpfr.clone(),
            prec: self.prec,
            extra_calls: self.extra_calls.clone(),
            extra_form_integrands: self
                .extra_calls
                .iter()
                .map(|(diag, _)| {
                    Library::open(&format!("{}/lib/libFORM_sg_{}.so", self.base_path, diag))
                        .expect("Could not open library or load symbols")
                })
                .collect(),
        }
    }
}

#[derive(Clone, Deserialize)]
pub struct FORMIntegrand {
    call_signature: Option<FORMIntegrandCallSignature>,
}

#[derive(Clone, Deserialize)]
pub struct MultiChannelingBasis {
    pub channel_id: usize,
    pub cutkosky_cut_id: i8,
    pub defining_propagators: Vec<String>,
    pub defining_propagator_masses: Vec<f64>,
    pub signatures: Vec<(Vec<i8>, Vec<i8>)>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Propagator {
    pub name: String,
    #[serde(rename = "PDG")]
    pub pdg: isize,
    pub signature: (Vec<i8>, Vec<i8>),
    pub m_squared: f64,
    pub vertices: (usize, usize),
    pub power: usize,
}

#[derive(Clone, Deserialize)]
pub struct SquaredTopology {
    pub name: String,
    pub n_loops: usize,
    pub n_incoming_momenta: usize,
    #[serde(skip_deserializing)]
    pub e_cm_squared: f128,
    pub overall_numerator: f64,
    #[serde(skip_deserializing)]
    pub external_momenta: Vec<LorentzVector<f128>>,
    pub cutkosky_cuts: Vec<CutkoskyCuts>,
    pub analytical_result_real: Option<f64>,
    pub analytical_result_imag: Option<f64>,
    #[serde(skip_deserializing)]
    pub settings: Settings,
    #[serde(skip_deserializing)]
    pub rotation_matrix: [[f128; 3]; 3],
    pub default_fixed_cut_momenta: (Vec<LorentzVector<f64>>, Vec<LorentzVector<f64>>),
    #[serde(default)]
    pub multi_channeling_bases: Vec<MultiChannelingBasis>,
    #[serde(default)]
    pub multi_channeling_lmb_bases: Vec<MultiChannelingBasis>,
    #[serde(default)]
    pub optimal_channel_ids: Option<Vec<usize>>,
    #[serde(skip_deserializing)]
    pub multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f128>>, Vec<f128>)>,
    #[serde(rename = "FORM_integrand")]
    pub form_integrand: FORMIntegrand,
    #[serde(skip_deserializing)]
    pub is_stability_check_topo: bool,
    pub propagators: Vec<Propagator>,
}

#[derive(Clone)]
pub struct SquaredTopologySet {
    pub name: String,
    pub e_cm_squared: f128,
    pub topologies: Vec<SquaredTopology>,
    additional_topologies: Vec<Vec<(SquaredTopology, Vec<(Vec<i8>, Vec<i8>)>)>>,
    pub multiplicity: Vec<f64>,
    pub settings: Settings,
    pub rotation_matrix: [[f128; 3]; 3],
    pub multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f128>>, Vec<f128>)>,
    pub stability_check_topologies: Vec<Vec<SquaredTopology>>,
    pub is_stability_check_topo: bool,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetAdditionalTopology {
    pub name: String,
    pub this_lmb_to_defining_lmb: Vec<(Vec<i8>, Vec<i8>)>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetTopology {
    pub name: String,
    pub multiplicity: f64,
    #[serde(default, rename = "additional_LMBs")]
    pub additional_lmbs: Vec<SquaredTopologySetAdditionalTopology>,
    #[serde(default)]
    pub stability_check_topologies: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetInput {
    pub name: String,
    topologies: Vec<SquaredTopologySetTopology>,
}

impl SquaredTopologySet {
    pub fn from_one(mut squared_topology: SquaredTopology) -> SquaredTopologySet {
        squared_topology.generate_multi_channeling_channels();

        SquaredTopologySet {
            name: squared_topology.name.clone(),
            settings: squared_topology.settings.clone(),
            e_cm_squared: squared_topology.e_cm_squared,
            rotation_matrix: squared_topology.rotation_matrix.clone(),
            multi_channeling_channels: squared_topology.multi_channeling_channels.clone(),
            topologies: vec![squared_topology],
            additional_topologies: vec![vec![]],
            multiplicity: vec![1.],
            stability_check_topologies: vec![vec![]],
            is_stability_check_topo: false,
        }
    }

    pub fn from_file(filename: &str, settings: &Settings) -> Result<SquaredTopologySet, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open squared topology set file {}", filename))
            .suggestion("Does the path exist?")?;

        let squared_topology_set_input: SquaredTopologySetInput = serde_yaml::from_reader(f)
            .wrap_err("Could not parse squared topology set file")
            .suggestion("Is it a correct yaml file")?;

        let mut topologies: Vec<SquaredTopology> = vec![];
        let mut additional_topologies: Vec<Vec<(SquaredTopology, Vec<(Vec<i8>, Vec<i8>)>)>> =
            vec![];
        let mut multiplicity: Vec<f64> = vec![];
        let mut stability_topologies = vec![];

        for topo in squared_topology_set_input.topologies {
            let filename = std::path::Path::new(&filename)
                .with_file_name(topo.name)
                .with_extension("yaml");
            let squared_topology = SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                .wrap_err("Could not load subtopology file")?;

            let mut additional_topologies_for_topo = vec![];
            for t in topo.additional_lmbs {
                let filename = std::path::Path::new(&filename)
                    .with_file_name(t.name)
                    .with_extension("yaml");
                let additional_squared_topology =
                    SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                        .wrap_err("Could not load subtopology file")?;
                additional_topologies_for_topo
                    .push((additional_squared_topology, t.this_lmb_to_defining_lmb));
            }
            additional_topologies.push(additional_topologies_for_topo);

            let mut stability_topologies_for_topo = vec![];
            for t in topo.stability_check_topologies {
                let filename = std::path::Path::new(&filename)
                    .with_file_name(t)
                    .with_extension("yaml");
                let mut additional_squared_topology =
                    SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                        .wrap_err("Could not load subtopology file")?;
                additional_squared_topology.is_stability_check_topo = true;
                stability_topologies_for_topo.push(additional_squared_topology);
            }
            stability_topologies.push(stability_topologies_for_topo);

            topologies.push(squared_topology);
            multiplicity.push(topo.multiplicity);
        }

        let rotation_matrix = [
            [f128::one(), f128::zero(), f128::zero()],
            [f128::zero(), f128::one(), f128::zero()],
            [f128::zero(), f128::zero(), f128::one()],
        ];

        let mut sts = SquaredTopologySet {
            name: squared_topology_set_input.name,
            e_cm_squared: topologies[0].e_cm_squared,
            topologies,
            additional_topologies,
            rotation_matrix,
            settings: settings.clone(),
            multiplicity,
            multi_channeling_channels: vec![],
            stability_check_topologies: stability_topologies,
            is_stability_check_topo: false,
        };
        sts.create_multi_channeling_channels();
        Ok(sts)
    }

    fn create_multi_channeling_channels(&mut self) {
        let mut multi_channeling_channels = vec![];
        let mut max_loops = 0;
        for t in &mut self.topologies {
            t.generate_multi_channeling_channels();

            if t.n_loops > max_loops {
                multi_channeling_channels.clear();
                max_loops = t.n_loops;
            }
            if t.n_loops == max_loops {
                multi_channeling_channels.extend(t.multi_channeling_channels.iter().cloned());
            }
        }

        // filter duplicate channels between supergraphs
        let mut unique_multi_channeling_channels: Vec<(
            Vec<i8>,
            Vec<i8>,
            Vec<LorentzVector<f128>>,
            Vec<f128>,
        )> = vec![];

        'mc_loop: for c in multi_channeling_channels {
            for uc in &unique_multi_channeling_channels {
                if uc.0 == c.0 && uc.1 == c.1 {
                    // test for true equality
                    if uc
                        .2
                        .iter()
                        .zip(c.2.iter())
                        .all(|(s1, s2)| (s1 - s2).spatial_squared() < Into::<f128>::into(1.0e-15))
                    {
                        continue 'mc_loop;
                    }
                }
            }
            unique_multi_channeling_channels.push(c);
        }

        self.multi_channeling_channels = unique_multi_channeling_channels;
    }

    /// Create a rotated version of this squared topology. The axis needs to be normalized.
    fn rotate(&self, angle: f128, axis: (f128, f128, f128)) -> SquaredTopologySet {
        let rotated_topologies: Vec<_> = self
            .topologies
            .iter()
            .map(|t| {
                let mut stab_t = t.rotate(angle, axis);
                stab_t.is_stability_check_topo = true;
                stab_t
            })
            .collect();
        let mut c = self.multi_channeling_channels.clone();

        let rot_matrix = &rotated_topologies[0].rotation_matrix;
        for (_, _, shifts, _) in &mut c {
            for shift in shifts.iter_mut() {
                let old_shift = shift.clone();
                shift.x = rot_matrix[0][0] * old_shift.x
                    + rot_matrix[0][1] * old_shift.y
                    + rot_matrix[0][2] * old_shift.z;
                shift.y = rot_matrix[1][0] * old_shift.x
                    + rot_matrix[1][1] * old_shift.y
                    + rot_matrix[1][2] * old_shift.z;
                shift.z = rot_matrix[2][0] * old_shift.x
                    + rot_matrix[2][1] * old_shift.y
                    + rot_matrix[2][2] * old_shift.z;
            }
        }

        SquaredTopologySet {
            name: self.name.clone(),
            multiplicity: self.multiplicity.clone(),
            e_cm_squared: self.e_cm_squared,
            settings: self.settings.clone(),
            rotation_matrix: rotated_topologies[0].rotation_matrix.clone(),
            additional_topologies: vec![vec![]; rotated_topologies.len()], // rotated versions of additional topologies are not supported
            topologies: rotated_topologies,
            multi_channeling_channels: c,
            stability_check_topologies: vec![vec![]; self.topologies.len()],
            is_stability_check_topo: true,
        }
    }

    pub fn get_maximum_loop_count(&self) -> usize {
        self.topologies
            .iter()
            .map(|t| t.n_loops - self.settings.cross_section.fixed_cut_momenta.len())
            .max()
            .unwrap()
    }

    pub fn create_caches<T: FloatLike>(&self) -> SquaredTopologyCache<T> {
        SquaredTopologyCache {
            topology_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
            deformation_vector_cache: vec![],
            deformation_jacobian_matrix: vec![],
            momenta: vec![],
            params: vec![],
            current_supergraph: 0,
            current_deformation_index: 0,
        }
    }

    pub fn print_info(&self, status_update_sender: &mut StatusUpdateSender) {
        /*
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Overall number of loops in all supergraphs: {}",
                self.n_loops
            )))
            .unwrap();
        */
        let n_topologies = self.topologies.len();

        let mut n_cutkosky_cuts_per_cut_cardinality = HashMap::new();
        let mut n_diagrams_per_loop = HashMap::new();
        let mut n_ltd_cuts_per_loop = HashMap::new();

        for topology in &self.topologies {
            for cutkosky_cuts in &topology.cutkosky_cuts {
                *n_cutkosky_cuts_per_cut_cardinality
                    .entry(cutkosky_cuts.cuts.len())
                    .or_insert(0) += 1;
                for diagram_set in &cutkosky_cuts.diagram_sets {
                    for d_info in &diagram_set.diagram_info {
                        *n_diagrams_per_loop.entry(d_info.graph.n_loops).or_insert(0) += 1;
                        *n_ltd_cuts_per_loop.entry(d_info.graph.n_loops).or_insert(0) += d_info
                            .graph
                            .ltd_cut_options
                            .iter()
                            .map(|v| v.len())
                            .sum::<usize>()
                            .max(1);
                    }
                }
            }
        }

        let mut tmp_vec: Vec<_> = n_ltd_cuts_per_loop.into_iter().collect();
        let tmp_sum: usize = tmp_vec.iter().map(|(_, v)| *v).sum();
        tmp_vec.sort_by(|a, b| a.0.cmp(&b.0));
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Loop count -> number of LTD cuts: {} total: {}",
                tmp_vec
                    .iter()
                    .map(|(k, v)| format!("{}->{}", k, v))
                    .collect::<Vec<String>>()
                    .join(", "),
                tmp_sum
            )))
            .unwrap();

        let mut tmp_vec: Vec<_> = n_diagrams_per_loop.into_iter().collect();
        let tmp_sum: usize = tmp_vec.iter().map(|(_, v)| *v).sum();
        tmp_vec.sort_by(|a, b| a.0.cmp(&b.0));
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Loop count -> number of subdiagrams: {} total: {}",
                tmp_vec
                    .iter()
                    .map(|(k, v)| format!("{}->{}", k, v))
                    .collect::<Vec<String>>()
                    .join(", "),
                tmp_sum
            )))
            .unwrap();

        let mut tmp_vec: Vec<(usize, usize)> =
            n_cutkosky_cuts_per_cut_cardinality.into_iter().collect();
        let tmp_sum: usize = tmp_vec.iter().map(|(_, v)| *v).sum();
        tmp_vec.sort_by(|a, b| a.0.cmp(&b.0));
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Cut cardinality -> number of Cutkosky cuts: {} total: {}",
                tmp_vec
                    .iter()
                    .map(|(k, v)| format!("{}->{}", k, v))
                    .collect::<Vec<String>>()
                    .join(", "),
                tmp_sum
            )))
            .unwrap();

        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Number of integration channels: {}",
                self.multi_channeling_channels.len()
            )))
            .unwrap();

        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Number of supergraphs: {}",
                n_topologies
            )))
            .unwrap();

        if n_topologies == 1 {
            if let Some(imag_res) = self.topologies[0].analytical_result_imag {
                if imag_res != 0.0 {
                    status_update_sender
                        .send(StatusUpdate::Message(format!(
                            "Target benchmark result (imag) : {:+e}",
                            imag_res
                        )))
                        .unwrap();
                }
            }
            if let Some(real_res) = self.topologies[0].analytical_result_real {
                if real_res != 0.0 {
                    status_update_sender
                        .send(StatusUpdate::Message(format!(
                            "Target benchmark result (real) : {:+e}",
                            real_res
                        )))
                        .unwrap();
                }
            }
        }
    }

    pub fn multi_channeling<'a, T: form_integrand::GetIntegrand + FloatLike>(
        &mut self,
        selected_topology: Option<usize>,
        sample: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCache<T>,
        mut event_manager: Option<&mut EventManager>,
    ) -> Complex<T> {
        let mut xrot = [T::zero(); MAX_SG_LOOP * 3];

        // obtain the sampled channel if one is provided
        let (selected_channel, x_orig) = match sample {
            IntegrandSample::Flat(_, x) => (None, x),
            IntegrandSample::Nested(x) => match x {
                havana::Sample::ContinuousGrid(_, x) => (None, x.as_slice()),
                havana::Sample::DiscreteGrid(_, selected_channel, sub_sample) => (
                    Some(selected_channel[0]),
                    match sub_sample.as_ref().unwrap().as_ref() {
                        havana::Sample::ContinuousGrid(_, x) => x.as_slice(),
                        _ => unreachable!(),
                    },
                ),
                havana::Sample::MultiChannel(_, _, _) => unimplemented!(),
            },
        };

        let x: ArrayVec<[T; MAX_SG_LOOP * 3]> =
            x_orig.iter().map(|x| Into::<T>::into(*x)).collect();

        let mut mc_channels = if let Some(t) = selected_topology {
            std::mem::replace(&mut self.topologies[t].multi_channeling_channels, vec![])
        } else {
            std::mem::replace(&mut self.multi_channeling_channels, vec![])
        };

        // paramaterize and consider the result in a channel basis
        let n_fixed = self.settings.cross_section.fixed_cut_momenta.len();
        let n_loops = x.len() / 3 + n_fixed;
        let mut k_channel = [LorentzVector::default(); MAX_SG_LOOP];
        for i in 0..n_loops {
            xrot[..x.len()].copy_from_slice(&x);

            let (l_energy, l_space) = if i < n_loops - n_fixed {
                // set the loop index to i + 1 so that we can also shift k
                (
                    T::zero(),
                    Topology::parameterize(
                        &xrot[i * 3..(i + 1) * 3],
                        T::convert_from(&self.e_cm_squared),
                        i,
                        &self.settings,
                    )
                    .0,
                )
            } else {
                let m = self.settings.cross_section.fixed_cut_momenta[i + n_fixed - n_loops].cast();
                (m.t, [m.x, m.y, m.z])
            };

            let rot = &self.rotation_matrix;
            k_channel[i] = LorentzVector::from_args(
                l_energy,
                T::convert_from(&rot[0][0]) * l_space[0]
                    + T::convert_from(&rot[0][1]) * l_space[1]
                    + T::convert_from(&rot[0][2]) * l_space[2],
                T::convert_from(&rot[1][0]) * l_space[0]
                    + T::convert_from(&rot[1][1]) * l_space[1]
                    + T::convert_from(&rot[1][2]) * l_space[2],
                T::convert_from(&rot[2][0]) * l_space[0]
                    + T::convert_from(&rot[2][1]) * l_space[1]
                    + T::convert_from(&rot[2][2]) * l_space[2],
            );
        }

        let mut k_lmb = [LorentzVector::default(); MAX_SG_LOOP];
        let mut k_other_channel = [LorentzVector::default(); MAX_SG_LOOP];
        let mut event_counter = 0;
        let mut result = Complex::zero();
        for (channel_id, (channel, _, channel_shift, channel_masses)) in
            mc_channels.iter().enumerate()
        {
            if let Some(selected_channel) = selected_channel {
                if selected_channel != channel_id {
                    continue;
                }
            }

            if let Some(selected_channel) = self.settings.general.multi_channeling_channel {
                if selected_channel != channel_id as isize {
                    continue;
                }
            }

            // transform to the loop momentum basis
            for (kk, r) in utils::zip_eq(&mut k_lmb[..n_loops], channel.chunks(n_loops)) {
                *kk = LorentzVector::default();

                for ((ko, s), shift) in
                    utils::zip_eq(utils::zip_eq(&k_channel[..n_loops], r), channel_shift)
                {
                    *kk += (ko - shift.cast()).multiply_sign(*s);
                }
            }

            // determine the normalization constant
            let mut normalization = T::zero();
            for (_, other_channel_inv, other_channel_shift, other_channel_masses) in &mc_channels {
                // transform from the loop momentum basis to the other channel basis
                for ((kk, r), shift) in utils::zip_eq(
                    utils::zip_eq(
                        &mut k_other_channel[..n_loops],
                        other_channel_inv.chunks(n_loops),
                    ),
                    other_channel_shift,
                ) {
                    *kk = shift.cast();
                    for (ko, s) in utils::zip_eq(&k_lmb[..n_loops], r) {
                        *kk += ko.multiply_sign(*s);
                    }
                }

                let mut inv_jac_para = T::one();
                for i in 0..n_loops {
                    let jac = if i < n_loops - n_fixed {
                        let rot = &self.rotation_matrix;

                        // undo the rotation
                        let rotated = LorentzVector::from_args(
                            T::zero(),
                            T::convert_from(&rot[0][0]) * k_other_channel[i].x
                                + T::convert_from(&rot[1][0]) * k_other_channel[i].y
                                + T::convert_from(&rot[2][0]) * k_other_channel[i].z,
                            T::convert_from(&rot[0][1]) * k_other_channel[i].x
                                + T::convert_from(&rot[1][1]) * k_other_channel[i].y
                                + T::convert_from(&rot[2][1]) * k_other_channel[i].z,
                            T::convert_from(&rot[0][2]) * k_other_channel[i].x
                                + T::convert_from(&rot[1][2]) * k_other_channel[i].y
                                + T::convert_from(&rot[2][2]) * k_other_channel[i].z,
                        );
                        if self.settings.general.multi_channeling_alpha < 0. {
                            Topology::inv_parametrize(
                                &rotated,
                                T::convert_from(&self.e_cm_squared),
                                i,
                                &self.settings,
                            )
                            .1
                        } else {
                            T::powf(
                                (rotated.spatial_squared()
                                    + T::convert_from(&other_channel_masses[i])
                                        * T::convert_from(&other_channel_masses[i]))
                                .sqrt(),
                                Into::<T>::into(self.settings.general.multi_channeling_alpha),
                            )
                            .inv()
                        }
                    } else {
                        if self.settings.general.multi_channeling_alpha < 0. {
                            (Into::<T>::into(2.) * <T as FloatConst>::PI())
                                .powi(4)
                                .inv()
                        } else {
                            Into::<T>::into(1.)
                        }
                    };
                    inv_jac_para *= jac;
                }

                normalization += inv_jac_para;
            }

            // evaluate the integrand in this channel
            let mut channel_result = Complex::zero();
            for (current_supergraph, (t, &m)) in
                utils::zip_eq(&mut self.topologies, &self.multiplicity).enumerate()
            {
                if let Some(tt) = selected_topology {
                    if current_supergraph != tt {
                        continue;
                    }
                }

                cache.current_supergraph = current_supergraph;

                if self.settings.general.debug > 0 {
                    println!(
                        "Evaluating supergraph {} for channel {}",
                        t.name, channel_id
                    );

                    print!("Point in channel: ");
                    for i in 0..n_loops {
                        let (x, _) = Topology::inv_parametrize(
                            &k_lmb[i],
                            T::convert_from(&self.e_cm_squared),
                            i,
                            &self.settings,
                        );
                        if i > 0 {
                            print!(", {}", x.iter().join(","));
                        } else {
                            print!("{}", x.iter().join(","));
                        }
                    }
                    println!("");
                }

                // undo the jacobian for unused dimensions
                let mut jac_correction = T::one();
                for i in t.n_loops..n_loops {
                    let (_, jac) = Topology::inv_parametrize(
                        &k_lmb[i],
                        T::convert_from(&self.e_cm_squared),
                        i,
                        &self.settings,
                    );
                    jac_correction *= jac;
                }
                if self.settings.general.multi_channeling_alpha >= 0. {
                    for i in 0..t.n_loops {
                        let (_, jac) = Topology::inv_parametrize(
                            &k_channel[i],
                            T::convert_from(&self.e_cm_squared),
                            i,
                            &self.settings,
                        );
                        jac_correction *= jac.inv();
                    }
                    for i in 0..t.n_loops {
                        jac_correction *= T::powf(
                            (k_channel[i].spatial_squared()
                                + T::convert_from(&channel_masses[i])
                                    * T::convert_from(&channel_masses[i]))
                            .sqrt(),
                            Into::<T>::into(self.settings.general.multi_channeling_alpha),
                        )
                        .inv();
                    }
                    for _i in t.n_loops..n_loops {
                        jac_correction *= (Into::<T>::into(2.) * <T as FloatConst>::PI()).powi(4);
                    }
                }

                channel_result += t.evaluate_mom(&k_lmb[..t.n_loops], cache, &mut event_manager)
                    / normalization
                    * jac_correction
                    * Into::<T>::into(m as f64);
                if let Some(em) = &mut event_manager {
                    if em.track_events {
                        for e in em.event_buffer[event_counter..].iter_mut() {
                            e.integrand *= m as f64 / normalization.to_f64().unwrap()
                                * jac_correction.to_f64().unwrap();
                        }
                        event_counter = em.event_buffer.len();
                    }
                }
            }

            result += channel_result;
        }

        if let Some(t) = selected_topology {
            std::mem::swap(
                &mut self.topologies[t].multi_channeling_channels,
                &mut mc_channels,
            );
        } else {
            std::mem::swap(&mut self.multi_channeling_channels, &mut mc_channels);
        };

        result
    }

    pub fn evaluate<'a, T: form_integrand::GetIntegrand + FloatLike>(
        &mut self,
        sample: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCache<T>,
        mut event_manager: Option<&mut EventManager>,
    ) -> Complex<T> {
        if self.settings.general.debug > 0 {
            println!("x-space point: {:?}", sample);
        }

        // obtain the sampled topology if one is provided
        let (selected_topology, sub_sample) = if self.topologies.len() > 1 {
            match sample {
                IntegrandSample::Flat(w, x) => (None, IntegrandSample::Flat(w, x)),
                IntegrandSample::Nested(x) => match x {
                    havana::Sample::ContinuousGrid(_, _) => (None, IntegrandSample::Nested(x)),
                    havana::Sample::DiscreteGrid(_, selected_topology, sub_sample) => (
                        Some(selected_topology[0]),
                        IntegrandSample::Nested(sub_sample.as_ref().unwrap().as_ref()),
                    ),
                    havana::Sample::MultiChannel(_, _, _) => unimplemented!(),
                },
            }
        } else {
            (None, sample)
        };

        // clear the deformation cache for non-stability topologies, unless
        // we are integrating an amplitude, since the deformation is fixed
        cache.current_deformation_index = 0;
        if !self.is_stability_check_topo && self.settings.cross_section.fixed_cut_momenta.is_empty()
        {
            cache.deformation_vector_cache.clear();
        }

        if self.settings.general.multi_channeling && !self.multi_channeling_channels.is_empty() {
            return self.multi_channeling(selected_topology, sub_sample, cache, event_manager);
        }

        let x: ArrayVec<[T; MAX_SG_LOOP * 3]> = sub_sample
            .to_flat()
            .iter()
            .map(|x| Into::<T>::into(*x))
            .collect();

        let mut xrot = [T::zero(); MAX_SG_LOOP * 3];

        // jointly parameterize all squared topologies
        let n_fixed = self.settings.cross_section.fixed_cut_momenta.len();
        let n_loops = x.len() / 3 + n_fixed;
        let mut k = [LorentzVector::default(); MAX_SG_LOOP];
        let mut para_jacs = [T::one(); MAX_SG_LOOP];
        let mut jac_para = T::one();
        for i in 0..n_loops {
            xrot[..x.len()].copy_from_slice(&x);

            let (l_energy, (l_space, jac)) = if i < n_loops - n_fixed {
                // set the loop index to i + 1 so that we can also shift k
                (
                    T::zero(),
                    Topology::parameterize(
                        &xrot[i * 3..(i + 1) * 3],
                        T::convert_from(&self.e_cm_squared), // NOTE: taking e_cm from the first graph
                        i,
                        &self.settings,
                    ),
                )
            } else {
                let m = self.settings.cross_section.fixed_cut_momenta[i + n_fixed - n_loops].cast();
                (
                    m.t,
                    (
                        [m.x, m.y, m.z],
                        (Into::<T>::into(2.) * <T as FloatConst>::PI()).powi(4),
                    ),
                )
            };

            // there could be some rounding here
            let rot = &self.rotation_matrix;
            k[i] = LorentzVector::from_args(
                l_energy,
                T::convert_from(&rot[0][0]) * l_space[0]
                    + T::convert_from(&rot[0][1]) * l_space[1]
                    + T::convert_from(&rot[0][2]) * l_space[2],
                T::convert_from(&rot[1][0]) * l_space[0]
                    + T::convert_from(&rot[1][1]) * l_space[1]
                    + T::convert_from(&rot[1][2]) * l_space[2],
                T::convert_from(&rot[2][0]) * l_space[0]
                    + T::convert_from(&rot[2][1]) * l_space[1]
                    + T::convert_from(&rot[2][2]) * l_space[2],
            );
            jac_para *= jac;
            para_jacs[i] = jac_para;
        }

        let mut result = Complex::zero();
        let mut event_counter = 0;
        for (current_supergraph, ((t, &m), extra_topos)) in utils::zip_eq(
            utils::zip_eq(&mut self.topologies, &self.multiplicity),
            &mut self.additional_topologies,
        )
        .enumerate()
        {
            if let Some(tt) = selected_topology {
                if current_supergraph != tt {
                    continue;
                }
            }

            cache.current_supergraph = current_supergraph;

            if self.settings.general.debug > 0 {
                println!("Evaluating supergraph {}", t.name);
            }

            // a squared topology may not use all loop variables, so we set the unused jacobians to 1
            let r = t.evaluate_mom(&k[..t.n_loops], cache, &mut event_manager)
                * para_jacs[t.n_loops - 1]
                * Into::<T>::into(m as f64);
            result += r;

            if self
                .settings
                .cross_section
                .compare_with_additional_topologies
            {
                let mut track_setting = false;
                if let Some(em) = event_manager.as_mut() {
                    track_setting = em.track_events;
                    em.track_events = false;
                }
                for (add_id, (topo, mat)) in extra_topos.iter_mut().enumerate() {
                    if self.settings.general.debug > 0 {
                        println!("Evaluating additional topology {}", add_id);
                    }

                    let mut add_ks = [LorentzVector::default(); MAX_SG_LOOP];
                    let external_momenta: ArrayVec<[LorentzVector<T>; MAX_SG_LOOP]> =
                        t.external_momenta.iter().map(|e| e.cast()).collect();

                    for (add_k, row) in add_ks.iter_mut().zip(mat) {
                        *add_k = SquaredTopology::evaluate_signature(
                            row,
                            &external_momenta,
                            &k[..t.n_loops],
                        );
                    }

                    let add_r = topo.evaluate_mom(&add_ks[..t.n_loops], cache, &mut event_manager)
                        * para_jacs[t.n_loops - 1]
                        * Into::<T>::into(m as f64);

                    if r.is_zero() || add_r.is_zero() {
                        if (r - add_r).norm_sqr()
                            > Into::<T>::into(1e-10) * T::convert_from(&self.e_cm_squared)
                        {
                            println!(
                            "Mismatch for {} between standard and additional topology {}: r={} vs add_r={}", t.name, add_id, r, add_r
                        );
                        }
                    } else if (r - add_r).norm_sqr() / (r.norm_sqr() + add_r.norm_sqr())
                        > Into::<T>::into(1e-10)
                    {
                        println!(
                        "Mismatch for {} between standard and additional topology {}: r={} vs add_r={}",
                        t.name, add_id, r, add_r
                    );
                    }
                }
                if let Some(em) = event_manager.as_mut() {
                    em.track_events = track_setting;
                }
            }

            if let Some(em) = &mut event_manager {
                if em.track_events {
                    for e in em.event_buffer[event_counter..].iter_mut() {
                        e.integrand *= para_jacs[t.n_loops - 1].to_f64().unwrap() * m as f64;
                    }
                    event_counter = em.event_buffer.len();
                }
            }
        }

        result
    }
}

pub trait DualWrapper<T> {
    fn from_real(real: T) -> Self;
    fn get_size(&self) -> usize;
    fn get_real(&self) -> T;
    fn get_der(&self, i: &[usize]) -> T;
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T;
    fn set_t(&mut self, val: T) -> &mut Self;
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>);
    fn from_zipped_slice(slice: &[T]) -> (Self, Self)
    where
        Self: Sized;
}

impl<T: Scalar + Zero + Copy> DualWrapper<T> for T {
    #[inline]
    fn from_real(real: T) -> Self {
        real
    }

    fn get_size(&self) -> usize {
        1
    }

    #[inline]
    fn get_real(&self) -> T {
        *self
    }

    #[cold]
    fn get_der(&self, _i: &[usize]) -> T {
        panic!("Bad index for dual");
    }

    #[cold]
    fn get_der_mut(&mut self, _i: &[usize]) -> &mut T {
        panic!("Bad index for dual");
    }

    #[inline]
    fn set_t(&mut self, _val: T) -> &mut Self {
        // there is no t
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[*self, *other]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (slice[0], slice[1])
    }
}

impl<T: Scalar + Zero + Copy> DualWrapper<T> for Hyperdual<T, 2> {
    #[inline]
    fn from_real(real: T) -> Self {
        Self::from_real(real)
    }

    fn get_size(&self) -> usize {
        2
    }

    #[inline]
    fn get_real(&self) -> T {
        self.real()
    }

    #[inline]
    fn get_der(&self, i: &[usize]) -> T {
        match i {
            [0] => self[1],
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T {
        match i {
            [0] => &mut self[1],
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn set_t(&mut self, val: T) -> &mut Self {
        self[1] = val;
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[self[0], other[0], self[1], other[1]]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (Self::new(slice[0], slice[2]), Self::new(slice[1], slice[3]))
    }
}

impl<T: Zero + Copy> DualWrapper<T> for Dualt2<T> {
    #[inline]
    fn from_real(real: T) -> Self {
        Self::from_real(real)
    }

    fn get_size(&self) -> usize {
        3
    }

    #[inline]
    fn get_real(&self) -> T {
        self.real
    }

    #[inline]
    fn get_der(&self, i: &[usize]) -> T {
        match i {
            [0] => self.ep_t,
            [0, 0] => self.ep_t2,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T {
        match i {
            [0] => &mut self.ep_t,
            [0, 0] => &mut self.ep_t2,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn set_t(&mut self, val: T) -> &mut Self {
        self.ep_t = val;
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[
            self.real,
            other.real,
            self.ep_t,
            other.ep_t,
            self.ep_t2,
            other.ep_t2,
        ]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (
            Dualt2::new(slice[0], slice[2], slice[4]),
            Dualt2::new(slice[1], slice[3], slice[5]),
        )
    }
}

impl<T: Zero + Copy> DualWrapper<T> for Dualkt2<T> {
    #[inline]
    fn from_real(real: T) -> Self {
        Self::from_real(real)
    }

    fn get_size(&self) -> usize {
        5
    }

    #[inline]
    fn get_real(&self) -> T {
        self.real
    }

    #[inline]
    fn get_der(&self, i: &[usize]) -> T {
        match i {
            [0] => self.ep_k0,
            [1] => self.ep_t,
            [0, 1] => self.ep_k0_t,
            [1, 1] => self.ep_t2,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T {
        match i {
            [0] => &mut self.ep_k0,
            [1] => &mut self.ep_t,
            [0, 1] => &mut self.ep_k0_t,
            [1, 1] => &mut self.ep_t2,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn set_t(&mut self, val: T) -> &mut Self {
        self.ep_t = val;
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[
            self.real,
            other.real,
            self.ep_k0,
            other.ep_k0,
            self.ep_t,
            other.ep_t,
            self.ep_k0_t,
            other.ep_k0_t,
            self.ep_t2,
            other.ep_t2,
        ]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (
            Self::new(slice[0], slice[2], slice[4], slice[6], slice[8]),
            Self::new(slice[1], slice[3], slice[5], slice[7], slice[9]),
        )
    }
}

impl<T: Zero + Copy> DualWrapper<T> for Dualt3<T> {
    #[inline]
    fn from_real(real: T) -> Self {
        Self::from_real(real)
    }

    fn get_size(&self) -> usize {
        4
    }

    #[inline]
    fn get_real(&self) -> T {
        self.real
    }

    #[inline]
    fn get_der(&self, i: &[usize]) -> T {
        match i {
            [0] => self.ep_t,
            [0, 0] => self.ep_t2,
            [0, 0, 0] => self.ep_t3,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T {
        match i {
            [0] => &mut self.ep_t,
            [0, 0] => &mut self.ep_t2,
            [0, 0, 0] => &mut self.ep_t3,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn set_t(&mut self, val: T) -> &mut Self {
        self.ep_t = val;
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[
            self.real,
            other.real,
            self.ep_t,
            other.ep_t,
            self.ep_t2,
            other.ep_t2,
            self.ep_t3,
            other.ep_t3,
        ]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (
            Self::new(slice[0], slice[2], slice[4], slice[6]),
            Self::new(slice[1], slice[3], slice[5], slice[7]),
        )
    }
}

impl<T: Zero + Copy> DualWrapper<T> for Dualkt3<T> {
    #[inline]
    fn from_real(real: T) -> Self {
        Self::from_real(real)
    }

    fn get_size(&self) -> usize {
        7
    }

    #[inline]
    fn get_real(&self) -> T {
        self.real
    }

    #[inline]
    fn get_der(&self, i: &[usize]) -> T {
        match i {
            [0] => self.ep_k0,
            [1] => self.ep_t,
            [0, 1] => self.ep_k0_t,
            [1, 1] => self.ep_t2,
            [0, 1, 1] => self.ep_k0_t2,
            [1, 1, 1] => self.ep_t3,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T {
        match i {
            [0] => &mut self.ep_k0,
            [1] => &mut self.ep_t,
            [0, 1] => &mut self.ep_k0_t,
            [1, 1] => &mut self.ep_t2,
            [0, 1, 1] => &mut self.ep_k0_t2,
            [1, 1, 1] => &mut self.ep_t3,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn set_t(&mut self, val: T) -> &mut Self {
        self.ep_t = val;
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[
            self.real,
            other.real,
            self.ep_k0,
            other.ep_k0,
            self.ep_t,
            other.ep_t,
            self.ep_k0_t,
            other.ep_k0_t,
            self.ep_t2,
            other.ep_t2,
            self.ep_k0_t2,
            other.ep_k0_t2,
            self.ep_t3,
            other.ep_t3,
        ]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (
            Self::new(
                slice[0], slice[2], slice[4], slice[6], slice[8], slice[10], slice[12],
            ),
            Self::new(
                slice[1], slice[3], slice[5], slice[7], slice[9], slice[11], slice[13],
            ),
        )
    }
}

impl<T: Zero + Copy> DualWrapper<T> for Dualklt3<T> {
    #[inline]
    fn from_real(real: T) -> Self {
        Self::from_real(real)
    }

    fn get_size(&self) -> usize {
        12
    }

    #[inline]
    fn get_real(&self) -> T {
        self.real
    }

    #[inline]
    fn get_der(&self, i: &[usize]) -> T {
        match i {
            [0] => self.ep_k0,
            [1] => self.ep_l0,
            [2] => self.ep_t,
            [0, 1] => self.ep_k0_l0,
            [0, 2] => self.ep_k0_t,
            [1, 2] => self.ep_l0_t,
            [2, 2] => self.ep_t2,
            [0, 1, 2] => self.ep_k0_l0_t,
            [0, 2, 2] => self.ep_k0_t2,
            [1, 2, 2] => self.ep_l0_t2,
            [2, 2, 2] => self.ep_t3,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn get_der_mut(&mut self, i: &[usize]) -> &mut T {
        match i {
            [0] => &mut self.ep_k0,
            [1] => &mut self.ep_l0,
            [2] => &mut self.ep_t,
            [0, 1] => &mut self.ep_k0_l0,
            [0, 2] => &mut self.ep_k0_t,
            [1, 2] => &mut self.ep_l0_t,
            [2, 2] => &mut self.ep_t2,
            [0, 1, 2] => &mut self.ep_k0_l0_t,
            [0, 2, 2] => &mut self.ep_k0_t2,
            [1, 2, 2] => &mut self.ep_l0_t2,
            [2, 2, 2] => &mut self.ep_t3,
            _ => panic!("Bad index for dual"),
        }
    }

    #[inline]
    fn set_t(&mut self, val: T) -> &mut Self {
        self.ep_t = val;
        self
    }

    #[inline]
    fn zip_flatten(&self, other: &Self, buffer: &mut Vec<T>) {
        buffer.extend(&[
            self.real,
            other.real,
            self.ep_k0,
            other.ep_k0,
            self.ep_l0,
            other.ep_l0,
            self.ep_t,
            other.ep_t,
            self.ep_k0_l0,
            other.ep_k0_l0,
            self.ep_k0_t,
            other.ep_k0_t,
            self.ep_l0_t,
            other.ep_l0_t,
            self.ep_t2,
            other.ep_t2,
            self.ep_k0_l0_t,
            other.ep_k0_l0_t,
            self.ep_k0_t2,
            other.ep_k0_t2,
            self.ep_l0_t2,
            other.ep_l0_t2,
            self.ep_t3,
            other.ep_t3,
        ]);
    }

    fn from_zipped_slice(slice: &[T]) -> (Self, Self) {
        (
            Self::new(
                slice[0], slice[2], slice[4], slice[6], slice[8], slice[10], slice[12], slice[14],
                slice[16], slice[18], slice[20], slice[22],
            ),
            Self::new(
                slice[1], slice[3], slice[5], slice[7], slice[9], slice[11], slice[13], slice[15],
                slice[17], slice[19], slice[21], slice[23],
            ),
        )
    }
}

#[derive(Default)]
pub struct SquaredTopologyCache<T: FloatLike> {
    topology_cache: Vec<Vec<Vec<Vec<LTDCache<T>>>>>,
    momenta: Vec<T>,
    params: Vec<T>,
    deformation_vector_cache: Vec<Vec<FixedDeformationLimit>>,
    deformation_jacobian_matrix: Vec<Complex<T>>,
    current_supergraph: usize,
    current_deformation_index: usize,
}

impl<T: FloatLike> SquaredTopologyCache<T> {
    #[inline]
    pub fn get_topology_cache(
        &mut self,
        cut_index: usize,
        diagram_set_index: usize,
        diagram_index: usize,
    ) -> &mut LTDCache<T> {
        &mut self.topology_cache[self.current_supergraph][cut_index][diagram_set_index]
            [diagram_index]
    }
}

/// Cache for squared topology sets.
#[derive(Default)]
pub struct SquaredTopologyCacheCollection {
    pub float_cache: SquaredTopologyCache<f64>,
    pub quad_cache: SquaredTopologyCache<f128>,
}

pub trait CachePrecisionSelector<T: FloatLike> {
    fn get(&mut self) -> &mut SquaredTopologyCache<T>;
}

impl CachePrecisionSelector<f64> for SquaredTopologyCacheCollection {
    #[inline]
    fn get(&mut self) -> &mut SquaredTopologyCache<f64> {
        &mut self.float_cache
    }
}

impl CachePrecisionSelector<f128> for SquaredTopologyCacheCollection {
    #[inline]
    fn get(&mut self) -> &mut SquaredTopologyCache<f128> {
        &mut self.quad_cache
    }
}

impl SquaredTopology {
    pub fn from_file(filename: &str, settings: &Settings) -> Result<SquaredTopology, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open squared topology file {}", filename))
            .suggestion("Does the path exist?")?;

        let mut squared_topo: SquaredTopology = serde_yaml::from_reader(f)
            .wrap_err("Could not parse squared topology file")
            .suggestion("Is it a correct yaml file")?;

        assert!(
            squared_topo.n_loops <= MAX_SG_LOOP,
            "MAX_SG_LOOP is too small: it should be at least {}",
            squared_topo.n_loops
        );

        squared_topo.settings = settings.clone();
        for cutkosky_cuts in &mut squared_topo.cutkosky_cuts {
            for cut in &mut cutkosky_cuts.cuts {
                if cut.m_squared == 0. {
                    cut.m_squared = settings.cross_section.small_mass_sq;
                }
            }

            for diagram_set in &mut cutkosky_cuts.diagram_sets {
                for d in &mut diagram_set.diagram_info {
                    d.graph.settings = settings.clone();
                    d.graph.process(false);

                    // update the UV mass
                    for ll in &mut d.graph.loop_lines {
                        for p in &mut ll.propagators {
                            if p.uv {
                                p.m_squared = settings.cross_section.m_uv_sq;
                            }
                            if p.m_squared == 0. {
                                p.m_squared = settings.cross_section.small_mass_sq;
                            }
                        }
                    }
                }
            }
        }

        // overwrite an empty list of fixed cut momenta in the topology settings if there is a default
        if settings.cross_section.fixed_cut_momenta.is_empty() {
            squared_topo.settings.cross_section.fixed_cut_momenta =
                squared_topo.default_fixed_cut_momenta.1.clone();
        }

        // set the external momenta and e_cm
        let incoming_momenta: Vec<_> = if settings.cross_section.incoming_momenta.is_empty() {
            &squared_topo.default_fixed_cut_momenta.0
        } else {
            &settings.cross_section.incoming_momenta
        }
        .iter()
        .map(|v| v.cast())
        .collect();

        let mut sum_incoming: LorentzVector<f128> = LorentzVector::default();
        for m in &incoming_momenta {
            sum_incoming += m;
        }
        squared_topo.e_cm_squared = sum_incoming.square().abs();
        squared_topo.external_momenta = incoming_momenta.clone();
        squared_topo.external_momenta.extend(incoming_momenta);

        debug_assert_eq!(
            squared_topo.external_momenta.len(),
            squared_topo.n_incoming_momenta * 2,
            "The number of external momenta is wrong."
        );

        squared_topo.rotation_matrix = [
            [f128::one(), f128::zero(), f128::zero()],
            [f128::zero(), f128::one(), f128::zero()],
            [f128::zero(), f128::zero(), f128::one()],
        ];

        let base_path = std::env::var("MG_NUMERATOR_PATH")
            .or_else(|_| {
                let mut pb = Path::new(filename)
                    .parent()
                    .ok_or("no parent")?
                    .parent()
                    .ok_or("no parent")?
                    .to_path_buf();
                let mut pb_alt = Path::new(filename)
                    .parent()
                    .ok_or("no parent")?
                    .parent()
                    .ok_or("no parent")?
                    .parent()
                    .ok_or("no parent")?
                    .to_path_buf();
                pb.push("lib");
                pb_alt.push("lib");
                if pb.exists() {
                    pb.pop();
                    Ok(pb.into_os_string().into_string().map_err(|_| "bad path")?)
                } else if pb_alt.exists() {
                    pb_alt.pop();
                    Ok(pb_alt
                        .into_os_string()
                        .into_string()
                        .map_err(|_| "bad path")?)
                } else {
                    Err("Cannot determine root folder")
                }
            })
            .expect("Cannot determine base folder from filename. Use MG_NUMERATOR_PATH");

        if let Some(cs) = &mut squared_topo.form_integrand.call_signature {
            cs.initialise(base_path);
        }

        Ok(squared_topo)
    }

    pub fn set_partial_fractioning(&mut self, enable: bool) {
        if enable {
            self.settings.cross_section.integrand_type = IntegrandType::PF;
        } else {
            self.settings.cross_section.integrand_type = IntegrandType::LTD;
        }
    }

    pub fn set_precision(&mut self, prec: usize) {
        if let Some(cs) = self.form_integrand.call_signature.as_mut() {
            cs.prec = prec;
        }
    }

    fn evaluate_signature<T: RealNumberLike + FromPrimitive, U: RealNumberLike + FromPrimitive>(
        signature: &(Vec<i8>, Vec<i8>),
        external_momenta: &[LorentzVector<T>],
        loop_momenta: &[LorentzVector<U>],
    ) -> LorentzVector<U>
    where
        U: From<T>,
    {
        let mut cut_momentum = LorentzVector::default();
        for (&sign, mom) in utils::zip_eq(&signature.0, loop_momenta) {
            if sign != 0 {
                // note: we allow for the sign to be any small integer
                cut_momentum += mom * U::from_i8(sign).unwrap();
            }
        }

        debug_assert!(signature.1.len() <= external_momenta.len());
        for (&sign, mom) in signature.1.iter().zip(external_momenta) {
            if sign != 0 {
                cut_momentum += mom.convert().multiply_sign(sign);
            }
        }
        cut_momentum
    }

    pub fn generate_multi_channeling_channels(&mut self) {
        let mut multi_channeling_channels: Vec<(
            Vec<i8>,
            Vec<i8>,
            Vec<LorentzVector<f128>>,
            Vec<f128>,
        )> = vec![];

        let multi_channeling_bases_to_consider = if self.settings.general.use_lmb_channels {
            &mut self.multi_channeling_lmb_bases
        } else {
            &mut self.multi_channeling_bases
        };

        for mcb in multi_channeling_bases_to_consider {
            if self.settings.general.use_optimal_channels {
                if let Some(selected_channel_ids) = &self.optimal_channel_ids {
                    if selected_channel_ids.len() > 0
                        && !selected_channel_ids.iter().any(|&i| i == mcb.channel_id)
                    {
                        continue;
                    }
                }
            }

            let mut cut_signatures_matrix = vec![];
            let mut lmb_to_cb_mat_i8 = vec![];
            for sig in mcb.signatures.clone() {
                for s in sig.0 {
                    cut_signatures_matrix.push(s as f64);
                    lmb_to_cb_mat_i8.push(s);
                }
            }
            let lmb_to_cb_mat = na::DMatrix::from_row_slice(
                mcb.signatures.len(),
                mcb.signatures.len(),
                &cut_signatures_matrix,
            );
            let c = lmb_to_cb_mat.try_inverse().unwrap();
            // note: this transpose is ONLY there because the iteration is column-wise instead of row-wise
            let cb_to_lmb_mat: Vec<i8> = c.transpose().iter().map(|x| *x as i8).collect();

            let mut shifts = vec![];
            for sig in mcb.signatures.clone() {
                shifts.push(utils::evaluate_signature(&sig.1, &self.external_momenta));
            }

            multi_channeling_channels.push((
                cb_to_lmb_mat.clone(),
                lmb_to_cb_mat_i8.clone(),
                shifts.clone(),
                mcb.defining_propagator_masses
                    .iter()
                    .map(|x| f128::from_f64(*x).unwrap())
                    .collect(),
            ));
        }

        self.multi_channeling_channels = multi_channeling_channels;
    }

    pub fn create_caches<T: FloatLike>(&self) -> Vec<Vec<Vec<LTDCache<T>>>> {
        let mut caches = vec![];
        for cutkosky_cuts in &self.cutkosky_cuts {
            let mut lim_cache = vec![];
            for diagram_set in &cutkosky_cuts.diagram_sets {
                let mut diag_cache = vec![];
                for d in &diagram_set.diagram_info {
                    diag_cache.push(LTDCache::<T>::new(&d.graph));
                }
                lim_cache.push(diag_cache);
            }
            caches.push(lim_cache);
        }
        caches
    }

    #[inline]
    /// Kahlen function.
    fn lambda<T: FloatLike>(s: T, sqr_mass_a: T, sqr_mass_b: T) -> T {
        s * s + sqr_mass_a * sqr_mass_a + sqr_mass_b * sqr_mass_b
            - Into::<T>::into(2.) * s * sqr_mass_a
            - Into::<T>::into(2.) * sqr_mass_b * sqr_mass_a
            - Into::<T>::into(2.) * s * sqr_mass_b
    }

    /// Solve the momentum conservation delta using Newton's method. In the case of massless propagators and an external momentum with 0 spatial part,
    /// it will take one step to find the solution. This function returns the scaling parameter and its Jacobian.
    pub fn find_scaling<T: FloatLike>(
        cutkosky_cuts: &CutkoskyCuts,
        external_momenta: &[LorentzVector<T>],
        loop_momenta: &[LorentzVector<T>],
        incoming_energy: T,
        debug_level: usize,
    ) -> Option<[(T, T); 2]> {
        let mut cut_info_cache = [(
            LorentzVector::default(),
            LorentzVector::default(),
            T::zero(),
            T::zero(),
        ); MAX_SG_LOOP + 1];

        // determine an overestimate of the t that solves the energy constraint
        // then -t and t should give two different solutions or do not converge
        let mut t_start = T::zero();
        let mut sum_k = T::zero();
        for (cut_info, cut) in cut_info_cache[..cutkosky_cuts.cuts.len()]
            .iter_mut()
            .zip(&cutkosky_cuts.cuts)
        {
            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            let shift = utils::evaluate_signature(&cut.signature.1, external_momenta);
            let k_norm_sq = k.spatial_squared();
            let k_dot_shift = k.spatial_dot(&shift);
            t_start += Float::abs(k_dot_shift) / k_norm_sq;
            sum_k += k_norm_sq.sqrt();

            *cut_info = (k, shift, k_norm_sq, k_dot_shift);
        }

        t_start += incoming_energy / sum_k;

        // find the two solutions
        let mut solutions = [(T::zero(), T::zero()); 2];
        for (i, &mut mut t) in [-t_start, t_start].iter_mut().enumerate() {
            for it in 0..20 {
                let mut f = incoming_energy;
                let mut df = T::zero();

                for (cut, cut_info) in cutkosky_cuts
                    .cuts
                    .iter()
                    .zip(&cut_info_cache[..cutkosky_cuts.cuts.len()])
                {
                    let energy = ((cut_info.0 * t + cut_info.1).spatial_squared()
                        + Into::<T>::into(cut.m_squared))
                    .sqrt();
                    f -= energy;
                    df -= (t * cut_info.2 + cut_info.3) / energy;
                }

                if debug_level > 4 {
                    println!("  | t{} finder: f={}, df={}, t={}", i, f, df, t);
                }

                if Float::abs(f) < T::epsilon() * Into::<T>::into(10.) * incoming_energy {
                    if debug_level > 2 {
                        println!("  | t{} = {}", i, t);
                    }

                    solutions[i] = (t, Float::abs(df).inv());
                    break;
                }

                if it == 19 {
                    if debug_level > 2 {
                        println!(
                            "  | no convergence after {} iterations: f={}, df={}, t={}",
                            i, f, df, t
                        );
                    }
                    return None;
                }

                t = t - f / df;
            }
        }
        if Float::abs(solutions[0].0) + Float::abs(solutions[1].0) == T::zero() {
            panic!(
                "Found exact zero solutions: {} for t={} and t={} for k={:?}, ext={:?}",
                solutions[0].0, -t_start, t_start, loop_momenta, external_momenta
            );
        }
        if Float::abs(solutions[0].0 - solutions[1].0)
            / (Float::abs(solutions[0].0) + Float::abs(solutions[1].0))
            < Into::<T>::into(1e-12)
        {
            panic!(
                "Found the same scaling solution twice: {} for t={} and t={} for k={:?}, ext={:?}",
                solutions[0].0, solutions[0].0, solutions[1].0, loop_momenta, external_momenta
            );
        }

        Some(solutions)
    }

    pub fn evaluate_mom<T: form_integrand::GetIntegrand + FloatLike>(
        &mut self,
        loop_momenta: &[LorentzVector<T>],
        cache: &mut SquaredTopologyCache<T>,
        event_manager: &mut Option<&mut EventManager>,
    ) -> Complex<T> {
        debug_assert_eq!(
            loop_momenta.len(),
            self.n_loops,
            "The number of loop momenta is wrong."
        );

        // find the index of the current event
        let event_counter = if let Some(em) = event_manager {
            em.event_buffer.len()
        } else {
            0
        };

        let external_momenta: ArrayVec<[LorentzVector<T>; MAX_SG_LOOP]> = self
            .external_momenta
            .iter()
            .map(|m| m.map(|c| T::convert_from(&c)))
            .collect();

        let mut raised_cut_powers: ArrayVec<[usize; MAX_SG_LOOP + 1]>;

        let mut result = Complex::<T>::zero();
        for cut_index in 0..self.cutkosky_cuts.len() {
            let cutkosky_cuts = &mut self.cutkosky_cuts[cut_index];
            raised_cut_powers = cutkosky_cuts
                .cuts
                .iter()
                .filter(|cc| cc.power > 1)
                .map(|cc| cc.power)
                .collect();

            if self.settings.general.debug >= 1 {
                println!(
                    "Cut {}:",
                    cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                );
            }

            let scaling_solutions = if self.settings.cross_section.do_rescaling
                && self.settings.cross_section.fixed_cut_momenta.is_empty()
            {
                let incoming_energy = external_momenta[..self.n_incoming_momenta]
                    .iter()
                    .map(|m| m.t)
                    .sum();
                SquaredTopology::find_scaling(
                    cutkosky_cuts,
                    &external_momenta[..self.external_momenta.len()],
                    loop_momenta,
                    incoming_energy,
                    self.settings.general.debug,
                )
            } else {
                // we do not scale, so give one positive solution that is just 1
                Some([(-T::one(), T::one()), (T::one(), T::one())])
            };

            if scaling_solutions.is_none() || scaling_solutions.unwrap()[1].0 < T::zero() {
                if self.settings.general.debug >= 1 {
                    println!(
                        "Phase space point has no solutions for cut {:?}:",
                        cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                    );
                }

                if let Some(em) = event_manager {
                    em.no_phase_space_counter += 1;
                }

                continue;
            }

            for &(scaling, _scaling_jac) in &scaling_solutions.unwrap() {
                if scaling < T::zero() {
                    // we ignore all negative solutions
                    continue;
                }

                macro_rules! select_eval {
                    ($( $c:pat, $dual:ty ),*) => (
                        match &raised_cut_powers[..] {
                            $($c => self.evaluate_cut::<T, $dual>(
                                &loop_momenta,
                                &external_momenta,
                                cache,
                                event_manager,
                                cut_index,
                                scaling,
                                None,
                            ),)+
                            _ => panic!("No supported dual for raised cut configuration"),
                        }
                    )
                }

                result += select_eval!(
                    [], T,
                    [2], Hyperdual<T, 2>,
                    [3], Dualt2<T>,
                    [2,2], Dualkt2<T>,
                    [4], Dualt3<T>,
                    [2,3], Dualkt3<T>,
                    [2,2,2], Dualklt3<T>);
            }
        }

        if let Some(em) = event_manager {
            if em.track_events {
                for e in em.event_buffer[event_counter..].iter_mut() {
                    e.integrand *= self.overall_numerator;
                }
            }
        }

        result *= Into::<T>::into(self.overall_numerator);

        if self.settings.general.debug >= 1 {
            println!("Final result = {:e}", result);
        }

        let test_zero_res = match self.settings.integrator.integrated_phase {
            IntegratedPhase::Real => result.re,
            IntegratedPhase::Imag => result.im,
            IntegratedPhase::Both => result.re,
        };
        if test_zero_res == T::zero() {
            if let Some(em) = event_manager {
                em.zero_eval_counter += 1;
            }
        }

        result
    }

    pub fn evaluate_cut<
        T: form_integrand::GetIntegrand + FloatLike,
        D: DualWrapper<T>
            + num::Num
            + num::Float
            + Field
            + RealNumberLike
            + From<T>
            + Signum
            + std::fmt::LowerExp
            + std::ops::Add<T, Output = D>
            + std::ops::AddAssign<D>
            + std::ops::Mul<T, Output = D>
            + std::ops::Div<T, Output = D>
            + FromPrimitive,
    >(
        &mut self,
        loop_momenta: &[LorentzVector<T>],
        external_momenta: &[LorentzVector<T>],
        cache: &mut SquaredTopologyCache<T>,
        event_manager: &mut Option<&mut EventManager>,
        cut_index: usize,
        scaling: T,
        mut deformation: Option<(&mut [LorentzVector<Complex<T>>], &mut Complex<T>)>,
    ) -> Complex<T> {
        let scaling = *D::from_real(scaling).set_t(T::one());
        let n_cuts = self.cutkosky_cuts[cut_index].cuts.len();

        let mut cut_energies_summed = D::from_real(T::zero());
        let mut scaling_result: Complex<D> = -Complex::one(); // take into account the winding number

        let mut cut_momenta: SmallVec<[LorentzVector<D>; MAX_SG_LOOP]> =
            (0..n_cuts).map(|_| LorentzVector::new()).collect();
        let mut rescaled_loop_momenta: SmallVec<[LorentzVector<D>; MAX_SG_LOOP]> =
            (0..self.n_loops).map(|_| LorentzVector::new()).collect();
        let mut subgraph_loop_momenta: SmallVec<[LorentzVector<D>; MAX_SG_LOOP]> =
            (0..self.n_loops).map(|_| LorentzVector::new()).collect();
        let mut k_def: SmallVec<[LorentzVector<Complex<D>>; MAX_SG_LOOP]> =
            (0..self.n_loops).map(|_| LorentzVector::new()).collect();
        let mut result_buffer = [T::zero(); MAX_DUAL_SIZE * 2];

        let t_prop_power = self.cutkosky_cuts[cut_index].cuts.last().unwrap().power;
        let max_extra_t_raisings: usize = self.cutkosky_cuts[cut_index].cuts[..n_cuts - 1]
            .iter()
            .map(|c| c.power - 1)
            .sum();
        let mut e_surface_expansion = D::zero();
        let first_cut_with_raising = self.cutkosky_cuts[cut_index]
            .cuts
            .iter()
            .enumerate()
            .find(|(_, c)| c.power > 1)
            .map(|(p, _)| p)
            .unwrap_or(0);

        // evaluate the cuts with the proper scaling
        for (cut_index, (cut_mom, cut)) in utils::zip_eq(
            &mut cut_momenta[..n_cuts],
            self.cutkosky_cuts[cut_index].cuts.iter(),
        )
        .enumerate()
        {
            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            let shift = utils::evaluate_signature(
                &cut.signature.1,
                &external_momenta[..self.external_momenta.len()],
            );

            *cut_mom = k.convert::<D>() * scaling + shift.convert::<D>();

            let energy = (cut_mom.spatial_squared() + Into::<T>::into(cut.m_squared)).sqrt();
            cut_energies_summed += energy;

            let k2 = k.spatial_squared();
            let kp = k.spatial_dot(&shift);
            let p2 = shift.spatial_squared();
            let m2 = Into::<T>::into(cut.m_squared);
            let e_inv = energy.get_real().inv();

            // take the n-th order expansion in t around t* of the E-surface
            e_surface_expansion += D::from_real((k2 * scaling.get_real() + kp) * e_inv);
            if t_prop_power + max_extra_t_raisings >= 2 {
                e_surface_expansion += (scaling - D::from_real(scaling.get_real()))
                    / Into::<T>::into(2.)
                    * (-kp * kp + k2 * (m2 + p2))
                    * e_inv.powi(3);
            }
            if t_prop_power + max_extra_t_raisings >= 3 {
                e_surface_expansion += (scaling - D::from_real(scaling.get_real())).powi(2)
                    / Into::<T>::into(2.)
                    * (kp * kp - k2 * (m2 + p2))
                    * (kp + k2 * scaling.get_real())
                    * e_inv.powi(5)
            }
            if t_prop_power + max_extra_t_raisings >= 4 {
                e_surface_expansion += (scaling - D::from_real(scaling.get_real())).powi(3)
                    / Into::<T>::into(8.)
                    * (kp * kp - k2 * (m2 + p2))
                    * (-Into::<T>::into(5.) * kp * kp
                        - Into::<T>::into(8.) * k2 * kp * scaling.get_real()
                        + k2 * (m2 + p2
                            - Into::<T>::into(4.) * k2 * scaling.get_real() * scaling.get_real()))
                    * e_inv.powi(7)
            }
            if t_prop_power + max_extra_t_raisings >= 5 {
                panic!(
                    "Unimplemented E-surface expansion for depth {}",
                    t_prop_power + max_extra_t_raisings
                );
            }

            if cut_index + 1 == n_cuts {
                cut_mom.t = (-cut_energies_summed + energy + cut_energies_summed.get_real())
                    .multiply_sign(cut.sign); // note: unused in practice

                e_surface_expansion = -e_surface_expansion.multiply_sign(cut.sign);

                let h_surface = (-cut_energies_summed
                    + energy * Into::<T>::into(2.0)
                    + cut_energies_summed.get_real())
                .multiply_sign(cut.sign);
                scaling_result *= num::Complex::new(
                    (h_surface * e_surface_expansion)
                        .powi(cut.power as i32)
                        .inv(),
                    D::zero(),
                );
            } else {
                let mut e1 = energy;

                if cut.power > 1 {
                    // tag the k0 derivative
                    *e1.get_der_mut(&[cut_index - first_cut_with_raising]) = T::one();
                    cut_energies_summed += e1 - energy;
                }

                cut_mom.t = e1.multiply_sign(cut.sign);

                // Construct the full H-surface, (2|k| + ep_k)^-p
                scaling_result *=
                    num::Complex::new((e1 + energy).powi(cut.power as i32).inv(), D::zero());
            }

            scaling_result *= num::Complex::new(
                D::zero(),
                D::from_real(-Into::<T>::into(2.0) * <T as FloatConst>::PI()),
            );
        }

        if self.settings.general.debug >= 1 {
            println!("  | 1/Es = {:.16e}", scaling_result);
            println!("  | q0 = {:.16e}", cut_energies_summed);
            println!("  | scaling = {:.16e}", scaling);
        }

        // h is any function that integrates to 1
        let (h, h_norm) = match self.settings.cross_section.normalising_function.name {
            NormalisingFunction::RightExponential => {
                // Only support center at 1 for now
                (
                    (-Float::powi(
                        scaling
                            / Into::<T>::into(
                                self.settings.cross_section.normalising_function.spread,
                            ),
                        2,
                    ))
                    .exp(),
                    (<T as FloatConst>::PI().sqrt() / Into::<T>::into(2.0)
                        * Into::<T>::into(self.settings.cross_section.normalising_function.spread)),
                )
            }
            NormalisingFunction::LeftRightExponential => {
                // Only support center and spread at 1 for now
                assert_eq!(self.settings.cross_section.normalising_function.center, 1.,
                        "For now the left-right exponential normalising function only support a center at 1.0.");

                let s = Into::<T>::into(self.settings.cross_section.normalising_function.spread);
                let n = (scaling.powi(2) + T::one()) * s / scaling;
                let e = (-n).exp();

                (
                    e,
                    if self.settings.cross_section.normalising_function.spread == 1. {
                        // FIXME: this is only 64 bit precise
                        /*T::from_str_radix("0.27973176363304485456919761407082", 10)
                        .ok()
                        .unwrap()*/
                        Into::<T>::into(0.27973176363304485456919761407082)
                    } else {
                        unimplemented!("No light-weight option to compute this number");
                    },
                )
            }
            NormalisingFunction::LeftRightPolynomial => {
                // Only support center and spread at 1 for now
                assert_eq!(self.settings.cross_section.normalising_function.center, 1.,
                        "For now the left-right polynomial normalising function only support a center at 1.0.");
                assert!(self.settings.cross_section.normalising_function.spread>1.,
                        "The left-right polynomial normalising function only support a spread larger than 1.0.");
                let sigma =
                    Into::<T>::into(self.settings.cross_section.normalising_function.spread);
                (
                    scaling.powf(D::from_real(sigma))
                        / (scaling.powf(D::from_real(Into::<T>::into(2.0) * sigma))
                            + Into::<T>::into(1.0)),
                    <T as FloatConst>::PI()
                        / (Into::<T>::into(2.0)
                            * sigma
                            * (<T as FloatConst>::PI() / (Into::<T>::into(2.0) * sigma)).cos()),
                )
            }
            NormalisingFunction::None => (Into::<T>::into(1.0).into(), Into::<T>::into(1.0).into()),
        };

        scaling_result *= Complex::new(
            scaling.powi(self.n_loops as i32 * 3) * h / h_norm,
            D::zero(),
        );

        // rescale the loop momenta
        for (rlm, lm) in utils::zip_eq(&mut rescaled_loop_momenta[..self.n_loops], loop_momenta) {
            *rlm = lm.convert::<D>() * scaling;
        }

        let mut constants = Complex::one();

        constants *= utils::powi(
            num::Complex::new(
                Into::<T>::into(1.)
                    / <T as Float>::powi(Into::<T>::into(2.) * <T as FloatConst>::PI(), 4),
                T::zero(),
            ),
            self.n_loops,
        );

        // multiply the flux factor
        constants /= if self.n_incoming_momenta == 2 {
            Into::<T>::into(2.)
                * (SquaredTopology::lambda(
                    (external_momenta[0] + external_momenta[1]).square().abs(),
                    external_momenta[0].square(),
                    external_momenta[1].square(),
                ))
                .sqrt()
        } else {
            let e_cm = external_momenta[0].square().sqrt().abs();
            e_cm * Into::<T>::into(2.0)
        };

        if self.settings.cross_section.picobarns {
            // return a weight in picobarns (from GeV^-2)
            constants *= Into::<T>::into(0.389379304e9);
        }

        scaling_result *= Complex::new(D::from_real(constants.re), D::from_real(constants.im));

        if self.settings.general.debug >= 2 {
            println!("  | constants={:.16e}", constants);
            println!("  | scaling part = {:.16e}", scaling_result);
            println!(
                "  | rescaled loop momenta = {:?}",
                &rescaled_loop_momenta[..self.n_loops]
            );
        }

        let e_cm_sq = if self.n_incoming_momenta == 2 {
            (external_momenta[0] + external_momenta[1]).square().abs()
        } else {
            external_momenta[0].square().abs()
        };

        if let Some(em) = event_manager {
            // TODO: performance
            let real_cut_momenta: SmallVec<[LorentzVector<T>; MAX_SG_LOOP]> = cut_momenta[..n_cuts]
                .iter()
                .map(|c| c.map(|x| x.get_real()))
                .collect();

            if !em.add_event(
                &self.external_momenta,
                &real_cut_momenta,
                &self.cutkosky_cuts[cut_index].cuts,
                &self.rotation_matrix,
            ) {
                // the event was cut
                return Complex::zero();
            }
        }

        // for the evaluation of the numerator we need complex loop momenta of the supergraph.
        // the order is: cut momenta, momenta graph 1, ... graph n
        for (kd, cut_mom) in k_def[..n_cuts - 1]
            .iter_mut()
            .zip(&cut_momenta[..n_cuts - 1])
        {
            *kd = cut_mom.to_complex(true);
        }

        // initialize the constants
        if cache.params.len() != 10 {
            cache.params = vec![
                Into::<T>::into(self.settings.cross_section.m_uv_sq).sqrt(),
                T::zero(),
                Into::<T>::into(self.settings.cross_section.mu_r_sq).sqrt(),
                T::zero(),
                Into::<T>::into(self.settings.cross_section.gs),
                T::zero(),
                Into::<T>::into(self.settings.cross_section.small_mass_sq),
                T::zero(),
                Into::<T>::into(self.settings.cross_section.uv_cutoff_scale_sq),
                T::zero(),
            ];
        }

        // now apply the same procedure for all uv limits
        let mut diag_and_num_contributions: Complex<D> = Complex::zero();
        let mut def_jacobian = Complex::one();

        // TODO: deprecate diagram sets
        let diagram_set = &mut self.cutkosky_cuts[cut_index].diagram_sets[0];

        for diagram_info in &mut diagram_set.diagram_info {
            let subgraph = &mut diagram_info.graph;

            if self.settings.general.deformation_strategy != DeformationStrategy::Fixed {
                subgraph.fixed_deformation.clear();
                continue;
            }

            // set the shifts, which are expressed in the cut basis
            for ll in &mut subgraph.loop_lines {
                for p in &mut ll.propagators {
                    p.q = SquaredTopology::evaluate_signature(
                        &p.parametric_shift,
                        &external_momenta[..self.external_momenta.len()],
                        &cut_momenta[..n_cuts],
                    )
                    .map(|x| x.to_f64().unwrap());
                }
            }

            subgraph.e_cm_squared = e_cm_sq.to_f64().unwrap();

            subgraph.update_ellipsoids();

            // the stability topologies will inherit the sources
            // amplitudes will also inherit the sources when they are computed once
            if !self.is_stability_check_topo
                && (cache.current_deformation_index >= cache.deformation_vector_cache.len()
                    || self.settings.cross_section.fixed_cut_momenta.is_empty())
            {
                subgraph.fixed_deformation = subgraph.determine_ellipsoid_overlap_structure(
                    true,
                    !self.settings.cross_section.fixed_cut_momenta.is_empty(),
                );

                cache
                    .deformation_vector_cache
                    .push(subgraph.fixed_deformation.clone());
            } else {
                subgraph.fixed_deformation =
                    cache.deformation_vector_cache[cache.current_deformation_index].clone();

                // rotate the fixed deformation vectors
                // FIXME: deformation vectors are given in f64
                let rot_matrix = &self.rotation_matrix;
                for d_lim in &mut subgraph.fixed_deformation {
                    for d in &mut d_lim.deformation_per_overlap {
                        for source in &mut d.deformation_sources {
                            let old_x = source.x;
                            let old_y = source.y;
                            let old_z = source.z;
                            source.x = rot_matrix[0][0].to_f64().unwrap() * old_x
                                + rot_matrix[0][1].to_f64().unwrap() * old_y
                                + rot_matrix[0][2].to_f64().unwrap() * old_z;
                            source.y = rot_matrix[1][0].to_f64().unwrap() * old_x
                                + rot_matrix[1][1].to_f64().unwrap() * old_y
                                + rot_matrix[1][2].to_f64().unwrap() * old_z;
                            source.z = rot_matrix[2][0].to_f64().unwrap() * old_x
                                + rot_matrix[2][1].to_f64().unwrap() * old_y
                                + rot_matrix[2][2].to_f64().unwrap() * old_z;
                        }
                    }
                }

                // update the excluded surfaces
                for ex in &mut subgraph.all_excluded_surfaces {
                    *ex = false;
                }

                for d_lim in &mut subgraph.fixed_deformation {
                    for d in &d_lim.deformation_per_overlap {
                        for surf_id in &d.excluded_surface_indices {
                            subgraph.all_excluded_surfaces[*surf_id] = true;
                        }
                    }
                }
            }
            cache.current_deformation_index += 1;

            if self.settings.general.debug > 0 {
                // check if the overlap structure makes sense
                subgraph.check_fixed_deformation();
            }
        }

        let diagram_set = &self.cutkosky_cuts[cut_index].diagram_sets[0];

        let free_loops = self.n_loops + 1 - n_cuts;
        let mat_row_length = 3 * free_loops;
        if self.settings.general.deformation_strategy == DeformationStrategy::Fixed {
            // reset jacobian matrix
            cache.deformation_jacobian_matrix.clear();
            cache
                .deformation_jacobian_matrix
                .resize(mat_row_length * mat_row_length, Complex::zero());

            def_jacobian = Complex::one();
        }

        // determine the cut momentum basis where every non-cut momentum has a dual
        let mut cut_momenta_deformation_dual: SmallVec<
            [LorentzVector<Hyperdual<T, { MAX_SG_LOOP * 3 }>>; MAX_SG_LOOP],
        > = cut_momenta[..n_cuts - 1]
            .iter()
            .map(|cm| cm.map(|x| Hyperdual::from_real(x.get_real())))
            .collect();

        let mut dual_index = 0;
        for diagram_info in diagram_set.diagram_info.iter() {
            let subgraph = &diagram_info.graph;
            for lmm in &subgraph.loop_momentum_map {
                let m = SquaredTopology::evaluate_signature(
                    lmm,
                    &external_momenta[..self.external_momenta.len()],
                    &rescaled_loop_momenta[..self.n_loops],
                );

                let mut slm = m.map(|x| Hyperdual::from_real(x.get_real()));
                slm.x[1 + dual_index * 3] = T::one();
                slm.y[1 + dual_index * 3 + 1] = T::one();
                slm.z[1 + dual_index * 3 + 2] = T::one();
                cut_momenta_deformation_dual.push(slm);
                dual_index += 1;
            }
        }

        let mut dampening_factor: Hyperdual<T, { MAX_SG_LOOP * 3 }> = Hyperdual::one();
        let mat = diagram_set.cb_to_lmb.as_ref().unwrap();

        if self.settings.general.deformation_strategy == DeformationStrategy::Fixed {
            // determine the soft dampening factor
            for p in &self.propagators {
                if p.m_squared > 0. || p.signature.0.iter().all(|s| *s == 0) {
                    continue;
                }

                let mut momentum = utils::evaluate_signature(&p.signature.1, &external_momenta)
                    .map(|x| Hyperdual::from_real(x));

                for i in 0..self.n_loops {
                    for j in 0..self.n_loops {
                        // note: transposed matrix
                        match mat[i + j * self.n_loops] * p.signature.0[j] {
                            0 => {}
                            s @ 1 | s @ -1 => {
                                let mut mom = cut_momenta_deformation_dual[i];
                                if i < n_cuts - 1 {
                                    mom -= utils::evaluate_signature(
                                        &self.cutkosky_cuts[cut_index].cuts[i].signature.1,
                                        &external_momenta[..self.external_momenta.len()],
                                    )
                                    .map(|x| Hyperdual::from_real(x));
                                }

                                if s == 1 {
                                    momentum += mom;
                                } else {
                                    momentum -= mom;
                                }
                            }
                            _ => unreachable!(),
                        }
                    }
                }

                let t = (momentum.spatial_squared() / T::convert_from(&self.e_cm_squared)).powf(
                    Hyperdual::from_real(Into::<T>::into(
                        self.settings.deformation.scaling.soft_dampening_power,
                    )),
                );

                let dampening = t
                    / (t + Into::<T>::into(
                        self.settings.deformation.scaling.soft_dampening_m.powi(2),
                    ));

                dampening_factor = dampening_factor.min(dampening);
            }

            // determine the collinear dampening factor
            for (i, c1) in self.cutkosky_cuts.iter().enumerate() {
                'next_cut: for (j, c2) in self.cutkosky_cuts.iter().enumerate() {
                    if i > j {
                        if c1.cuts.iter().any(|cc| c2.cuts.contains(cc)) {
                            let mut surface_eval = Hyperdual::zero();
                            for (diff_sign, collinear_edge) in c1
                                .cuts
                                .iter()
                                .filter_map(|cc1| {
                                    if !c2.cuts.contains(cc1) {
                                        Some((1, cc1))
                                    } else {
                                        None
                                    }
                                })
                                .chain(c2.cuts.iter().filter_map(|cc2| {
                                    if !c1.cuts.contains(cc2) {
                                        Some((-1, cc2))
                                    } else {
                                        None
                                    }
                                }))
                            {
                                if collinear_edge.m_squared > 0. {
                                    continue 'next_cut;
                                }

                                let mut momentum = LorentzVector::default();

                                for ii in 0..self.n_loops {
                                    for jj in 0..self.n_loops {
                                        // note: transposed matrix
                                        match mat[ii + jj * self.n_loops]
                                            * collinear_edge.signature.0[jj]
                                        {
                                            0 => {}
                                            s @ 1 | s @ -1 => {
                                                let mut mom = cut_momenta_deformation_dual[ii];
                                                if ii < n_cuts - 1 {
                                                    mom -= utils::evaluate_signature(
                                                        &self.cutkosky_cuts[cut_index].cuts[ii]
                                                            .signature
                                                            .1,
                                                        &external_momenta
                                                            [..self.external_momenta.len()],
                                                    )
                                                    .map(|x| Hyperdual::from_real(x));
                                                }

                                                if s == 1 {
                                                    momentum += mom;
                                                } else {
                                                    momentum -= mom;
                                                }
                                            }
                                            _ => unreachable!(),
                                        }
                                    }
                                }

                                if collinear_edge.sign * diff_sign == 1 {
                                    surface_eval += momentum.spatial_distance();
                                } else {
                                    surface_eval -= momentum.spatial_distance();
                                }
                            }

                            // TODO: normalize by the size the line segment?
                            let e = surface_eval * surface_eval
                                / (Into::<T>::into(
                                    self.settings.deformation.fixed.pinch_dampening_k_com,
                                ) * T::convert_from(&self.e_cm_squared));

                            let t = e.sqrt().powf(Hyperdual::from_real(Into::<T>::into(
                                self.settings.deformation.fixed.pinch_dampening_alpha,
                            )));

                            let dampening = t
                                / (t + Into::<T>::into(
                                    self.settings.deformation.scaling.pinch_dampening_m.powi(2),
                                ));

                            dampening_factor = dampening_factor.min(dampening);
                        }
                    }
                }
            }
        }

        // compute the deformation vectors for each subgraph
        let mut k_def_index = n_cuts - 1;
        let diagram_set = &mut self.cutkosky_cuts[cut_index].diagram_sets[0];
        for (diagram_index, diagram_info) in diagram_set.diagram_info.iter_mut().enumerate() {
            let subgraph = &diagram_info.graph;
            let subgraph_cache = cache.get_topology_cache(cut_index, 0, diagram_index);

            // do the loop momentum map, which is expressed in the loop momentum basis
            // the time component should not matter here
            for (slm, lmm) in utils::zip_eq(
                &mut subgraph_loop_momenta[..subgraph.n_loops],
                &subgraph.loop_momentum_map,
            ) {
                *slm = SquaredTopology::evaluate_signature(
                    lmm,
                    &external_momenta[..self.external_momenta.len()],
                    &rescaled_loop_momenta[..self.n_loops],
                );
            }

            let kappas = if self.settings.general.deformation_strategy == DeformationStrategy::Fixed
            {
                let s: SmallVec<[LorentzVector<T>; MAX_SG_LOOP]> = subgraph_loop_momenta
                    [..subgraph.n_loops]
                    .iter()
                    .map(|c| c.map(|x| x.get_real()))
                    .collect();

                subgraph.deform(&s, subgraph_cache, false).0
            } else {
                [LorentzVector::default(); crate::MAX_LOOP]
            };

            // add the contributions of the subgraph to the jacobian matrix, including the global dampening
            if self.settings.general.deformation_strategy == DeformationStrategy::Fixed {
                let mom_index = k_def_index + 1 - n_cuts;
                let sg_jac = std::mem::take(&mut subgraph_cache.deformation_jacobian_matrix);

                // get x and y offset
                let yoffset = mom_index * 3 * mat_row_length;
                let xoffset = mom_index * 3;

                for i in 0..3 * subgraph.n_loops {
                    // set cross terms
                    for j in 0..mat_row_length {
                        if j < xoffset || j > xoffset + 3 * subgraph.n_loops {
                            let val = Complex::new(
                                T::zero(),
                                kappas[i / 3][i % 3 + 1] * dampening_factor[1 + j],
                            );

                            if diagram_info.conjugate_deformation {
                                cache.deformation_jacobian_matrix
                                    [yoffset + i * mat_row_length + j] = -val;
                            } else {
                                cache.deformation_jacobian_matrix
                                    [yoffset + i * mat_row_length + j] = val;
                            }
                        }
                    }

                    for j in 0..3 * subgraph.n_loops {
                        let dest = &mut cache.deformation_jacobian_matrix
                            [yoffset + xoffset + i * mat_row_length + j];

                        *dest = sg_jac[i * 3 * subgraph.n_loops + j];
                        if i == j {
                            *dest -= Complex::one(); // undo the real part
                        }

                        *dest *= dampening_factor[0]; // multiply real part

                        // add derivatives of the global dampening
                        *dest += Complex::new(
                            T::zero(),
                            kappas[i / 3][i % 3 + 1] * dampening_factor[1 + mom_index * 3 + j],
                        );

                        if diagram_info.conjugate_deformation {
                            *dest = Complex::new(dest.re, -dest.im);
                        }

                        if i == j {
                            *dest += Complex::one();
                        }
                    }
                }

                cache
                    .get_topology_cache(cut_index, 0, diagram_index)
                    .deformation_jacobian_matrix = sg_jac;
            }

            for (lm, kappa) in subgraph_loop_momenta[..subgraph.n_loops]
                .iter()
                .zip(&kappas[..subgraph.n_loops])
            {
                // FIXME: this kappa needs rescaling
                k_def[k_def_index] = if diagram_info.conjugate_deformation {
                    // take the complex conjugate of the deformation
                    lm.map(|x| Complex::new(x, D::zero()))
                        - kappa.map(|x| Complex::new(D::zero(), (x * dampening_factor[0]).into()))
                } else {
                    lm.map(|x| Complex::new(x, D::zero()))
                        + kappa.map(|x| Complex::new(D::zero(), (x * dampening_factor[0]).into()))
                };

                if let Some((store_kappa, _)) = deformation.as_mut() {
                    store_kappa[k_def_index] =
                        k_def[k_def_index].map(|x| Complex::new(x.re.get_real(), x.im.get_real()));
                }

                k_def_index += 1;
            }
        }

        if self.settings.general.deformation_strategy == DeformationStrategy::Fixed {
            def_jacobian = utils::determinant(&cache.deformation_jacobian_matrix, 3 * free_loops);
        }

        if let Some((_, store_jac)) = deformation {
            *store_jac = def_jacobian;
            return Complex::zero(); // early return
        }

        cache.momenta.clear();
        for e in &external_momenta[..self.n_incoming_momenta] {
            D::from_real(e.t).zip_flatten(&D::zero(), &mut cache.momenta);
            D::from_real(e.x).zip_flatten(&D::zero(), &mut cache.momenta);
            D::from_real(e.y).zip_flatten(&D::zero(), &mut cache.momenta);
            D::from_real(e.z).zip_flatten(&D::zero(), &mut cache.momenta);
        }
        for l in &k_def[..self.n_loops] {
            l.t.re.zip_flatten(&l.t.im, &mut cache.momenta);
            l.x.re.zip_flatten(&l.x.im, &mut cache.momenta);
            l.y.re.zip_flatten(&l.y.im, &mut cache.momenta);
            l.z.re.zip_flatten(&l.z.im, &mut cache.momenta);
        }

        if let Some(em) = event_manager {
            if em.time_integrand_evaluation {
                em.integrand_evaluation_timing_start = Some(Instant::now());
            }
        }

        let mut res: Complex<D> = if let Some(call_signature) = &self.form_integrand.call_signature
        {
            let mut res = Complex::default();

            let (d, c) = (
                call_signature.id,
                if self.settings.cross_section.sum_diagram_sets {
                    1000 + cut_index
                } else {
                    diagram_set.id
                },
            );

            // FIXME: extra_calls not supported atm!
            assert!(call_signature.extra_calls.is_empty());

            for &(_diag, conf) in [(d, c)].iter().chain(&call_signature.extra_calls) {
                for r in &mut result_buffer {
                    *r = T::zero();
                }

                match self.settings.cross_section.integrand_type {
                    IntegrandType::LTD => {
                        if call_signature.prec == 0
                            || call_signature.prec == 16
                            || call_signature.prec == 32
                        {
                            T::get_integrand_ltd(
                                call_signature,
                                &cache.momenta,
                                &cache.params,
                                conf,
                                &mut result_buffer,
                            )
                        } else {
                            T::get_integrand_ltd_mpfr(
                                call_signature,
                                &cache.momenta,
                                &cache.params,
                                conf,
                                // precision is given in bits, so multiply by log(10)/log(2)
                                (call_signature.prec as f64 * 3.3219280948873624) as usize,
                                &mut result_buffer,
                            )
                        }
                    }
                    IntegrandType::PF => {
                        if call_signature.prec == 0
                            || call_signature.prec == 16
                            || call_signature.prec == 32
                        {
                            T::get_integrand_pf(
                                call_signature,
                                &cache.momenta,
                                &cache.params,
                                conf,
                                &mut result_buffer,
                            )
                        } else {
                            T::get_integrand_pf_mpfr(
                                call_signature,
                                &cache.momenta,
                                &cache.params,
                                conf,
                                // precision is given in bits, so multiply by log(10)/log(2)
                                (call_signature.prec as f64 * 3.3219280948873624) as usize,
                                &mut result_buffer,
                            )
                        }
                    }
                };

                let br = D::from_zipped_slice(&result_buffer);
                res += Complex::new(br.0, br.1);
            }
            res
        } else {
            panic!("No call signature for FORM integrand, but FORM integrand mode is enabled");
        };

        if let Some(em) = event_manager {
            if let Some(s) = em.integrand_evaluation_timing_start.take() {
                em.integrand_evaluation_timing += Instant::now().duration_since(s).as_nanos();
            }
        }

        if self.settings.general.debug >= 1 {
            println!("  | res diagram set {} = {:.16e}", diagram_set.id, res);
        }

        res *= Complex::new(D::from_real(def_jacobian.re), D::from_real(def_jacobian.im));

        diag_and_num_contributions += res;

        diag_and_num_contributions *= scaling_result;

        let raised_cut_powers: ArrayVec<[usize; MAX_SG_LOOP + 1]> = self.cutkosky_cuts[cut_index]
            .cuts
            .iter()
            .filter(|cc| cc.power > 1)
            .map(|cc| cc.power)
            .collect();

        let cut_result: Complex<T> =
            match &raised_cut_powers[..] {
                [] => Complex::new(
                    diag_and_num_contributions.re.get_real(),
                    diag_and_num_contributions.im.get_real(),
                ),
                [2] => Complex::new(
                    diag_and_num_contributions.re.get_der(&[0]),
                    diag_and_num_contributions.im.get_der(&[0]),
                ),
                [3] => Complex::new(
                    diag_and_num_contributions.re.get_der(&[0, 0]),
                    diag_and_num_contributions.im.get_der(&[0, 0]),
                ),
                [2, 2] => {
                    let r01 = Complex::new(
                        diag_and_num_contributions.re.get_der(&[0, 1]),
                        diag_and_num_contributions.im.get_der(&[0, 1]),
                    );

                    // compute the extra factor of the derivative of the t-propagator E-surface in a loop momentum energy
                    let one_extra_derivative = diag_and_num_contributions
                        * (e_surface_expansion.inv().multiply_sign(
                            -self.cutkosky_cuts[cut_index].cuts.last().unwrap().sign,
                        ) * -Into::<T>::into(2.));

                    let r11 = Complex::new(
                        one_extra_derivative.re.get_der(&[1, 1]),
                        one_extra_derivative.im.get_der(&[1, 1]),
                    );

                    r01 + r11
                }
                [4] => Complex::new(
                    diag_and_num_contributions.re.get_der(&[0, 0, 0]),
                    diag_and_num_contributions.im.get_der(&[0, 0, 0]),
                ),
                [2, 3] => {
                    let r011 = Complex::new(
                        diag_and_num_contributions.re.get_der(&[0, 1, 1]),
                        diag_and_num_contributions.im.get_der(&[0, 1, 1]),
                    );

                    let one_extra_derivative = diag_and_num_contributions
                        * (e_surface_expansion.inv().multiply_sign(
                            -self.cutkosky_cuts[cut_index].cuts.last().unwrap().sign,
                        ) * -Into::<T>::into(3.));

                    let r111 = Complex::new(
                        one_extra_derivative.re.get_der(&[1, 1, 1]),
                        one_extra_derivative.im.get_der(&[1, 1, 1]),
                    );

                    r011 + r111
                }
                [2, 2, 2] => {
                    let r012 = Complex::new(
                        diag_and_num_contributions.re.get_der(&[0, 1, 2]),
                        diag_and_num_contributions.im.get_der(&[0, 1, 2]),
                    );

                    let one_extra_derivative = diag_and_num_contributions
                        * (e_surface_expansion.inv().multiply_sign(
                            -self.cutkosky_cuts[cut_index].cuts.last().unwrap().sign,
                        ) * -Into::<T>::into(2.));

                    let two_extra_derivatives = diag_and_num_contributions
                        * (e_surface_expansion.inv().powi(2) * Into::<T>::into(6.));

                    let r022 = Complex::new(
                        one_extra_derivative.re.get_der(&[0, 2, 2]),
                        one_extra_derivative.im.get_der(&[0, 2, 2]),
                    );
                    let r122 = Complex::new(
                        one_extra_derivative.re.get_der(&[1, 2, 2]),
                        one_extra_derivative.im.get_der(&[1, 2, 2]),
                    );

                    let r222 = Complex::new(
                        two_extra_derivatives.re.get_der(&[2, 2, 2]),
                        two_extra_derivatives.im.get_der(&[2, 2, 2]),
                    );

                    r012 + r022 + r122 + r222
                }
                _ => unreachable!(),
            };

        if let Some(em) = event_manager {
            // set the event result
            if em.track_events {
                em.event_buffer.last_mut().unwrap().integrand = Complex::new(
                    cut_result.re.to_f64().unwrap(),
                    cut_result.im.to_f64().unwrap(),
                );
            }
        }

        if self.settings.general.debug >= 1 {
            println!("  | scaling res = {:e}", cut_result);
        }

        cut_result
    }

    /// Create a rotated version of this squared topology. The axis needs to be normalized.
    fn rotate(&self, angle: f128, axis: (f128, f128, f128)) -> SquaredTopology {
        let (sin_t, cos_t) = angle.sin_cos();
        let cos_t_bar = f128::one() - cos_t;

        let rot_matrix: [[f128; 3]; 3] = [
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
        rotated_topology.name += "_rot";
        rotated_topology.rotation_matrix = rot_matrix.clone();

        // remove the SOCP problem allocations from the rotated topology as we will always inherit them
        for cc in &mut rotated_topology.cutkosky_cuts {
            for ds in &mut cc.diagram_sets {
                for di in &mut ds.diagram_info {
                    di.graph.socp_problem = SOCPProblem::default();
                }
            }
        }

        for e in &mut rotated_topology.external_momenta {
            let old_x = e.x;
            let old_y = e.y;
            let old_z = e.z;
            e.x = rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z;
            e.y = rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z;
            e.z = rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z;
        }

        rotated_topology
    }
}

impl IntegrandImplementation for SquaredTopologySet {
    type Cache = SquaredTopologyCacheCollection;

    fn create_stability_check(&self, num_checks: usize) -> Vec<SquaredTopologySet> {
        let mut stability_topologies = vec![];

        let mut rng = rand::thread_rng();
        for i in 0..num_checks {
            if self.stability_check_topologies.iter().all(|x| i < x.len()) {
                let mut sts = SquaredTopologySet {
                    name: self.name.clone(),
                    e_cm_squared: self.e_cm_squared,
                    topologies: self
                        .stability_check_topologies
                        .iter()
                        .map(|sct| sct[i].clone())
                        .collect(),
                    additional_topologies: vec![vec![]; self.topologies.len()],
                    multiplicity: self.multiplicity.clone(),
                    settings: self.settings.clone(),
                    rotation_matrix: self.rotation_matrix.clone(),
                    multi_channeling_channels: vec![],
                    stability_check_topologies: vec![vec![]; self.topologies.len()],
                    is_stability_check_topo: true,
                };
                sts.create_multi_channeling_channels();
                stability_topologies.push(sts);
            } else {
                // we don't have enough stability topologies, so we pad with rotations
                let angle =
                    f128::from_f64(rng.gen::<f64>() * 2.).unwrap() * <f128 as FloatConst>::PI();
                let mut rv = (
                    f128::from_f64(rng.gen()).unwrap(),
                    f128::from_f64(rng.gen()).unwrap(),
                    f128::from_f64(rng.gen()).unwrap(),
                ); // rotation axis
                let inv_norm = (rv.0 * rv.0 + rv.1 * rv.1 + rv.2 * rv.2).sqrt().inv();
                rv = (rv.0 * inv_norm, rv.1 * inv_norm, rv.2 * inv_norm);

                stability_topologies.push(self.rotate(angle, rv))
            }
        }
        stability_topologies
    }

    fn get_target(&self) -> Option<Complex<f64>> {
        let mut target = Complex::zero();

        for (t, m) in self.topologies.iter().zip(&self.multiplicity) {
            if t.analytical_result_real.is_some() || t.analytical_result_imag.is_some() {
                target += Complex::new(
                    t.analytical_result_real.unwrap_or(0.),
                    t.analytical_result_imag.unwrap_or(0.),
                ) * *m;
            }
        }

        if !target.is_zero() {
            Some(target)
        } else {
            None
        }
    }

    fn set_partial_fractioning(&mut self, enable: bool) {
        for t in self.topologies.iter_mut() {
            t.set_partial_fractioning(enable);
        }
    }

    fn create_grid(&self) -> Grid {
        if self.topologies.len() > 1 {
            Grid::DiscreteGrid(DiscreteGrid::new(
                &[self.topologies.len()],
                self.topologies
                    .iter()
                    .map(|t| {
                        if self.settings.general.multi_channeling
                            && !t.multi_channeling_channels.is_empty()
                            && self.settings.general.multi_channeling_channel.is_none()
                        {
                            // construct the discrete grid for multi-channeling
                            // TODO: now we can replace the loop counts by the actual loop counts needed per topology
                            Grid::DiscreteGrid(DiscreteGrid::new(
                                &[t.multi_channeling_channels.len()],
                                (0..t.multi_channeling_channels.len())
                                    .map(|_| {
                                        Grid::ContinuousGrid(ContinuousGrid::new(
                                            3 * t.n_loops,
                                            self.settings.integrator.n_bins,
                                            self.settings.integrator.min_samples_for_update,
                                        ))
                                    })
                                    .collect(),
                                self.settings.integrator.max_discrete_bin_probability_ratio,
                            ))
                        } else {
                            Grid::ContinuousGrid(ContinuousGrid::new(
                                3 * self.get_maximum_loop_count(), // t.n_loops,
                                self.settings.integrator.n_bins,
                                self.settings.integrator.min_samples_for_update,
                            ))
                        }
                    })
                    .collect(),
                self.settings.integrator.max_discrete_bin_probability_ratio,
            ))
        } else {
            if self.settings.general.multi_channeling
                && !self.multi_channeling_channels.is_empty()
                && self.settings.general.multi_channeling_channel.is_none()
            {
                // construct the discrete grid for multi-channeling
                // TODO: now we can replace the loop counts by the actual loop counts needed per topology
                Grid::DiscreteGrid(DiscreteGrid::new(
                    &[self.multi_channeling_channels.len()],
                    (0..self.multi_channeling_channels.len())
                        .map(|_| {
                            Grid::ContinuousGrid(ContinuousGrid::new(
                                3 * self.get_maximum_loop_count(),
                                self.settings.integrator.n_bins,
                                self.settings.integrator.min_samples_for_update,
                            ))
                        })
                        .collect(),
                    self.settings.integrator.max_discrete_bin_probability_ratio,
                ))
            } else {
                Grid::ContinuousGrid(ContinuousGrid::new(
                    3 * self.get_maximum_loop_count(),
                    self.settings.integrator.n_bins,
                    self.settings.integrator.min_samples_for_update,
                ))
            }
        }
    }

    #[inline]
    fn evaluate_f64<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCacheCollection,
        event_manager: Option<&mut EventManager>,
    ) -> Complex<f64> {
        self.evaluate(x, cache.get(), event_manager)
    }

    #[inline]
    fn evaluate_f128<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCacheCollection,
        event_manager: Option<&mut EventManager>,
    ) -> Complex<f128> {
        self.evaluate(x, cache.get(), event_manager)
    }

    #[inline]
    fn create_cache(&self) -> SquaredTopologyCacheCollection {
        SquaredTopologyCacheCollection {
            float_cache: SquaredTopologyCache {
                topology_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
                deformation_vector_cache: vec![],
                deformation_jacobian_matrix: vec![],
                momenta: vec![],
                params: vec![],
                current_supergraph: 0,
                current_deformation_index: 0,
            },
            quad_cache: SquaredTopologyCache {
                topology_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
                deformation_vector_cache: vec![],
                deformation_jacobian_matrix: vec![],
                momenta: vec![],
                params: vec![],
                current_supergraph: 0,
                current_deformation_index: 0,
            },
        }
    }

    #[inline]
    fn set_precision(&mut self, prec: usize) {
        for t in &mut self.topologies {
            t.set_precision(prec);
        }
    }
}
