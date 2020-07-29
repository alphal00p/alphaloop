#[cfg(feature = "python_api")]
#[macro_use]
extern crate cpython;
extern crate arrayvec;
extern crate dual_num;
extern crate mpolynomial;
#[macro_use]
extern crate itertools;
extern crate color_eyre;
extern crate eyre;
extern crate fnv;
extern crate num;
extern crate serde;
extern crate serde_yaml;
extern crate vector;

#[cfg(feature = "python_api")]
use cpython::{exc, PyErr, PyResult};
#[cfg(feature = "python_api")]
use std::cell::RefCell;
extern crate colored;
extern crate cuba;
extern crate dlopen;
extern crate f128;
extern crate libc;
extern crate nalgebra as na;
extern crate num_traits;
extern crate rand;
extern crate scs;
extern crate thousands;
extern crate tui;
#[macro_use]
extern crate dlopen_derive;

use color_eyre::{Help, Report};
use eyre::WrapErr;

use num_traits::{Float, FloatConst, FromPrimitive, Num, Signed};
#[cfg(feature = "python_api")]
use num_traits::{One, ToPrimitive, Zero};
use utils::Signum;
use vector::{Field, RealNumberLike};

pub const MAX_LOOP: usize = 4;
pub const MAX_PROP: usize = 10;

#[allow(non_camel_case_types)]
#[cfg(feature = "use_f128")]
pub type float = f128::f128;
#[allow(non_camel_case_types)]
#[cfg(not(feature = "use_f128"))]
pub type float = f64;

pub trait FloatLike:
    From<f64>
    + Num
    + FromPrimitive
    + Float
    + Field
    + RealNumberLike
    + mpolynomial::Field
    + mpolynomial::RealNumberLike
    + Signed
    + FloatConst
    + std::fmt::LowerExp
    + 'static
    + Signum
{
}

impl FloatLike for f64 {}
impl FloatLike for f128::f128 {}

pub mod amplitude;
pub mod cts;
pub mod dashboard;
pub mod gamma_chain;
pub mod integrand;
pub mod ltd;
pub mod observables;
pub mod partial_fractioning;
pub mod squared_topologies;
pub mod topologies;
pub mod utils;

#[cfg(feature = "python_api")]
use arrayvec::ArrayVec;
#[cfg(feature = "python_api")]
use num::Complex;
use serde::Deserialize;
use std::fmt;
use std::fs::File;
pub use vector::LorentzVector;

#[derive(Debug, Copy, Default, Clone, PartialEq, Deserialize)]
pub struct Scaling {
    pub lambda: f64,
    pub softmin_sigma: f64,
    pub expansion_threshold: f64,
    pub branch_cut_m: f64,
    pub branch_cut_alpha: f64,
    pub source_branch_cut_m: f64,
    pub branch_cut_check: bool,
    pub source_branch_cut_threshold: f64,
    pub source_branch_cut_multiplier: f64,
    pub expansion_check_strategy: ExpansionCheckStrategy,
    pub pole_check_strategy: PoleCheckStrategy,
    pub theta_r_out: f64,
    pub theta_r_in: f64,
    pub theta_c: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum DeformationStrategy {
    #[serde(rename = "additive")]
    Additive,
    #[serde(rename = "constant")]
    Constant,
    #[serde(rename = "fixed")]
    Fixed,
    #[serde(rename = "none")]
    None,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum ExpansionCheckStrategy {
    #[serde(rename = "first_lambda_order")]
    FirstLambdaOrder,
    #[serde(rename = "full_lambda_dependence")]
    FullLambdaDependence,
    #[serde(rename = "magic_fudge")]
    MagicFudge,
    #[serde(rename = "magic_fudge_with_min")]
    MagicFudgeWithMin,
    #[serde(rename = "ratio")]
    Ratio,
    #[serde(rename = "none")]
    None,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum NormalisingFunction {
    #[serde(rename = "right_exponential")]
    RightExponential,
    #[serde(rename = "left_right_exponential")]
    LeftRightExponential,
    #[serde(rename = "none")]
    None,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum PoleCheckStrategy {
    #[serde(rename = "exact")]
    Exact,
    #[serde(rename = "real_solution")]
    RealSolution,
    #[serde(rename = "tangent_check")]
    TangentCheck,
    #[serde(rename = "none")]
    None,
}

impl Default for PoleCheckStrategy {
    fn default() -> PoleCheckStrategy {
        PoleCheckStrategy::RealSolution
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum OverallDeformationScaling {
    #[serde(rename = "constant")]
    Constant,
    #[serde(rename = "linear")]
    Linear,
    #[serde(rename = "sigmoid")]
    Sigmoid,
    #[serde(rename = "exp_dampening")]
    ExpDampening,
}

impl Default for OverallDeformationScaling {
    fn default() -> OverallDeformationScaling {
        OverallDeformationScaling::Constant
    }
}

impl From<&str> for DeformationStrategy {
    fn from(s: &str) -> Self {
        match s {
            "additive" => DeformationStrategy::Additive,
            "constant" => DeformationStrategy::Constant,
            "fixed" => DeformationStrategy::Fixed,
            "none" => DeformationStrategy::None,
            _ => panic!("Unknown deformation strategy {}", s),
        }
    }
}

impl From<&str> for NormalisingFunction {
    fn from(s: &str) -> Self {
        match s {
            "right_exponential" => NormalisingFunction::RightExponential,
            "left_right_exponential" => NormalisingFunction::LeftRightExponential,
            "none" => NormalisingFunction::None,
            _ => panic!("Unknown normalising function {}", s),
        }
    }
}

impl fmt::Display for NormalisingFunction {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NormalisingFunction::RightExponential => {
                write!(f, "exponential dampening to the right")
            }
            NormalisingFunction::LeftRightExponential => {
                write!(f, "exponential dampening to the left and right")
            }
            NormalisingFunction::None => write!(f, "No normalising function"),
        }
    }
}

impl From<&str> for ExpansionCheckStrategy {
    fn from(s: &str) -> Self {
        match s {
            "first_lambda_order" => ExpansionCheckStrategy::FirstLambdaOrder,
            "full_lambda_dependence" => ExpansionCheckStrategy::FullLambdaDependence,
            "ratio" => ExpansionCheckStrategy::Ratio,
            "magic_fudge" => ExpansionCheckStrategy::MagicFudge,
            "magic_fudge_with_min" => ExpansionCheckStrategy::MagicFudgeWithMin,
            "none" => ExpansionCheckStrategy::None,
            _ => panic!("Unknown expansion check strategy {}", s),
        }
    }
}

impl fmt::Display for ExpansionCheckStrategy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ExpansionCheckStrategy::FirstLambdaOrder => write!(f, "first lambda order"),
            ExpansionCheckStrategy::FullLambdaDependence => write!(f, "full lambda dependence"),
            ExpansionCheckStrategy::MagicFudge => write!(f, "magic fudge"),
            ExpansionCheckStrategy::MagicFudgeWithMin => write!(f, "magic fudge with min"),
            ExpansionCheckStrategy::Ratio => write!(f, "ratio"),
            ExpansionCheckStrategy::None => write!(f, "none"),
        }
    }
}

impl From<&str> for PoleCheckStrategy {
    fn from(s: &str) -> Self {
        match s {
            "exact" => PoleCheckStrategy::Exact,
            "real_solution" => PoleCheckStrategy::RealSolution,
            "tangent_check" => PoleCheckStrategy::TangentCheck,
            "none" => PoleCheckStrategy::None,
            _ => panic!("Unknown pole check strategy {}", s),
        }
    }
}

impl fmt::Display for PoleCheckStrategy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            PoleCheckStrategy::Exact => write!(f, "exact"),
            PoleCheckStrategy::RealSolution => write!(f, "real solution"),
            PoleCheckStrategy::TangentCheck => write!(f, "tangent check"),
            PoleCheckStrategy::None => write!(f, "none"),
        }
    }
}

impl Default for ExpansionCheckStrategy {
    fn default() -> ExpansionCheckStrategy {
        ExpansionCheckStrategy::Ratio
    }
}

impl Default for NormalisingFunction {
    fn default() -> NormalisingFunction {
        NormalisingFunction::LeftRightExponential
    }
}

impl fmt::Display for DeformationStrategy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DeformationStrategy::Additive => write!(f, "additive"),
            DeformationStrategy::Constant => write!(f, "constant"),
            DeformationStrategy::Fixed => write!(f, "fixed"),
            DeformationStrategy::None => write!(f, "none"),
        }
    }
}

impl Default for DeformationStrategy {
    fn default() -> DeformationStrategy {
        DeformationStrategy::Additive
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum IntegratedPhase {
    #[serde(rename = "real")]
    Real,
    #[serde(rename = "imag")]
    Imag,
    #[serde(rename = "both")]
    Both,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum AdditiveMode {
    #[serde(rename = "exponential")]
    Exponential,
    #[serde(rename = "softmin")]
    SoftMin,
    #[serde(rename = "hyperbolic")]
    Hyperbolic,
    #[serde(rename = "unity")]
    Unity,
}

impl Default for AdditiveMode {
    fn default() -> AdditiveMode {
        AdditiveMode::Exponential
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct DeformationAdditiveSettings {
    pub mode: AdditiveMode,
    pub a_ij: f64,
    pub a_ijs: Vec<f64>,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct DeformationFixedSettings {
    pub mode: AdditiveMode,
    #[serde(rename = "M_ij")]
    pub m_ij: f64,
    pub delta: f64,
    pub sigma: f64,
    pub m_ijs: Vec<f64>,
    pub use_heuristic_centers: bool,
    pub local: bool,
    pub a_ijs: Vec<f64>,
    pub dampen_on_pinch: bool,
    pub dampen_on_pinch_after_lambda: bool,
    pub normalize_per_source: bool,
    pub normalisation_of_subspace_components: bool,
    pub normalisation_per_number_of_sources: bool,
    pub include_normal_source: bool,
    pub source_dampening_factor: f64,
    pub maximize_radius: bool,
}

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(non_snake_case)]
pub struct DeformationSettings {
    pub scaling: Scaling,
    pub normalize_on_E_surfaces_m: f64,
    pub overall_scaling: OverallDeformationScaling,
    pub overall_scaling_constant: f64,
    pub lambdas: Vec<f64>,
    pub additive: DeformationAdditiveSettings,
    pub fixed: DeformationFixedSettings,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum ParameterizationMode {
    #[serde(rename = "cartesian")]
    Cartesian,
    #[serde(rename = "spherical")]
    Spherical,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum ParameterizationMapping {
    #[serde(rename = "log")]
    Log,
    #[serde(rename = "linear")]
    Linear,
}

impl Default for ParameterizationMapping {
    fn default() -> ParameterizationMapping {
        ParameterizationMapping::Log
    }
}

impl Default for ParameterizationMode {
    fn default() -> ParameterizationMode {
        ParameterizationMode::Spherical
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct ParameterizationSettings {
    pub mode: ParameterizationMode,
    pub mapping: ParameterizationMapping,
    pub b: f64,
    pub input_rescaling: Vec<Vec<(f64, f64)>>,
    pub shifts: Vec<(f64, f64, f64, f64)>,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum ObservableMode {
    Jet1PT,
    #[serde(rename = "cross_section")]
    CrossSection,
}

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(non_snake_case)]
pub struct JetSliceSettings {
    pub min_jets: usize,
    pub max_jets: usize,
    pub min_j1pt: f64,
    pub max_j1pt: f64,
    pub dR: f64,
    pub min_jpt: f64,
    pub use_fastjet: bool,
}

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(non_snake_case)]
pub struct Jet1PTSettings {
    pub x_min: f64,
    pub x_max: f64,
    pub dR: f64,
    pub min_jpt: f64,
    pub n_bins: usize,
    pub write_to_file: bool,
    pub filename: String,
    pub use_fastjet: bool,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum SelectorMode {
    #[serde(rename = "jet")]
    Jet,
}

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(non_snake_case)]
pub struct PhaseSpaceSelectorSettings {
    pub active_selectors: Vec<SelectorMode>,
    pub jet: JetSliceSettings,
}

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(non_snake_case)]
pub struct ObservableSettings {
    pub active_observables: Vec<ObservableMode>,
    pub Jet1PT: Jet1PTSettings,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct GeneralSettings {
    pub multi_channeling: bool,
    pub multi_channeling_channel: Option<isize>,
    pub multi_channeling_including_massive_propagators: bool,
    pub log_file_prefix: String,
    pub res_file_prefix: String,
    pub screen_log_core: Option<usize>,
    pub log_points_to_screen: bool,
    pub log_stats_to_screen: bool,
    pub log_quad_upgrade: bool,
    pub derive_overlap_structure: bool,
    pub deformation_strategy: DeformationStrategy,
    pub mu_uv_sq_re_im: Vec<f64>,
    pub use_ct: bool,
    pub use_collinear_ct: bool,
    pub cut_filter: Vec<usize>,
    pub topology: String,
    pub partial_fractioning_threshold: f64,
    pub partial_fractioning_multiloop: bool,
    #[serde(default)]
    pub use_amplitude: bool,
    pub amplitude: String,
    pub unstable_point_warning_percentage: f64,
    pub numerical_threshold: f64,
    pub force_f128: bool,
    pub relative_precision_f64: f64,
    pub relative_precision_f128: f64,
    pub absolute_precision: f64,
    pub force_f128_for_large_weight_threshold: f64,
    pub stability_nudge_size: f64,
    pub numerical_instability_check: bool,
    pub minimal_precision_for_returning_result: f64,
    pub minimal_precision_to_skip_further_checks: f64,
    pub num_f64_samples: usize,
    pub num_f128_samples: usize,
    pub integration_statistics: bool,
    pub statistics_interval: usize,
    pub debug: usize,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct IntegratorSettings {
    pub internal_parallelization: bool,
    pub dashboard: bool,
    pub integrator: Integrator,
    pub n_vec: usize,
    pub n_increase: usize,
    pub n_max: usize,
    pub n_start: usize,
    pub eps_rel: f64,
    pub eps_abs: f64,
    pub border: f64,
    pub maxpass: usize,
    pub maxchisq: f64,
    pub mindeviation: f64,
    pub n_new: usize,
    pub n_min: usize,
    pub flatness: f64,
    pub seed: i32,
    pub integrated_phase: IntegratedPhase,
    pub state_filename_prefix: Option<String>,
    pub survey_n_points: usize,
    pub survey_n_iterations: usize,
    pub refine_n_runs: usize,
    pub refine_n_points: usize,
    pub load_from_state_file: bool,
    pub keep_state_file: bool,
    pub reset_vegas_integrator: bool,
    pub use_only_last_sample: bool,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum Integrator {
    #[serde(rename = "vegas")]
    Vegas,
    #[serde(rename = "suave")]
    Suave,
    #[serde(rename = "cuhre")]
    Cuhre,
    #[serde(rename = "divonne")]
    Divonne,
}

impl Default for IntegratorSettings {
    fn default() -> IntegratorSettings {
        IntegratorSettings {
            internal_parallelization: false,
            dashboard: false,
            integrator: Integrator::Vegas,
            n_increase: 0,
            n_vec: 1,
            n_start: 10000,
            n_max: 10000000,
            eps_rel: 1e-4,
            eps_abs: 0.,
            border: 1e-12,
            maxpass: 5,
            maxchisq: 10.,
            mindeviation: 0.25,
            n_new: 1000,
            n_min: 2,
            flatness: 50.,
            seed: 1,
            integrated_phase: IntegratedPhase::Real,
            state_filename_prefix: None,
            survey_n_points: 0,
            survey_n_iterations: 0,
            refine_n_runs: 0,
            refine_n_points: 0,
            load_from_state_file: false,
            keep_state_file: false,
            reset_vegas_integrator: true,
            use_only_last_sample: false,
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct NormalisingFunctionSettings {
    pub name: NormalisingFunction,
    pub spread: f64,
    pub center: f64,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum NumeratorSource {
    #[serde(rename = "yaml")]
    Yaml,
    #[serde(rename = "FORM")]
    Form,
    #[serde(rename = "FORM_integrand")]
    FormIntegrand,
}

impl Default for NumeratorSource {
    fn default() -> Self {
        Self::Yaml
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct CrossSectionSettings {
    pub picobarns: bool,
    pub do_rescaling: bool,
    #[serde(rename = "NormalisingFunction")]
    pub normalising_function: NormalisingFunctionSettings,
    pub inherit_deformation_for_uv_counterterm: bool,
    pub numerator_source: NumeratorSource,
    pub compare_with_additional_topologies: bool,
    pub m_uv_sq: f64,
    pub mu_r_sq: f64,
    pub gs: f64,
    pub small_mass_sq: f64,
    pub incoming_momenta: Vec<LorentzVector<f64>>,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct Settings {
    #[serde(rename = "General")]
    pub general: GeneralSettings,
    #[serde(rename = "Integrator")]
    pub integrator: IntegratorSettings,
    #[serde(rename = "Deformation")]
    pub deformation: DeformationSettings,
    #[serde(rename = "CrossSection")]
    pub cross_section: CrossSectionSettings,
    #[serde(rename = "Parameterization")]
    pub parameterization: ParameterizationSettings,
    #[serde(rename = "Observables")]
    pub observables: ObservableSettings,
    #[serde(rename = "Selectors")]
    pub selectors: PhaseSpaceSelectorSettings,
}

impl Settings {
    pub fn from_file(filename: &str) -> Result<Settings, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open settings file {}", filename))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .wrap_err("Could not parse settings file")
            .suggestion("Is it a correct yaml file")
    }
}

// add bindings to the generated python module
#[cfg(feature = "python_api")]
py_module_initializer!(ltd, initltd, PyInit_ltd, |py, m| {
    m.add(py, "__doc__", "LTD")?;
    m.add_class::<LTD>(py)?;
    m.add_class::<CrossSection>(py)?;
    Ok(())
});

#[cfg(feature = "python_api")]
py_class!(class CrossSection |py| {
    data squared_topology: RefCell<squared_topologies::SquaredTopology>;
    data integrand: RefCell<integrand::Integrand<squared_topologies::SquaredTopologySet>>;
    data caches: RefCell<squared_topologies::SquaredTopologyCache<f64>>;
    data caches_f128: RefCell<squared_topologies::SquaredTopologyCache<f128::f128>>;
    data dashboard: RefCell<dashboard::Dashboard>;

    def __new__(_cls, squared_topology_file: &str, settings_file: &str)
    -> PyResult<CrossSection> {
        use integrand::IntegrandImplementation;

        let settings = Settings::from_file(settings_file).unwrap();
        let squared_topology = squared_topologies::SquaredTopology::from_file(squared_topology_file, &settings).unwrap();
        let squared_topology_set = squared_topologies::SquaredTopologySet::from_one(squared_topology.clone());
        let squared_topologies::SquaredTopologyCacheCollection {float_cache, quad_cache} = squared_topology_set.create_cache();
        let dashboard = dashboard::Dashboard::minimal_dashboard();
        let integrand = integrand::Integrand::new(squared_topology.n_loops,
            squared_topology_set, settings.clone(), true, dashboard.status_update_sender.clone(), 0);

        CrossSection::create_instance(py, RefCell::new(squared_topology), RefCell::new(integrand), RefCell::new(float_cache), RefCell::new(quad_cache),
            RefCell::new(dashboard))
    }

    def evaluate_integrand(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.integrand(py).borrow_mut().evaluate(&x, 1.0, 1);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def evaluate(&self, loop_momenta: Vec<LorentzVector<f64>>) -> PyResult<(f64, f64)> {
        let mut integrand = self.integrand(py).borrow_mut();
        let res = self.squared_topology(py).borrow_mut().evaluate_mom(&loop_momenta, &mut *self.caches(py).borrow_mut(), &mut Some(&mut integrand.event_manager));
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def get_scaling(&self, loop_momenta: Vec<LorentzVector<f64>>, cut_index: usize) -> PyResult<Option<((f64, f64), (f64, f64))>> {
        let squared_topo = &mut self.squared_topology(py).borrow_mut();

        let incoming_energy = squared_topo.external_momenta[..squared_topo.n_incoming_momenta]
        .iter()
        .map(|m| m.t)
        .sum();

        let mut moms : ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut ext = Vec::with_capacity(squared_topo.external_momenta.len());
        for e in &squared_topo.external_momenta[..squared_topo.external_momenta.len()] {
            ext.push(e.cast());
        }

        let cutkosky_cuts = &squared_topo.cutkosky_cuts[cut_index];

        let scaling = squared_topologies::SquaredTopology::find_scaling(
            cutkosky_cuts,
            &ext,
            &moms[..squared_topo.n_loops],
            f128::f128::from_f64(incoming_energy).unwrap(),
            squared_topo.settings.general.debug,
        );

        Ok(scaling.map(|x| ((x[0].0.to_f64().unwrap(), x[0].1.to_f64().unwrap()), (x[1].0.to_f64().unwrap(), x[1].1.to_f64().unwrap()))))
    }

    def evaluate_f128(&self, loop_momenta: Vec<LorentzVector<f64>>) -> PyResult<(f64, f64)> {
        let mut moms : ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut integrand = self.integrand(py).borrow_mut();
        let res = self.squared_topology(py).borrow_mut().evaluate_mom(&moms, &mut *self.caches_f128(py).borrow_mut(), &mut Some(&mut integrand.event_manager));
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def evaluate_cut(&self, loop_momenta: Vec<LorentzVector<f64>>, cut_index: usize, scaling: f64, scaling_jac: f64, diagram_set: Option<usize> = None) -> PyResult<(f64, f64)> {
        let mut squared_topology = self.squared_topology(py).borrow_mut();

        let mut external_momenta: ArrayVec<[LorentzVector<float>; MAX_LOOP]> = squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];
        let cache = &mut *self.caches(py).borrow_mut();

        let mut integrand = self.integrand(py).borrow_mut();

        let n_loops = squared_topology.n_loops;
        let res = squared_topology.evaluate_cut(
            &loop_momenta,
            &mut cut_momenta,
            &mut external_momenta,
            &mut rescaled_loop_momenta,
            &mut subgraph_loop_momenta,
            &mut k_def[..n_loops],
            cache,
            &mut Some(&mut integrand.event_manager),
            cut_index,
            scaling,
            scaling_jac,
            diagram_set,
        );

        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

   def evaluate_cut_f128(&self, loop_momenta: Vec<LorentzVector<f64>>, cut_index: usize, scaling: f64, scaling_jac: f64, diagram_set: Option<usize> = None) -> PyResult<(f64, f64)> {
        let mut moms : ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut squared_topology = self.squared_topology(py).borrow_mut();

        let mut external_momenta: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];
        let cache = &mut *self.caches_f128(py).borrow_mut();

        let mut integrand = self.integrand(py).borrow_mut();

        let n_loops = squared_topology.n_loops;
        let res = squared_topology.evaluate_cut(
            &moms,
            &mut cut_momenta,
            &mut external_momenta,
            &mut rescaled_loop_momenta,
            &mut subgraph_loop_momenta,
            &mut k_def[..n_loops],
            cache,
            &mut Some(&mut integrand.event_manager),
            cut_index,
            f128::f128::from_f64(scaling).unwrap(),
            f128::f128::from_f64(scaling_jac).unwrap(),
            diagram_set,
        );

        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def get_cut_deformation(&self, loop_momenta: Vec<LorentzVector<f64>>, cut_index: usize, diagram_set: Option<usize> = None) -> PyResult<Vec<LorentzVector<Complex<f64>>>> {
        let mut moms : ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut squared_topology = self.squared_topology(py).borrow_mut();

        let mut external_momenta: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];
        let cache = &mut *self.caches_f128(py).borrow_mut();

        let mut integrand = self.integrand(py).borrow_mut();

        let incoming_energy = squared_topology.external_momenta[..squared_topology.n_incoming_momenta]
        .iter()
        .map(|m| m.t)
        .sum();

        let cutkosky_cuts = &squared_topology.cutkosky_cuts[cut_index];

        let scaling = squared_topologies::SquaredTopology::find_scaling(
            cutkosky_cuts,
            &external_momenta,
            &moms[..squared_topology.n_loops],
            f128::f128::from_f64(incoming_energy).unwrap(),
            squared_topology.settings.general.debug,
        ).ok_or_else(
            || PyErr::new::<exc::TypeError, _>(
                py,
                "No scaling could be obtained",
            )
        )?;

        let max_cuts = squared_topology.n_loops + 1;
        let num_cuts = cutkosky_cuts.cuts.len();
        squared_topology.evaluate_cut(
            &moms,
            &mut cut_momenta,
            &mut external_momenta,
            &mut rescaled_loop_momenta,
            &mut subgraph_loop_momenta,
            &mut k_def[..max_cuts],
            cache,
            &mut Some(&mut integrand.event_manager),
            cut_index,
            scaling[1].0,
            scaling[1].1,
            diagram_set,
        );

        let def = k_def[num_cuts - 1..max_cuts - 1].iter().map(|k| k.map(|d| Complex::new(d.re.into(), d.im.into()))).collect();
        Ok(def)
    }

    def parameterize(&self, x: Vec<f64>, loop_index: usize, e_cm_squared: f64) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::parameterize::<float>(&x, e_cm_squared, loop_index, &self.squared_topology(py).borrow().settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def parameterize_f128(&self, x: Vec<f64>, loop_index: usize, e_cm_squared: f64) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::parameterize::<f128::f128>(&x, e_cm_squared, loop_index, &self.squared_topology(py).borrow().settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def inv_parameterize(&self, loop_momentum: LorentzVector<f64>, loop_index: usize, e_cm_squared: f64) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::inv_parametrize::<float>(&loop_momentum, e_cm_squared, loop_index, &self.squared_topology(py).borrow().settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def inv_parameterize_f128(&self, loop_momentum: LorentzVector<f64>, loop_index: usize, e_cm_squared: f64) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::inv_parametrize::<f128::f128>(&loop_momentum.cast(), e_cm_squared, loop_index, &self.squared_topology(py).borrow().settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }
});

#[cfg(feature = "python_api")]
py_class!(class LTD |py| {
    data topo: RefCell<topologies::Topology>;
    data integrand: RefCell<integrand::Integrand<topologies::Topology>>;
    data cache: RefCell<topologies::LTDCache<float>>;
    data cache_f128: RefCell<topologies::LTDCache<f128::f128>>;
    data dashboard: RefCell<dashboard::Dashboard>;

    def __new__(_cls, topology_file: &str, top_name: &str, amplitude_file: &str, amp_name: &str, settings_file: &str)
    -> PyResult<LTD> {
        let mut settings = Settings::from_file(settings_file).unwrap();
        let mut topologies = topologies::Topology::from_file(topology_file, &settings).unwrap();
        let mut topo = topologies.remove(top_name).expect("Unknown topology");
        topo.process(true);
        //Skip amplitude with no name is set
        let mut amp= amplitude::Amplitude::default();
        if amp_name == "" {
            settings.general.use_amplitude = false;
        } else {
            let mut amplitudes = amplitude::Amplitude::from_file(amplitude_file).unwrap();
            settings.general.use_amplitude = true;
            amp=amplitudes.remove(amp_name).expect("Unknown amplitude");
            assert_eq!(amp.topology,top_name);
            amp.process(&settings.general);
        }

        topo.amplitude = amp.clone();

        let cache = topologies::LTDCache::<float>::new(&topo);
        let cache_f128 = topologies::LTDCache::<f128::f128>::new(&topo);
        let dashboard = dashboard::Dashboard::minimal_dashboard();
        let integrand = integrand::Integrand::new(topo.n_loops, topo.clone(), settings.clone(), false, dashboard.status_update_sender.clone(), 0);

        LTD::create_instance(py, RefCell::new(topo), RefCell::new(integrand), RefCell::new(cache), RefCell::new(cache_f128),
            RefCell::new(dashboard))
    }

    def __copy__(&self) -> PyResult<LTD> {
        let topo = self.topo(py).borrow();
        let settings = self.integrand(py).borrow().settings.clone();
        let dashboard = dashboard::Dashboard::minimal_dashboard();
        let integrand = integrand::Integrand::new(topo.n_loops, topo.clone(), settings, false, dashboard.status_update_sender.clone(), 0);
        let cache = topologies::LTDCache::<float>::new(&topo);
        let cache_f128 = topologies::LTDCache::<f128::f128>::new(&topo);
        LTD::create_instance(py, RefCell::new(topo.clone()), RefCell::new(integrand), RefCell::new(cache), RefCell::new(cache_f128), RefCell::new(dashboard))
    }

   def evaluate_integrand(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.integrand(py).borrow_mut().evaluate(&x, 1.0, 1);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
   }

   def evaluate(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.topo(py).borrow_mut().evaluate::<float>(&x,
            &mut self.cache(py).borrow_mut());
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

   def evaluate_f128(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.topo(py).borrow_mut().evaluate::<f128::f128>(&x,
            &mut self.cache_f128(py).borrow_mut());
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def parameterize(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let t = self.topo(py).borrow();
        let (x, jac) = topologies::Topology::parameterize::<float>(&x, t.e_cm_squared, loop_index, &t.settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def parameterize_f128(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let t = self.topo(py).borrow();
        let (x, jac) = topologies::Topology::parameterize::<f128::f128>(&x, t.e_cm_squared, loop_index, &t.settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def inv_parameterize(&self, loop_momentum: LorentzVector<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let t = self.topo(py).borrow();
        let (x, jac) = topologies::Topology::inv_parametrize::<float>(&loop_momentum, t.e_cm_squared, loop_index, &t.settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def inv_parameterize_f128(&self, loop_momentum: LorentzVector<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let t = self.topo(py).borrow();
        let (x, jac) = topologies::Topology::inv_parametrize::<f128::f128>(&loop_momentum.cast(), t.e_cm_squared, loop_index, &t.settings);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def evaluate_cut(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(float::from_f64(l[0].0).unwrap(), float::from_f64(l[0].1).unwrap()),
                Complex::new(float::from_f64(l[1].0).unwrap(), float::from_f64(l[1].1).unwrap()),
                Complex::new(float::from_f64(l[2].0).unwrap(), float::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        if topo.compute_complex_cut_energies(&moms, &mut cache).is_err() {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = topo.settings.general.partial_fractioning_threshold >= 0.0
            && moms.len() > 0
            && moms[0].real().spatial_distance()
                > topo.settings.general.partial_fractioning_threshold;

        // Prepare numerator
        topo.numerator.evaluate_reduced_in_lb(&moms, 0, &mut cache, 0, true, use_partial_fractioning);

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        match topo.evaluate_cut::<float>(&mut moms, &topo.numerator, cut, mat, &mut cache, true, 0) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.))
        }
    }

    def evaluate_cut_f128(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(f128::f128::from_f64(l[0].0).unwrap(), f128::f128::from_f64(l[0].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[1].0).unwrap(), f128::f128::from_f64(l[1].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[2].0).unwrap(), f128::f128::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        if topo.compute_complex_cut_energies(&moms, &mut cache).is_err() {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = topo.settings.general.partial_fractioning_threshold >= 0.0
        && moms.len() > 0
        && moms[0].real().spatial_distance()
            > f128::f128::from_f64(topo.settings.general.partial_fractioning_threshold).unwrap();

        // Prepare numerator
        topo.numerator.evaluate_reduced_in_lb(&moms, 0, &mut cache, 0, true, use_partial_fractioning);

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        match topo.evaluate_cut::<f128::f128>(&mut moms, &topo.numerator, cut, mat, &mut cache, true, 0) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.))
        }
    }

    def evaluate_cut_ct(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(float::from_f64(l[0].0).unwrap(), float::from_f64(l[0].1).unwrap()),
                Complex::new(float::from_f64(l[1].0).unwrap(), float::from_f64(l[1].1).unwrap()),
                Complex::new(float::from_f64(l[2].0).unwrap(), float::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        if topo.compute_complex_cut_energies(&moms, &mut cache).is_err() {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = topo.settings.general.partial_fractioning_threshold >= 0.0
        && moms.len() > 0
        && moms[0].real().spatial_distance()
            > topo.settings.general.partial_fractioning_threshold;

        // Prepare numerator
        topo.numerator.evaluate_reduced_in_lb(&moms, 0, &mut cache, 0, true, use_partial_fractioning); //NOTE: Only necessary when k_vec is changed

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];
        let v = topo.evaluate_cut::<float>(&mut moms, &topo.numerator, cut, mat, &mut cache, true, 0).unwrap();
        // get the loop line result from the cache if possible
        let r = 2.0 * cache.complex_cut_energies[cut_index];
        let ct = topo.counterterm::<float>(&moms[..topo.n_loops], r, cut_index, &mut cache);

        let res = v * (ct + float::one());
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def evaluate_cut_ct_f128(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(f128::f128::from_f64(l[0].0).unwrap(), f128::f128::from_f64(l[0].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[1].0).unwrap(), f128::f128::from_f64(l[1].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[2].0).unwrap(), f128::f128::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        if topo.compute_complex_cut_energies(&moms, &mut cache).is_err() {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = topo.settings.general.partial_fractioning_threshold >= 0.0
        && moms.len() > 0
        && moms[0].real().spatial_distance()
            > f128::f128::from_f64(topo.settings.general.partial_fractioning_threshold).unwrap();

        // Prepare numerator
        topo.numerator.evaluate_reduced_in_lb(&moms, 0, &mut cache, 0, true, use_partial_fractioning); //NOTE: Only necessary when k_vec is changed

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];
        let v = topo.evaluate_cut::<f128::f128>(&mut moms, &topo.numerator, cut, mat, &mut cache, true, 0).unwrap();
        // get the loop line result from the cache if possible
        let r = cache.complex_cut_energies[cut_index] * f128::f128::from_f64(2.0).unwrap();
        let ct = topo.counterterm::<f128::f128>(&moms[..topo.n_loops], r, cut_index, &mut cache);

        match v*(Complex::new(f128::f128::from_f64(1.0).unwrap(), f128::f128::zero())+ct){
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    }

   def evaluate_amplitude_cut(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let mut topo = self.topo(py).borrow_mut();
        topo.settings.general.use_amplitude = true;
        let mut cache = self.cache(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(0.0, 0.0),
                Complex::new(l[0].0, l[0].1),
                Complex::new(l[1].0, l[1].1),
                Complex::new(l[2].0, l[2].1)));
        }

        // FIXME: recomputed every time
        if topo.compute_complex_cut_energies(&moms, &mut cache).is_err() {
            return Ok((0., 0.));
        }

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];
        match topo.evaluate_amplitude_cut::<float>(&mut moms, cut, mat, &mut cache).unwrap() {
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    }

   def evaluate_amplitude_cut_f128(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let mut topo = self.topo(py).borrow_mut();
        topo.settings.general.use_amplitude = true;
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(f128::f128::from_f64(l[0].0).unwrap(), f128::f128::from_f64(l[0].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[1].0).unwrap(), f128::f128::from_f64(l[1].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[2].0).unwrap(), f128::f128::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        if topo.compute_complex_cut_energies(&moms, &mut cache).is_err() {
            return Ok((0., 0.));
        }

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];
        match topo.evaluate_amplitude_cut::<f128::f128>(&mut moms, cut, mat, &mut cache).unwrap() {
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    }

       def get_loop_momentum_energies(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<Vec<(f64, f64)>> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(float::from_f64(l[0].0).unwrap(), float::from_f64(l[0].1).unwrap()),
                Complex::new(float::from_f64(l[1].0).unwrap(), float::from_f64(l[1].1).unwrap()),
                Complex::new(float::from_f64(l[2].0).unwrap(), float::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        topo.compute_complex_cut_energies(&moms, &mut cache).unwrap_or_else(|_| {println!("On-shell propagator detected");});

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        topo.set_loop_momentum_energies::<float>(&mut moms, cut, mat, &cache);

        let mut res = Vec::with_capacity(moms.len());
        for l in moms {
            res.push((l.t.re.to_f64().unwrap(), l.t.im.to_f64().unwrap()));
        }

        Ok(res)
    }

    def get_loop_momentum_energies_f128(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<Vec<(f64, f64)>> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms : ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(f128::f128::from_f64(l[0].0).unwrap(), f128::f128::from_f64(l[0].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[1].0).unwrap(), f128::f128::from_f64(l[1].1).unwrap()),
                Complex::new(f128::f128::from_f64(l[2].0).unwrap(), f128::f128::from_f64(l[2].1).unwrap())));
        }

        // FIXME: recomputed every time
        topo.compute_complex_cut_energies(&moms, &mut cache).unwrap_or_else(|_| {println!("On-shell propagator detected");});

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        topo.set_loop_momentum_energies::<f128::f128>(&mut moms, cut, mat, &cache);

        let mut res = Vec::with_capacity(moms.len());
        for l in moms {
            res.push((l.t.re.to_f64().unwrap(), l.t.im.to_f64().unwrap()));
        }

        Ok(res)
    }

    def deform(&self, loop_momenta: Vec<Vec<f64>>) -> PyResult<(Vec<(f64, f64, f64)>, f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms = Vec::with_capacity(loop_momenta.len());
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(float::zero(),
                float::from_f64(l[0]).unwrap(),
                float::from_f64(l[1]).unwrap(),
                float::from_f64(l[2]).unwrap()));
        }

        let (res, jac) = topo.deform::<float>(&moms, &mut cache);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..topo.n_loops].iter() {
            r.push((x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), x[3].to_f64().unwrap()));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }

    def deform_f128(&self, loop_momenta: Vec<Vec<f64>>) -> PyResult<(Vec<(f64, f64, f64)>, f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms = Vec::with_capacity(loop_momenta.len());
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(f128::f128::zero(),
                f128::f128::from_f64(l[0]).unwrap(),
                f128::f128::from_f64(l[1]).unwrap(),
                f128::f128::from_f64(l[2]).unwrap()));
        }

        let (res, jac) = topo.deform::<f128::f128>(&moms, &mut cache);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..topo.n_loops].iter() {
            r.push((x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), x[3].to_f64().unwrap()));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }
});
