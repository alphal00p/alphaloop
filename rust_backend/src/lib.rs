#[macro_use]
extern crate itertools;
#[cfg(feature = "python_api")]
use dualklt3::Dualklt3;
#[cfg(feature = "python_api")]
use dualkt2::Dualkt2;
#[cfg(feature = "python_api")]
use dualkt3::Dualkt3;
#[cfg(feature = "python_api")]
use dualt2::Dualt2;
#[cfg(feature = "python_api")]
use dualt3::Dualt3;
#[cfg(feature = "python_api")]
use hyperdual::Hyperdual;
use num::ToPrimitive;
#[cfg(feature = "python_api")]
use pyo3::prelude::*;
#[cfg(feature = "python_api")]
use smallvec::smallvec;
extern crate nalgebra as na;

use color_eyre::{Help, Report};
use eyre::WrapErr;

use lorentz_vector::{Field, RealNumberLike};
use num_traits::{Float, FloatConst, FromPrimitive, Num, Signed};
use utils::Signum;

#[cfg(not(feature = "higher_loops"))]
pub const MAX_LOOP: usize = 4;
#[cfg(feature = "higher_loops")]
pub const MAX_LOOP: usize = 6;

pub trait FloatConvertFrom<U> {
    fn convert_from(x: &U) -> Self;
}

impl FloatConvertFrom<f64> for f64 {
    fn convert_from(x: &f64) -> f64 {
        *x
    }
}

impl FloatConvertFrom<f128::f128> for f64 {
    fn convert_from(x: &f128::f128) -> f64 {
        (*x).to_f64().unwrap()
    }
}

impl FloatConvertFrom<f128::f128> for f128::f128 {
    fn convert_from(x: &f128::f128) -> f128::f128 {
        *x
    }
}

impl FloatConvertFrom<f64> for f128::f128 {
    fn convert_from(x: &f64) -> f128::f128 {
        f128::f128::from_f64(*x).unwrap()
    }
}

pub trait FloatLike:
    From<f64>
    + FloatConvertFrom<f64>
    + FloatConvertFrom<f128::f128>
    + Num
    + FromPrimitive
    + Float
    + Field
    + RealNumberLike
    + Signed
    + FloatConst
    + std::fmt::LowerExp
    + 'static
    + Signum
{
}

impl FloatLike for f64 {}
impl FloatLike for f128::f128 {}

pub mod dashboard;
pub mod dualklt3;
pub mod dualkt2;
pub mod dualkt3;
pub mod dualt2;
pub mod dualt3;
pub mod integrand;
pub mod ltd;
pub mod observables;
pub mod squared_topologies;
pub mod topologies;
pub mod utils;

#[cfg(feature = "python_api")]
use arrayvec::ArrayVec;
pub use lorentz_vector::LorentzVector;
#[cfg(feature = "python_api")]
use num::Complex;
use serde::Deserialize;
use std::fmt;
use std::fs::File;

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
    pub soft_dampening_power: f64,
    pub theta_r_out: f64,
    pub theta_r_in: f64,
    pub theta_c: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum DeformationStrategy {
    #[serde(rename = "constant")]
    Constant,
    #[serde(rename = "fixed")]
    Fixed,
    #[serde(rename = "none")]
    None,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum IRHandling {
    #[serde(rename = "dismiss_point")]
    DismissPoint,
    #[serde(rename = "dismiss_deformation")]
    DismissDeformation,
    #[serde(rename = "none")]
    None,
}

impl From<&str> for IRHandling {
    fn from(s: &str) -> Self {
        match s {
            "dismiss_point" => IRHandling::DismissPoint,
            "dismiss_deformation" => IRHandling::DismissDeformation,
            "none" => IRHandling::None,
            _ => panic!("Unknown IR handling strategy {}", s),
        }
    }
}

impl fmt::Display for IRHandling {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            IRHandling::DismissPoint => write!(f, "dismiss_point"),
            IRHandling::DismissDeformation => write!(f, "dismiss_deformation"),
            IRHandling::None => write!(f, "none"),
        }
    }
}

impl Default for IRHandling {
    fn default() -> IRHandling {
        IRHandling::None
    }
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
    #[serde(rename = "left_right_polynomial")]
    LeftRightPolynomial,
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
            "left_right_polynomial" => NormalisingFunction::LeftRightPolynomial,
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
            NormalisingFunction::LeftRightPolynomial => {
                write!(f, "polynomial dampening to the left and right")
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
            DeformationStrategy::Constant => write!(f, "constant"),
            DeformationStrategy::Fixed => write!(f, "fixed"),
            DeformationStrategy::None => write!(f, "none"),
        }
    }
}

impl Default for DeformationStrategy {
    fn default() -> DeformationStrategy {
        DeformationStrategy::None
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
    pub pinch_dampening_alpha: f64,
    pub pinch_dampening_k_com: f64,
    pub pinch_dampening_k_shift: f64,
    pub ir_handling_strategy: IRHandling,
    pub ir_alpha: f64,
    pub ir_k_com: f64,
    pub ir_k_shift: f64,
    pub ir_interpolation_length: f64,
    pub ir_threshold: f64,
    pub ir_beta_energy: f64,
    pub ir_beta_ellipse: f64,
    pub ir_beta_pinch: f64,
    pub normalize_per_source: bool,
    pub normalisation_of_subspace_components: bool,
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

#[derive(Debug, Clone, Default, Deserialize)]
#[allow(non_snake_case)]
pub struct AFBSettings {
    pub x_min: f64,
    pub x_max: f64,
    pub n_bins: usize,
    pub write_to_file: bool,
    pub filename: String,
}

#[derive(Debug, Clone, Deserialize)]
#[allow(non_snake_case)]
pub struct RangeFilterSettings {
    pub pdgs: Vec<isize>,
    pub filter: FilterQuantity,
    pub min_value: f64,
    pub max_value: f64,
}

#[derive(Debug, Copy, Clone, Deserialize, PartialEq)]
pub enum FilterQuantity {
    #[serde(rename = "E")]
    Energy,
    #[serde(rename = "CosThetaP")]
    CosThetaP,
    #[serde(rename = "pT")]
    PT,
}

#[derive(Debug, Clone, Deserialize)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum PhaseSpaceSelectorSettings {
    #[serde(rename = "jet")]
    Jet(JetSliceSettings),
    #[serde(rename = "ranged")]
    RangeFilter(RangeFilterSettings),
}

#[derive(Debug, Clone, Deserialize)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum ObservableSettings {
    #[serde(rename = "jet1pt")]
    Jet1PT(Jet1PTSettings),
    AFB(AFBSettings),
    #[serde(rename = "cross_section")]
    CrossSection,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct StabilityCheckSettings {
    pub n_samples: usize,
    pub prec: usize,
    pub use_pf: bool,
    pub relative_precision: f64,
    pub escalate_for_large_weight_threshold: f64,
    pub minimal_precision_to_skip_further_checks: f64,
    pub accepted_radius_range_in_x_space: (f64, f64),
}

fn default_as_false() -> bool {
    false
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct GeneralSettings {
    pub multi_channeling: bool,
    pub multi_channeling_channel: Option<isize>,
    pub multi_channeling_alpha: f64,
    pub use_optimal_channels: bool,
    #[serde(default = "default_as_false")]
    pub use_lmb_channels: bool,
    pub res_file_prefix: String,
    pub derive_overlap_structure: bool,
    pub deformation_strategy: DeformationStrategy,
    pub mu_uv_sq_re_im: Vec<f64>,
    pub topology: String,
    pub numerical_threshold: f64,
    pub stability_checks: Vec<StabilityCheckSettings>,
    pub absolute_precision: f64,
    pub stability_nudge_size: f64,
    pub minimal_precision_for_returning_result: f64,
    pub debug: usize,
}

#[derive(Debug, Clone, Deserialize)]
pub struct IntegratorSettings {
    pub internal_parallelization: bool,
    pub dashboard: bool,
    pub quiet_mode: bool,
    pub show_plot: bool,
    pub integrator: Integrator,
    pub n_vec: usize,
    pub n_increase: usize,
    pub n_max: usize,
    pub n_start: usize,
    pub max_discrete_bin_probability_ratio: f64,
    pub min_samples_for_update: usize,
    pub n_bins: usize,
    pub train_on_avg: bool,
    pub learning_rate: f64,
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
    #[serde(rename = "havana")]
    Havana,
}

impl Default for IntegratorSettings {
    fn default() -> IntegratorSettings {
        IntegratorSettings {
            internal_parallelization: false,
            dashboard: false,
            quiet_mode: false,
            show_plot: false,
            integrator: Integrator::Vegas,
            n_increase: 0,
            n_vec: 1,
            n_start: 10000,
            n_max: 10000000,
            max_discrete_bin_probability_ratio: 5.,
            min_samples_for_update: 1000,
            n_bins: 128,
            train_on_avg: false,
            learning_rate: 1.5,
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
pub enum IntegrandType {
    #[serde(rename = "LTD")]
    LTD,
    #[serde(rename = "PF")]
    PF,
}

impl Default for IntegrandType {
    fn default() -> Self {
        Self::PF
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct CrossSectionSettings {
    pub picobarns: bool,
    pub do_rescaling: bool,
    #[serde(rename = "NormalisingFunction")]
    pub normalising_function: NormalisingFunctionSettings,
    pub inherit_deformation_for_uv_counterterm: bool,
    pub integrand_type: IntegrandType,
    pub compare_with_additional_topologies: bool,
    pub m_uv_sq: f64,
    pub mu_r_sq: f64,
    pub gs: f64,
    pub small_mass_sq: f64,
    pub uv_cutoff_scale_sq: f64,
    pub incoming_momenta: Vec<LorentzVector<f64>>,
    pub fixed_cut_momenta: Vec<LorentzVector<f64>>,
    pub sum_diagram_sets: bool,
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
    pub observables: Vec<ObservableSettings>,
    #[serde(rename = "Selectors")]
    pub selectors: Vec<PhaseSpaceSelectorSettings>,
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

// create a Havana submodule
#[cfg(feature = "python_api")]
fn init_havana(m: &PyModule) -> PyResult<()> {
    m.add_class::<havana::bindings::HavanaWrapper>()?;
    m.add_class::<havana::bindings::SampleWrapper>()?;
    m.add_class::<havana::bindings::GridConstructor>()?;
    m.add_class::<havana::bindings::ContinuousGridConstructor>()?;
    m.add_class::<havana::bindings::DiscreteGridConstructor>()?;
    Ok(())
}

// add bindings to the generated python module
#[cfg(feature = "python_api")]
#[pymodule]
fn ltd(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PythonCrossSection>()?;

    let submod = PyModule::new(py, "havana")?;
    init_havana(submod)?;
    m.add_submodule(submod)?;

    // workaround for https://github.com/PyO3/pyo3/issues/759
    py.run(
        "\
import sys
sys.modules['ltd.havana'] = havana
    ",
        None,
        Some(m.dict()),
    )?;

    Ok(())
}

#[cfg(feature = "python_api")]
#[pyclass(name = "CrossSection")]
struct PythonCrossSection {
    squared_topology: squared_topologies::SquaredTopology,
    integrand: integrand::Integrand<squared_topologies::SquaredTopologySet>,
    caches: squared_topologies::SquaredTopologyCache<f64>,
    caches_f128: squared_topologies::SquaredTopologyCache<f128::f128>,
    _dashboard: dashboard::Dashboard,
}

#[cfg(feature = "python_api")]
#[pymethods]
impl PythonCrossSection {
    #[new]
    #[args(cross_section_set = false)]
    fn new(
        squared_topology_file: &str,
        settings_file: &str,
        cross_section_set: Option<bool>,
    ) -> PythonCrossSection {
        use integrand::IntegrandImplementation;

        let settings = Settings::from_file(settings_file).unwrap();
        let (squared_topology, squared_topology_set) = match cross_section_set {
            Some(false) => {
                let local_squared_topology = squared_topologies::SquaredTopology::from_file(
                    squared_topology_file,
                    &settings,
                )
                .unwrap();
                let local_squared_topology_clone = local_squared_topology.clone();
                (
                    local_squared_topology,
                    squared_topologies::SquaredTopologySet::from_one(local_squared_topology_clone),
                )
            }
            Some(true) => {
                let local_squared_topology_set = squared_topologies::SquaredTopologySet::from_file(
                    squared_topology_file,
                    &settings,
                )
                .unwrap();
                (
                    local_squared_topology_set.topologies[0].clone(),
                    local_squared_topology_set,
                )
            }
            _ => unreachable!(),
        };

        let squared_topologies::SquaredTopologyCacheCollection {
            float_cache,
            quad_cache,
        } = squared_topology_set.create_cache();
        let dashboard = dashboard::Dashboard::minimal_dashboard(settings.integrator.quiet_mode);
        let integrand = integrand::Integrand::new(
            squared_topology_set.get_maximum_loop_count(),
            squared_topology_set,
            settings.clone(),
            true,
            dashboard.status_update_sender.clone(),
            0,
            None,
        );

        PythonCrossSection {
            squared_topology,
            integrand,
            caches: float_cache,
            caches_f128: quad_cache,
            _dashboard: dashboard,
        }
    }

    #[args(
        havana_updater_re = "None",
        havana_updater_im = "None",
        train_on_avg = "None",
        real_phase = "None"
    )]
    fn evaluate_integrand_havana(
        &mut self,
        py: Python,
        havana_sampler: &mut havana::bindings::HavanaWrapper,
        mut havana_updater_re: Option<&mut havana::bindings::HavanaWrapper>,
        mut havana_updater_im: Option<&mut havana::bindings::HavanaWrapper>,
        real_phase: Option<bool>,
    ) -> PyResult<()> {
        for (i, s) in havana_sampler.samples.iter().enumerate() {
            let integrand_sample = self.integrand.evaluate(
                integrand::IntegrandSample::Nested(s),
                s.get_weight(),
                match &havana_sampler.grid {
                    havana::Grid::ContinuousGrid(cg) => cg.accumulator.cur_iter + 1,
                    havana::Grid::DiscreteGrid(dg) => dg.accumulator.cur_iter + 1,
                },
            );

            let phase = if let Some(selected_real_phase) = real_phase {
                if selected_real_phase {
                    IntegratedPhase::Real
                } else {
                    IntegratedPhase::Imag
                }
            } else {
                self.integrand.settings.integrator.integrated_phase
            };

            match (&mut havana_updater_re, &mut havana_updater_im) {
                (Some(a_havana_updater_re), Some(a_havana_updater_im)) => {
                    a_havana_updater_re
                        .grid
                        .add_training_sample(s, integrand_sample.re)
                        .unwrap();
                    a_havana_updater_im
                        .grid
                        .add_training_sample(s, integrand_sample.im)
                        .unwrap();
                }
                (Some(a_havana_updater_re), _) => {
                    a_havana_updater_re
                        .grid
                        .add_training_sample(s, integrand_sample.re)
                        .unwrap();
                    if phase == IntegratedPhase::Imag {
                        havana_sampler
                            .grid
                            .add_training_sample(s, integrand_sample.im)
                            .unwrap();
                    }
                }
                (_, Some(a_havana_updater_im)) => {
                    if phase == IntegratedPhase::Real {
                        havana_sampler
                            .grid
                            .add_training_sample(s, integrand_sample.re)
                            .unwrap();
                    }
                    a_havana_updater_im
                        .grid
                        .add_training_sample(s, integrand_sample.im)
                        .unwrap();
                }
                _ => {
                    let f = match phase {
                        IntegratedPhase::Real => integrand_sample.re,
                        IntegratedPhase::Imag => integrand_sample.im,
                        IntegratedPhase::Both => unimplemented!(),
                    };
                    havana_sampler.grid.add_training_sample(s, f).unwrap();
                }
            }

            // periodically check for ctrl-c
            if i % 1000 == 0 {
                py.check_signals()?;
            }
        }

        Ok(())
    }

    #[args(sg_ids = "None", channel_ids = "None")]
    fn evaluate_integrand_batch(
        &mut self,
        py: Python,
        xs: Vec<Vec<f64>>,
        sg_ids: Option<Vec<usize>>,
        channel_ids: Option<Vec<usize>>,
    ) -> PyResult<Vec<(f64, f64)>> {
        use havana;
        use integrand::IntegrandSample;
        use itertools::izip;

        let samples = match (sg_ids, channel_ids) {
            (Some(selected_sg_ids), Some(selected_channel_ids)) => {
                let mut all_samples = vec![];
                for (sg_id, channel_id, x) in izip!(selected_sg_ids, selected_channel_ids, xs) {
                    all_samples.push(havana::Sample::DiscreteGrid(
                        1.,
                        smallvec![sg_id],
                        Some(Box::new(havana::Sample::DiscreteGrid(
                            1.,
                            smallvec![channel_id],
                            Some(Box::new(havana::Sample::ContinuousGrid(1., x))),
                        ))),
                    ))
                }
                all_samples
            }
            (Some(selected_sg_ids), _) => {
                let mut all_samples = vec![];
                for (sg_id, x) in izip!(selected_sg_ids, xs) {
                    all_samples.push(havana::Sample::DiscreteGrid(
                        1.,
                        smallvec![sg_id],
                        Some(Box::new(havana::Sample::ContinuousGrid(1., x))),
                    ))
                }
                all_samples
            }
            (_, Some(_)) => {
                panic!("Cannot sum over SGs for specific integration channels.")
            }
            (_, _) => {
                let mut all_samples = vec![];
                for x in xs {
                    all_samples.push(havana::Sample::ContinuousGrid(1., x));
                }
                all_samples
            }
        };

        let mut all_res = vec![];
        for (i, sample) in samples.iter().enumerate() {
            let res = self
                .integrand
                .evaluate(IntegrandSample::Nested(&sample), 1.0, 1);
            all_res.push((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()));

            // periodically check for ctrl-c
            if i % 1000 == 0 {
                py.check_signals()?;
            }
        }

        Ok(all_res)
    }

    fn evaluate_integrand(&mut self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self
            .integrand
            .evaluate(integrand::IntegrandSample::Flat(1., &x), 1.0, 1);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    fn evaluate(&mut self, loop_momenta: Vec<LorentzVector<f64>>) -> PyResult<(f64, f64)> {
        let res = self.squared_topology.evaluate_mom(
            &loop_momenta,
            &mut self.caches,
            &mut Some(&mut self.integrand.event_manager),
        );
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    fn get_scaling(
        &mut self,
        loop_momenta: Vec<LorentzVector<f64>>,
        cut_index: usize,
    ) -> PyResult<Option<((f64, f64), (f64, f64))>> {
        let incoming_energy: f128::f128 = self.squared_topology.external_momenta
            [..self.squared_topology.n_incoming_momenta]
            .iter()
            .map(|m| m.t)
            .sum();

        let mut moms: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut ext = Vec::with_capacity(self.squared_topology.external_momenta.len());
        for e in
            &self.squared_topology.external_momenta[..self.squared_topology.external_momenta.len()]
        {
            ext.push(e.cast());
        }

        let cutkosky_cuts = &self.squared_topology.cutkosky_cuts[cut_index];

        let scaling = squared_topologies::SquaredTopology::find_scaling(
            cutkosky_cuts,
            &ext,
            &moms[..self.squared_topology.n_loops],
            incoming_energy,
            self.squared_topology.settings.general.debug,
        );

        Ok(scaling.map(|x| {
            (
                (x[0].0.to_f64().unwrap(), x[0].1.to_f64().unwrap()),
                (x[1].0.to_f64().unwrap(), x[1].1.to_f64().unwrap()),
            )
        }))
    }

    fn evaluate_f128(&mut self, loop_momenta: Vec<LorentzVector<f64>>) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let res = self.squared_topology.evaluate_mom(
            &moms,
            &mut self.caches_f128,
            &mut Some(&mut self.integrand.event_manager),
        );
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    #[args(diagram_set = "None")]
    fn evaluate_cut(
        &mut self,
        loop_momenta: Vec<LorentzVector<f64>>,
        cut_index: usize,
        scaling: f64,
        _scaling_jac: f64,
        _diagram_set: Option<usize>,
    ) -> PyResult<(f64, f64)> {
        let (use_pf, prec) = self
            .squared_topology
            .settings
            .general
            .stability_checks
            .first()
            .map(|sc| (sc.use_pf, sc.prec))
            .unwrap();
        self.squared_topology.set_partial_fractioning(use_pf);
        self.squared_topology.set_precision(prec);

        let external_momenta: ArrayVec<[LorentzVector<f64>; MAX_LOOP]> = self
            .squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let raised_cut_powers: ArrayVec<[usize; MAX_LOOP + 4]> =
            self.squared_topology.cutkosky_cuts[cut_index]
                .cuts
                .iter()
                .filter(|cc| cc.power > 1)
                .map(|cc| cc.power)
                .collect();

        let res = match &raised_cut_powers[..] {
            [] => self.squared_topology.evaluate_cut::<f64, f64>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                None,
            ),
            [2] => self
                .squared_topology
                .evaluate_cut::<f64, Hyperdual<f64, 2>>(
                    &loop_momenta,
                    &external_momenta,
                    &mut self.caches,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    scaling,
                    None,
                ),
            [3] => self.squared_topology.evaluate_cut::<f64, Dualt2<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                None,
            ),
            [2, 2] => self.squared_topology.evaluate_cut::<f64, Dualkt2<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                None,
            ),
            [4] => self.squared_topology.evaluate_cut::<f64, Dualt3<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                None,
            ),
            [2, 3] => self.squared_topology.evaluate_cut::<f64, Dualkt3<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                None,
            ),
            [2, 2, 2] => self.squared_topology.evaluate_cut::<f64, Dualklt3<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                None,
            ),
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "No scaling could be obtained",
                ));
            }
        };

        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    #[args(diagram_set = "None")]
    fn evaluate_cut_f128(
        &mut self,
        loop_momenta: Vec<LorentzVector<f64>>,
        cut_index: usize,
        scaling: f64,
        _scaling_jac: f64,
        _diagram_set: Option<usize>,
    ) -> PyResult<(f64, f64)> {
        let (use_pf, prec) = self
            .squared_topology
            .settings
            .general
            .stability_checks
            .last()
            .map(|sc| (sc.use_pf, sc.prec))
            .unwrap();
        self.squared_topology.set_partial_fractioning(use_pf);
        self.squared_topology.set_precision(prec);

        let mut moms: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let external_momenta: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = self
            .squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let raised_cut_powers: ArrayVec<[usize; MAX_LOOP + 4]> =
            self.squared_topology.cutkosky_cuts[cut_index]
                .cuts
                .iter()
                .filter(|cc| cc.power > 1)
                .map(|cc| cc.power)
                .collect();

        let res = match &raised_cut_powers[..] {
            [] => self
                .squared_topology
                .evaluate_cut::<f128::f128, f128::f128>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            [2] => self
                .squared_topology
                .evaluate_cut::<f128::f128, Hyperdual<f128::f128, 2>>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            [3] => self
                .squared_topology
                .evaluate_cut::<f128::f128, Dualt2<f128::f128>>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            [2, 2] => self
                .squared_topology
                .evaluate_cut::<f128::f128, Dualkt2<f128::f128>>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            [4] => self
                .squared_topology
                .evaluate_cut::<f128::f128, Dualt3<f128::f128>>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            [2, 3] => self
                .squared_topology
                .evaluate_cut::<f128::f128, Dualkt3<f128::f128>>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            [2, 2, 2] => self
                .squared_topology
                .evaluate_cut::<f128::f128, Dualklt3<f128::f128>>(
                    &moms,
                    &external_momenta,
                    &mut self.caches_f128,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    f128::f128::from_f64(scaling).unwrap(),
                    None,
                ),
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "No scaling could be obtained",
                ));
            }
        };

        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    #[args(diagram_set = "None")]
    fn get_cut_deformation(
        &mut self,
        loop_momenta: Vec<LorentzVector<f64>>,
        cut_index: usize,
        scaling: f64,
        _diagram_set: Option<usize>,
    ) -> PyResult<(Vec<LorentzVector<Complex<f64>>>, (f64, f64))> {
        let external_momenta: ArrayVec<[LorentzVector<f64>; MAX_LOOP]> = self
            .squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut def_jacobian = Complex::default();
        let mut deformation = vec![LorentzVector::default(); self.squared_topology.n_loops];

        let raised_cut_powers: ArrayVec<[usize; MAX_LOOP + 4]> =
            self.squared_topology.cutkosky_cuts[cut_index]
                .cuts
                .iter()
                .filter(|cc| cc.power > 1)
                .map(|cc| cc.power)
                .collect();

        match &raised_cut_powers[..] {
            [] => self.squared_topology.evaluate_cut::<f64, f64>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                Some((&mut deformation, &mut def_jacobian)),
            ),
            [2] => self
                .squared_topology
                .evaluate_cut::<f64, Hyperdual<f64, 2>>(
                    &loop_momenta,
                    &external_momenta,
                    &mut self.caches,
                    &mut Some(&mut self.integrand.event_manager),
                    cut_index,
                    scaling,
                    Some((&mut deformation, &mut def_jacobian)),
                ),
            [3] => self.squared_topology.evaluate_cut::<f64, Dualt2<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                Some((&mut deformation, &mut def_jacobian)),
            ),
            [2, 2] => self.squared_topology.evaluate_cut::<f64, Dualkt2<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                Some((&mut deformation, &mut def_jacobian)),
            ),
            [4] => self.squared_topology.evaluate_cut::<f64, Dualt3<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                Some((&mut deformation, &mut def_jacobian)),
            ),
            [2, 3] => self.squared_topology.evaluate_cut::<f64, Dualkt3<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                Some((&mut deformation, &mut def_jacobian)),
            ),
            [2, 2, 2] => self.squared_topology.evaluate_cut::<f64, Dualklt3<f64>>(
                &loop_momenta,
                &external_momenta,
                &mut self.caches,
                &mut Some(&mut self.integrand.event_manager),
                cut_index,
                scaling,
                Some((&mut deformation, &mut def_jacobian)),
            ),
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "No scaling could be obtained",
                ));
            }
        };

        Ok((deformation, (def_jacobian.re, def_jacobian.im)))
    }

    fn parameterize(
        &mut self,
        x: Vec<f64>,
        loop_index: usize,
        e_cm_squared: f64,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::parameterize::<f64>(
            &x,
            e_cm_squared,
            loop_index,
            &self.squared_topology.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn parameterize_f128(
        &mut self,
        x: Vec<f64>,
        loop_index: usize,
        e_cm_squared: f64,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let x: Vec<_> = x
            .iter()
            .map(|v| f128::f128::from_f64(*v).unwrap())
            .collect();
        let e_cm_squared = f128::f128::from_f64(e_cm_squared).unwrap();

        let (x, jac) = topologies::Topology::parameterize::<f128::f128>(
            &x,
            e_cm_squared,
            loop_index,
            &self.squared_topology.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn inv_parameterize(
        &mut self,
        loop_momentum: LorentzVector<f64>,
        loop_index: usize,
        e_cm_squared: f64,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::inv_parametrize::<f64>(
            &loop_momentum,
            e_cm_squared,
            loop_index,
            &self.squared_topology.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn inv_parameterize_f128(
        &mut self,
        loop_momentum: LorentzVector<f64>,
        loop_index: usize,
        e_cm_squared: f64,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let e_cm_squared = f128::f128::from_f64(e_cm_squared).unwrap();

        let (x, jac) = topologies::Topology::inv_parametrize::<f128::f128>(
            &loop_momentum.cast(),
            e_cm_squared,
            loop_index,
            &self.squared_topology.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }
}
