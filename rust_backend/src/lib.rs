#[macro_use]
extern crate itertools;
#[cfg(feature = "python_api")]
use pyo3::prelude::*;
#[macro_use]
extern crate dlopen_derive;
extern crate nalgebra as na;

use color_eyre::{Help, Report};
use eyre::WrapErr;

use lorentz_vector::{Field, RealNumberLike};
use num_traits::{Float, FloatConst, FromPrimitive, Num, Signed};
#[cfg(feature = "python_api")]
use num_traits::{One, ToPrimitive, Zero};
use utils::Signum;

#[cfg(not(feature = "higher_loops"))]
pub const MAX_LOOP: usize = 4;
#[cfg(feature = "higher_loops")]
pub const MAX_LOOP: usize = 6;

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
pub struct StabilityCheckSettings {
    pub n_samples: usize,
    pub use_f128: bool,
    pub use_pf: bool,
    pub relative_precision: f64,
    pub escalate_for_large_weight_threshold: f64,
    pub minimal_precision_to_skip_further_checks: f64,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct GeneralSettings {
    pub multi_channeling: bool,
    pub multi_channeling_channel: Option<isize>,
    pub multi_channeling_including_massive_propagators: bool,
    pub res_file_prefix: String,
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
    pub numerical_threshold: f64,
    pub stability_checks: Vec<StabilityCheckSettings>,
    pub absolute_precision: f64,
    pub stability_nudge_size: f64,
    pub minimal_precision_for_returning_result: f64,
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
    pub min_probability_per_bin: f64,
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
            integrator: Integrator::Vegas,
            n_increase: 0,
            n_vec: 1,
            n_start: 10000,
            n_max: 10000000,
            min_probability_per_bin: 0.1,
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
    pub numerator_source: NumeratorSource,
    pub integrand_type: IntegrandType,
    pub compare_with_additional_topologies: bool,
    pub m_uv_sq: f64,
    pub mu_r_sq: f64,
    pub gs: f64,
    pub small_mass_sq: f64,
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
#[pymodule]
fn ltd(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PythonLTD>()?;
    m.add_class::<PythonCrossSection>()?;
    Ok(())
}

#[cfg(feature = "python_api")]
#[pyclass(name = CrossSection)]
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
    fn new(squared_topology_file: &str, settings_file: &str) -> PythonCrossSection {
        use integrand::IntegrandImplementation;

        let settings = Settings::from_file(settings_file).unwrap();
        let squared_topology =
            squared_topologies::SquaredTopology::from_file(squared_topology_file, &settings)
                .unwrap();
        let squared_topology_set =
            squared_topologies::SquaredTopologySet::from_one(squared_topology.clone());
        let squared_topologies::SquaredTopologyCacheCollection {
            float_cache,
            quad_cache,
        } = squared_topology_set.create_cache();
        let dashboard = dashboard::Dashboard::minimal_dashboard();
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
        let incoming_energy = self.squared_topology.external_momenta
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
            f128::f128::from_f64(incoming_energy).unwrap(),
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
        scaling_jac: f64,
        diagram_set: Option<usize>,
    ) -> PyResult<(f64, f64)> {
        let mut external_momenta: ArrayVec<[LorentzVector<float>; MAX_LOOP]> = self
            .squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];

        let n_loops = self.squared_topology.n_loops;
        let res = self.squared_topology.evaluate_cut(
            &loop_momenta,
            &mut cut_momenta,
            &mut external_momenta,
            &mut rescaled_loop_momenta,
            &mut subgraph_loop_momenta,
            &mut k_def[..n_loops],
            &mut self.caches,
            &mut Some(&mut self.integrand.event_manager),
            cut_index,
            scaling,
            scaling_jac,
            diagram_set,
        );

        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    #[args(diagram_set = "None")]
    fn evaluate_cut_f128(
        &mut self,
        loop_momenta: Vec<LorentzVector<f64>>,
        cut_index: usize,
        scaling: f64,
        scaling_jac: f64,
        diagram_set: Option<usize>,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut external_momenta: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = self
            .squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];
        let n_loops = self.squared_topology.n_loops;
        let res = self.squared_topology.evaluate_cut(
            &moms,
            &mut cut_momenta,
            &mut external_momenta,
            &mut rescaled_loop_momenta,
            &mut subgraph_loop_momenta,
            &mut k_def[..n_loops],
            &mut self.caches_f128,
            &mut Some(&mut self.integrand.event_manager),
            cut_index,
            f128::f128::from_f64(scaling).unwrap(),
            f128::f128::from_f64(scaling_jac).unwrap(),
            diagram_set,
        );

        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    #[args(diagram_set = "None")]
    fn get_cut_deformation(
        &mut self,
        loop_momenta: Vec<LorentzVector<f64>>,
        cut_index: usize,
        diagram_set: Option<usize>,
    ) -> PyResult<Vec<LorentzVector<Complex<f64>>>> {
        let mut moms: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP + 4]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(l.cast());
        }

        let mut external_momenta: ArrayVec<[LorentzVector<f128::f128>; MAX_LOOP]> = self
            .squared_topology
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];
        let incoming_energy = self.squared_topology.external_momenta
            [..self.squared_topology.n_incoming_momenta]
            .iter()
            .map(|m| m.t)
            .sum();

        let cutkosky_cuts = &self.squared_topology.cutkosky_cuts[cut_index];

        let scaling = squared_topologies::SquaredTopology::find_scaling(
            cutkosky_cuts,
            &external_momenta,
            &moms[..self.squared_topology.n_loops],
            f128::f128::from_f64(incoming_energy).unwrap(),
            self.squared_topology.settings.general.debug,
        )
        .ok_or_else(|| pyo3::exceptions::ValueError::py_err("No scaling could be obtained"))?;

        let max_cuts = self.squared_topology.n_loops + 1;
        let num_cuts = cutkosky_cuts.cuts.len();
        self.squared_topology.evaluate_cut(
            &moms,
            &mut cut_momenta,
            &mut external_momenta,
            &mut rescaled_loop_momenta,
            &mut subgraph_loop_momenta,
            &mut k_def[..max_cuts],
            &mut self.caches_f128,
            &mut Some(&mut self.integrand.event_manager),
            cut_index,
            scaling[1].0,
            scaling[1].1,
            diagram_set,
        );

        let def = k_def[num_cuts - 1..max_cuts - 1]
            .iter()
            .map(|k| k.map(|d| Complex::new(d.re.into(), d.im.into())))
            .collect();
        Ok(def)
    }

    fn parameterize(
        &mut self,
        x: Vec<f64>,
        loop_index: usize,
        e_cm_squared: f64,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::parameterize::<float>(
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
        let (x, jac) = topologies::Topology::inv_parametrize::<float>(
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

#[cfg(feature = "python_api")]
#[pyclass(name = LTD)]
struct PythonLTD {
    topo: topologies::Topology,
    integrand: integrand::Integrand<topologies::Topology>,
    cache: topologies::LTDCache<float>,
    cache_f128: topologies::LTDCache<f128::f128>,
    _dashboard: dashboard::Dashboard,
}

#[cfg(feature = "python_api")]
#[pymethods]
impl PythonLTD {
    #[new]
    fn new(
        topology_file: &str,
        top_name: &str,
        amplitude_file: &str,
        amp_name: &str,
        settings_file: &str,
    ) -> PyResult<PythonLTD> {
        let mut settings = Settings::from_file(settings_file).unwrap();
        let mut topologies = topologies::Topology::from_file(topology_file, &settings).unwrap();
        let mut topo = topologies.remove(top_name).expect("Unknown topology");
        topo.process(true);
        //Skip amplitude with no name is set
        let mut amp = amplitude::Amplitude::default();
        if amp_name == "" {
            settings.general.use_amplitude = false;
        } else {
            let mut amplitudes = amplitude::Amplitude::from_file(amplitude_file).unwrap();
            settings.general.use_amplitude = true;
            amp = amplitudes.remove(amp_name).expect("Unknown amplitude");
            assert_eq!(amp.topology, top_name);
            amp.process(&settings.general);
        }

        topo.amplitude = amp.clone();

        let cache = topologies::LTDCache::<float>::new(&topo);
        let cache_f128 = topologies::LTDCache::<f128::f128>::new(&topo);
        let dashboard = dashboard::Dashboard::minimal_dashboard();
        let integrand = integrand::Integrand::new(
            topo.n_loops,
            topo.clone(),
            settings.clone(),
            false,
            dashboard.status_update_sender.clone(),
            0,
            None,
        );

        Ok(PythonLTD {
            topo,
            integrand,
            cache,
            cache_f128,
            _dashboard: dashboard,
        })
    }

    fn __copy__(&self) -> PyResult<PythonLTD> {
        let dashboard = dashboard::Dashboard::minimal_dashboard();
        let integrand = integrand::Integrand::new(
            self.topo.n_loops,
            self.topo.clone(),
            self.topo.settings.clone(),
            false,
            dashboard.status_update_sender.clone(),
            0,
            None,
        );
        let cache = topologies::LTDCache::<float>::new(&self.topo);
        let cache_f128 = topologies::LTDCache::<f128::f128>::new(&self.topo);
        Ok(PythonLTD {
            topo: self.topo.clone(),
            integrand,
            cache,
            cache_f128,
            _dashboard: dashboard,
        })
    }

    fn evaluate_integrand(&mut self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self
            .integrand
            .evaluate(integrand::IntegrandSample::Flat(1., &x), 1.0, 1);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    fn evaluate(&mut self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self
            .topo
            .evaluate::<float>(integrand::IntegrandSample::Flat(1., &x), &mut self.cache);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    fn evaluate_f128(&mut self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.topo.evaluate::<f128::f128>(
            integrand::IntegrandSample::Flat(1., &x),
            &mut self.cache_f128,
        );
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    fn parameterize(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::parameterize::<float>(
            &x,
            self.topo.e_cm_squared,
            loop_index,
            &self.topo.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn parameterize_f128(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::parameterize::<f128::f128>(
            &x,
            self.topo.e_cm_squared,
            loop_index,
            &self.topo.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn inv_parameterize(
        &self,
        loop_momentum: LorentzVector<f64>,
        loop_index: usize,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::inv_parametrize::<float>(
            &loop_momentum,
            self.topo.e_cm_squared,
            loop_index,
            &self.topo.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn inv_parameterize_f128(
        &self,
        loop_momentum: LorentzVector<f64>,
        loop_index: usize,
    ) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = topologies::Topology::inv_parametrize::<f128::f128>(
            &loop_momentum.cast(),
            self.topo.e_cm_squared,
            loop_index,
            &self.topo.settings,
        );
        Ok((
            x[0].to_f64().unwrap(),
            x[1].to_f64().unwrap(),
            x[2].to_f64().unwrap(),
            jac.to_f64().unwrap(),
        ))
    }

    fn evaluate_cut(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(
                    float::from_f64(l[0].0).unwrap(),
                    float::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    float::from_f64(l[1].0).unwrap(),
                    float::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    float::from_f64(l[2].0).unwrap(),
                    float::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        if self
            .topo
            .populate_ltd_cache(&moms, &mut self.cache)
            .is_err()
        {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = self.topo.settings.general.partial_fractioning_threshold
            >= 0.0
            && moms.len() > 0
            && moms[0].real().spatial_distance()
                > self.topo.settings.general.partial_fractioning_threshold;

        // Prepare numerator
        self.topo.numerator.evaluate_reduced_in_lb(
            &moms,
            0,
            &mut self.cache,
            0,
            true,
            use_partial_fractioning,
        );

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];

        match self.topo.evaluate_cut::<float>(
            &mut moms,
            &self.topo.numerator,
            cut,
            mat,
            &mut self.cache,
            true,
            0,
        ) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.)),
        }
    }

    fn evaluate_cut_f128(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(
                    f128::f128::from_f64(l[0].0).unwrap(),
                    f128::f128::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[1].0).unwrap(),
                    f128::f128::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[2].0).unwrap(),
                    f128::f128::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        if self
            .topo
            .populate_ltd_cache(&moms, &mut self.cache_f128)
            .is_err()
        {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = self.topo.settings.general.partial_fractioning_threshold
            >= 0.0
            && moms.len() > 0
            && moms[0].real().spatial_distance()
                > f128::f128::from_f64(self.topo.settings.general.partial_fractioning_threshold)
                    .unwrap();

        // Prepare numerator
        self.topo.numerator.evaluate_reduced_in_lb(
            &moms,
            0,
            &mut self.cache_f128,
            0,
            true,
            use_partial_fractioning,
        );

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];

        match self.topo.evaluate_cut::<f128::f128>(
            &mut moms,
            &self.topo.numerator,
            cut,
            mat,
            &mut self.cache_f128,
            true,
            0,
        ) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.)),
        }
    }

    fn evaluate_cut_ct(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(
                    float::from_f64(l[0].0).unwrap(),
                    float::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    float::from_f64(l[1].0).unwrap(),
                    float::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    float::from_f64(l[2].0).unwrap(),
                    float::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        if self
            .topo
            .populate_ltd_cache(&moms, &mut self.cache)
            .is_err()
        {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = self.topo.settings.general.partial_fractioning_threshold
            >= 0.0
            && moms.len() > 0
            && moms[0].real().spatial_distance()
                > self.topo.settings.general.partial_fractioning_threshold;

        // Prepare numerator
        self.topo.numerator.evaluate_reduced_in_lb(
            &moms,
            0,
            &mut self.cache,
            0,
            true,
            use_partial_fractioning,
        ); //NOTE: Only necessary when k_vec is changed

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];
        let v = self
            .topo
            .evaluate_cut::<float>(
                &mut moms,
                &self.topo.numerator,
                cut,
                mat,
                &mut self.cache,
                true,
                0,
            )
            .unwrap();
        // get the loop line result from the cache if possible
        let r = 2.0 * self.cache.complex_cut_energies[cut_index];
        let ct = self.topo.counterterm::<float>(
            &moms[..self.topo.n_loops],
            r,
            cut_index,
            &mut self.cache,
        );

        let res = v * (ct + float::one());
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    fn evaluate_cut_ct_f128(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(
                    f128::f128::from_f64(l[0].0).unwrap(),
                    f128::f128::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[1].0).unwrap(),
                    f128::f128::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[2].0).unwrap(),
                    f128::f128::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        if self
            .topo
            .populate_ltd_cache(&moms, &mut self.cache_f128)
            .is_err()
        {
            return Ok((0., 0.));
        }

        let use_partial_fractioning = self.topo.settings.general.partial_fractioning_threshold
            >= 0.0
            && moms.len() > 0
            && moms[0].real().spatial_distance()
                > f128::f128::from_f64(self.topo.settings.general.partial_fractioning_threshold)
                    .unwrap();

        // Prepare numerator
        self.topo.numerator.evaluate_reduced_in_lb(
            &moms,
            0,
            &mut self.cache_f128,
            0,
            true,
            use_partial_fractioning,
        ); //NOTE: Only necessary when k_vec is changed

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];
        let v = self
            .topo
            .evaluate_cut::<f128::f128>(
                &mut moms,
                &self.topo.numerator,
                cut,
                mat,
                &mut self.cache_f128,
                true,
                0,
            )
            .unwrap();
        // get the loop line result from the cache if possible
        let r =
            self.cache_f128.complex_cut_energies[cut_index] * f128::f128::from_f64(2.0).unwrap();
        let ct = self.topo.counterterm::<f128::f128>(
            &moms[..self.topo.n_loops],
            r,
            cut_index,
            &mut self.cache_f128,
        );

        match v * (Complex::new(f128::f128::from_f64(1.0).unwrap(), f128::f128::zero()) + ct) {
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    }

    fn evaluate_amplitude_cut(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(0.0, 0.0),
                Complex::new(l[0].0, l[0].1),
                Complex::new(l[1].0, l[1].1),
                Complex::new(l[2].0, l[2].1),
            ));
        }

        // FIXME: recomputed every time
        if self
            .topo
            .populate_ltd_cache(&moms, &mut self.cache)
            .is_err()
        {
            return Ok((0., 0.));
        }

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];
        match self
            .topo
            .evaluate_amplitude_cut::<float>(&mut moms, cut, mat, &mut self.cache)
            .unwrap()
        {
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    }

    fn evaluate_amplitude_cut_f128(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<(f64, f64)> {
        let mut moms: ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(
                    f128::f128::from_f64(l[0].0).unwrap(),
                    f128::f128::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[1].0).unwrap(),
                    f128::f128::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[2].0).unwrap(),
                    f128::f128::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        if self
            .topo
            .populate_ltd_cache(&moms, &mut self.cache_f128)
            .is_err()
        {
            return Ok((0., 0.));
        }

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];
        match self
            .topo
            .evaluate_amplitude_cut::<f128::f128>(&mut moms, cut, mat, &mut self.cache_f128)
            .unwrap()
        {
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    }

    fn get_loop_momentum_energies(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<Vec<(f64, f64)>> {
        let mut moms: ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(
                    float::from_f64(l[0].0).unwrap(),
                    float::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    float::from_f64(l[1].0).unwrap(),
                    float::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    float::from_f64(l[2].0).unwrap(),
                    float::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        self.topo
            .populate_ltd_cache(&moms, &mut self.cache)
            .unwrap_or_else(|_| {
                println!("On-shell propagator detected");
            });

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];

        self.topo
            .set_loop_momentum_energies::<float>(&mut moms, cut, mat, &self.cache);

        let mut res = Vec::with_capacity(moms.len());
        for l in moms {
            res.push((l.t.re.to_f64().unwrap(), l.t.im.to_f64().unwrap()));
        }

        Ok(res)
    }

    fn get_loop_momentum_energies_f128(
        &mut self,
        loop_momenta: Vec<Vec<(f64, f64)>>,
        cut_structure_index: usize,
        cut_index: usize,
    ) -> PyResult<Vec<(f64, f64)>> {
        let mut moms: ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(f128::f128::zero(), f128::f128::zero()),
                Complex::new(
                    f128::f128::from_f64(l[0].0).unwrap(),
                    f128::f128::from_f64(l[0].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[1].0).unwrap(),
                    f128::f128::from_f64(l[1].1).unwrap(),
                ),
                Complex::new(
                    f128::f128::from_f64(l[2].0).unwrap(),
                    f128::f128::from_f64(l[2].1).unwrap(),
                ),
            ));
        }

        // FIXME: recomputed every time
        self.topo
            .populate_ltd_cache(&moms, &mut self.cache_f128)
            .unwrap_or_else(|_| {
                println!("On-shell propagator detected");
            });

        let mat = &self.topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &self.topo.ltd_cut_options[cut_structure_index][cut_index];

        self.topo
            .set_loop_momentum_energies::<f128::f128>(&mut moms, cut, mat, &self.cache_f128);

        let mut res = Vec::with_capacity(moms.len());
        for l in moms {
            res.push((l.t.re.to_f64().unwrap(), l.t.im.to_f64().unwrap()));
        }

        Ok(res)
    }

    fn deform(
        &mut self,
        loop_momenta: Vec<Vec<f64>>,
    ) -> PyResult<(Vec<(f64, f64, f64)>, f64, f64)> {
        let mut moms = Vec::with_capacity(loop_momenta.len());
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                float::zero(),
                float::from_f64(l[0]).unwrap(),
                float::from_f64(l[1]).unwrap(),
                float::from_f64(l[2]).unwrap(),
            ));
        }

        let (res, jac) = self.topo.deform::<float>(&moms, &mut self.cache);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..self.topo.n_loops].iter() {
            r.push((
                x[1].to_f64().unwrap(),
                x[2].to_f64().unwrap(),
                x[3].to_f64().unwrap(),
            ));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }

    fn deform_f128(
        &mut self,
        loop_momenta: Vec<Vec<f64>>,
    ) -> PyResult<(Vec<(f64, f64, f64)>, f64, f64)> {
        let mut moms = Vec::with_capacity(loop_momenta.len());
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                f128::f128::zero(),
                f128::f128::from_f64(l[0]).unwrap(),
                f128::f128::from_f64(l[1]).unwrap(),
                f128::f128::from_f64(l[2]).unwrap(),
            ));
        }

        let (res, jac) = self.topo.deform::<f128::f128>(&moms, &mut self.cache_f128);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..self.topo.n_loops].iter() {
            r.push((
                x[1].to_f64().unwrap(),
                x[2].to_f64().unwrap(),
                x[3].to_f64().unwrap(),
            ));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }
}
