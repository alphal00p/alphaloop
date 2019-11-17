#![recursion_limit = "128"]
#[macro_use]
extern crate cpython;
extern crate arrayvec;
extern crate dual_num;
#[macro_use]
extern crate itertools;
extern crate num;
extern crate serde;
extern crate serde_yaml;
extern crate vector;
use cpython::PyResult;
use std::cell::RefCell;
extern crate colored;
extern crate cuba;
extern crate disjoint_sets;
extern crate f128;
extern crate nalgebra as na;
extern crate num_traits;
extern crate rand;

use num_traits::{Float, FloatConst, FromPrimitive, Num, One, ToPrimitive, Zero};
use utils::Signum;
use vector::{Field, RealNumberLike};

pub const MAX_LOOP: usize = 4;

#[allow(non_camel_case_types)]
#[cfg(feature = "use_f128")]
pub type float = f128::f128;
#[allow(non_camel_case_types)]
#[cfg(not(feature = "use_f128"))]
pub type float = f64;

pub trait FloatLike
where
    Self: From<f64>,
    Self: Num,
    Self: FromPrimitive,
    Self: Float,
    Self: Field,
    Self: RealNumberLike,
    Self: num_traits::Signed,
    Self: FloatConst,
    Self: std::fmt::LowerExp,
    Self: num_traits::float::FloatCore,
    Self: 'static,
    Self: Signum,
{
}

impl FloatLike for f64 {}
impl FloatLike for f128::f128 {}

pub mod amplitude;
pub mod cts;
pub mod gamma_chain;
pub mod integrand;
pub mod ltd;
pub mod topologies;
pub mod utils;

use arrayvec::ArrayVec;
use num::Complex;
use serde::Deserialize;
use std::fmt;
use std::fs::File;
pub use vector::LorentzVector;

use cpython::PythonObject;
use cpython::PythonObjectWithCheckedDowncast;
use cpython::{PyFloat, PyList, PyString, PyTuple, Python};

#[derive(Debug, Copy, Default, Clone, PartialEq, Deserialize)]
pub struct Scaling {
    pub lambda: f64,
    pub softmin_sigma: f64,
    pub expansion_check: bool,
    pub expansion_threshold: f64,
    pub branch_cut_m: f64,
    pub source_branch_cut_m: f64,
    pub positive_cut_check: bool,
    pub cut_propagator_check: bool,
    pub non_cut_propagator_check: bool,
    pub skip_hyperboloids: bool,
    pub source_branch_cut_threshold: f64,
    pub source_branch_cut_multiplier: f64,
    pub expansion_check_strategy: ExpansionCheckStrategy,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum DeformationStrategy {
    #[serde(rename = "additive")]
    Additive,
    #[serde(rename = "duals")]
    Duals,
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
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum OverallDeformationScaling {
    #[serde(rename = "constant")]
    Constant,
    #[serde(rename = "linear")]
    Linear,
    #[serde(rename = "sigmoid")]
    Sigmoid,
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
            "duals" => DeformationStrategy::Duals,
            "constant" => DeformationStrategy::Constant,
            "fixed" => DeformationStrategy::Fixed,
            "none" => DeformationStrategy::None,
            _ => panic!("Unknown deformation strategy {}", s),
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
        }
    }
}

impl Default for ExpansionCheckStrategy {
    fn default() -> ExpansionCheckStrategy {
        ExpansionCheckStrategy::Ratio
    }
}

impl fmt::Display for DeformationStrategy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DeformationStrategy::Additive => write!(f, "additive"),
            DeformationStrategy::Duals => write!(f, "duals"),
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
    pub local: bool,
    pub a_ijs: Vec<f64>,
    pub normalize_per_source: bool,
    pub normalisation_of_subspace_components: bool,
    pub normalisation_per_number_of_sources: bool,
    pub include_normal_source: bool,
    pub source_dampening_factor: f64,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct DeformationSettings {
    pub scaling: Scaling,
    pub normalize_on_E_surfaces_m: f64,
    pub overall_scaling: OverallDeformationScaling,
    pub overall_scaling_constant: f64,
    pub lambdas: Vec<f64>,
    pub additive: DeformationAdditiveSettings,
    pub fixed: DeformationFixedSettings,
    pub max_iterations: usize,
    pub stability_threshold: f64,
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
pub struct GeneralSettings {
    pub multi_channeling: bool,
    pub multi_channeling_channel: Option<isize>,
    pub log_file_prefix: String,
    pub res_file_prefix: String,
    pub screen_log_core: Option<usize>,
    pub log_points_to_screen: bool,
    pub log_stats_to_screen: bool,
    pub deformation_strategy: DeformationStrategy,
    pub python_numerator: Option<String>,
    pub mu_uv_sq_re_im: Vec<f64>,
    pub use_ct: bool,
    pub use_collinear_ct: bool,
    pub cut_filter: Vec<usize>,
    pub topology: String,
    #[serde(default)]
    pub use_amplitude: bool,
    pub amplitude: String,
    pub unstable_point_warning_percentage: f64,
    pub numerical_threshold: f64,
    pub relative_precision: f64,
    pub absolute_precision: f64,
    pub numerical_instability_check: bool,
    pub minimal_precision_for_returning_result: f64,
    pub num_digits_different_for_inconsistency: f64,
    pub integration_statistics: bool,
    pub statistics_interval: usize,
    pub debug: usize,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct IntegratorSettings {
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
    pub keep_state_file: bool,
    pub reset_vegas_integrator: bool,
    pub use_only_last_sample: bool,
}

#[derive(Debug, Clone, Deserialize)]
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
            keep_state_file: false,
            reset_vegas_integrator: true,
            use_only_last_sample: false,
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct Settings {
    #[serde(rename = "General")]
    pub general: GeneralSettings,
    #[serde(rename = "Integrator")]
    pub integrator: IntegratorSettings,
    #[serde(rename = "Deformation")]
    pub deformation: DeformationSettings,
    #[serde(rename = "Parameterization")]
    pub parameterization: ParameterizationSettings,
}

impl Settings {
    pub fn from_file(filename: &str) -> Settings {
        let f = File::open(filename).expect("Could not open settings file");
        serde_yaml::from_reader(f).unwrap()
    }
}

pub struct PythonNumerator {
    gil: cpython::GILGuard,
    module: cpython::PyModule,
    buffer: PyList,
    num_loops: usize,
}

impl fmt::Debug for PythonNumerator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "python numerator")
    }
}

impl PythonNumerator {
    pub fn new(module: &str, num_loops: usize) -> PythonNumerator {
        let gil = Python::acquire_gil();

        let (module, buffer) = {
            let py = gil.python();

            let sys = py.import("sys").unwrap();
            PyList::downcast_from(py, sys.get(py, "path").unwrap())
                .unwrap()
                .insert_item(py, 0, PyString::new(py, ".").into_object());

            let module = py.import(module).unwrap();
            let buffer = PyList::new(
                py,
                &(0..num_loops)
                    .map(|_| {
                        PyList::new(
                            py,
                            &[
                                PyList::new(py, &[py.None(), py.None()]).into_object(),
                                PyList::new(py, &[py.None(), py.None()]).into_object(),
                                PyList::new(py, &[py.None(), py.None()]).into_object(),
                                PyList::new(py, &[py.None(), py.None()]).into_object(),
                            ],
                        )
                        .into_object()
                    })
                    .collect::<Vec<_>>(),
            );
            (module, buffer)
        };

        PythonNumerator {
            gil,
            module,
            buffer,
            num_loops,
        }
    }

    pub fn evaluate_numerator<T: FloatLike>(
        &self,
        k_def: &[LorentzVector<Complex<T>>],
    ) -> Complex<T> {
        let py = self.gil.python();

        // convert the vector to tuples
        for (i, k) in k_def[..self.num_loops].iter().enumerate() {
            let tup = PyList::downcast_from(py, self.buffer.get_item(py, i)).unwrap();

            for j in 0..4 {
                let tup1 = PyList::downcast_from(py, tup.get_item(py, j)).unwrap();
                tup1.set_item(
                    py,
                    0,
                    PyFloat::new(py, k[j].re.to_f64().unwrap()).into_object(),
                );
                tup1.set_item(
                    py,
                    1,
                    PyFloat::new(py, k[j].im.to_f64().unwrap()).into_object(),
                );
            }
        }

        let res = self
            .module
            .call(py, "numerator", (&self.buffer,), None)
            .unwrap();
        let res_tup = PyTuple::downcast_from(py, res).unwrap();

        Complex::new(
            T::from_f64(
                res_tup
                    .get_item(py, 0)
                    .extract::<PyFloat>(py)
                    .unwrap()
                    .value(py),
            )
            .unwrap(),
            T::from_f64(
                res_tup
                    .get_item(py, 1)
                    .extract::<PyFloat>(py)
                    .unwrap()
                    .value(py),
            )
            .unwrap(),
        )
    }
}

// add bindings to the generated python module
py_module_initializer!(ltd, initltd, PyInit_ltd, |py, m| {
    m.add(py, "__doc__", "LTD")?;
    m.add_class::<LTD>(py)?;
    Ok(())
});

py_class!(class LTD |py| {
    data topo: RefCell<topologies::Topology>;
    data cache: RefCell<topologies::LTDCache<float>>;
    data cache_f128: RefCell<topologies::LTDCache<f128::f128>>;

    def __new__(_cls, topology_file: &str, top_name: &str, amplitude_file: &str, amp_name: &str, settings_file: &str)
    -> PyResult<LTD> {
        let mut settings = Settings::from_file(settings_file);
        let mut topologies = topologies::Topology::from_file(topology_file, &settings);
        let mut topo = topologies.remove(top_name).expect("Unknown topology");
        topo.process();
        //Skip amplitude with no name is set
        let mut amp= amplitude::Amplitude::default();
        if amp_name == "" {
            settings.general.use_amplitude = false;
        } else {
            let mut amplitudes = amplitude::Amplitude::from_file(amplitude_file);
            settings.general.use_amplitude = true;
            amp=amplitudes.remove(amp_name).expect("Unknown amplitude");
            assert_eq!(amp.topology,top_name);
            amp.process(&settings.general);
        }

        topo.amplitude = amp.clone();

        let cache = topologies::LTDCache::<float>::new(&topo);
        let cache_f128 = topologies::LTDCache::<f128::f128>::new(&topo);

        LTD::create_instance(py, RefCell::new(topo), RefCell::new(cache), RefCell::new(cache_f128))
    }

    def __copy__(&self) -> PyResult<LTD> {
        let topo = self.topo(py).borrow();
        let cache = topologies::LTDCache::<float>::new(&topo);
        let cache_f128 = topologies::LTDCache::<f128::f128>::new(&topo);
        LTD::create_instance(py, RefCell::new(topo.clone()), RefCell::new(cache), RefCell::new(cache_f128))
    }

   def evaluate(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let (_, _k_def, _jac_para, _jac_def, res) = self.topo(py).borrow().evaluate::<float>(&x,
            &mut self.cache(py).borrow_mut(), &None);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    } 

   def evaluate_f128(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let (_, _k_def, _jac_para, _jac_def, res) = self.topo(py).borrow().evaluate::<f128::f128>(&x,
            &mut self.cache_f128(py).borrow_mut(), &None);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }
    
    def parameterize(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = self.topo(py).borrow().parameterize::<float>(&x, loop_index);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def parameterize_f128(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = self.topo(py).borrow().parameterize::<f128::f128>(&x, loop_index);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def inv_parameterize(&self, loop_momentum: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let mom = LorentzVector::from_args(
                0.,
                loop_momentum[0],
                loop_momentum[1],
                loop_momentum[2]);

        let (x, jac) = self.topo(py).borrow().inv_parametrize::<float>(&mom, loop_index);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def inv_parameterize_f128(&self, loop_momentum: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let mom = LorentzVector::from_args(
                0.,
                loop_momentum[0],
                loop_momentum[1],
                loop_momentum[2]);

        let (x, jac) = self.topo(py).borrow().inv_parametrize::<f128::f128>(&mom, loop_index);
        Ok((x[0].to_f64().unwrap(), x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), jac.to_f64().unwrap()))
    }

    def evaluate_cut(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms = ArrayVec::new();
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

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        match topo.evaluate_cut::<float>(&mut moms, cut, mat, &mut cache) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.))
        }
    }

    def evaluate_cut_f128(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms = ArrayVec::new();
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

        match topo.evaluate_cut::<f128::f128>(&mut moms, cut, mat, &mut cache) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.))
        }
    }

    def evaluate_cut_ct(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms = ArrayVec::new();
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

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];
        let v = topo.evaluate_cut::<float>(&mut moms, cut, mat, &mut cache).unwrap();
        // get the loop line result from the cache if possible
        let r = 2.0 * cache.complex_cut_energies[cut_index];
        let ct = topo.counterterm::<float>(&moms[..topo.n_loops], r, cut_index, &mut cache);

        let res = v * (ct + float::one());
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def evaluate_cut_ct_f128(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms = ArrayVec::new();
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
        let v = topo.evaluate_cut::<f128::f128>(&mut moms, cut, mat, &mut cache).unwrap();
        // get the loop line result from the cache if possible
        let r = cache.complex_cut_energies[cut_index] * f128::f128::from_f64(2.0).unwrap();
        let ct = topo.counterterm::<f128::f128>(&moms[..topo.n_loops], r, cut_index, &mut cache);

        match v*(Complex::new(f128::f128::from_f64(1.0).unwrap(), f128::f128::zero())+ct){
            res => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
        }
    } 

   def evaluate_amplitude_cut(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();
        let mut cache = self.cache(py).borrow_mut();

        let mut moms = ArrayVec::new();
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
        let topo = self.topo(py).borrow();
        let mut cache = self.cache_f128(py).borrow_mut();

        let mut moms = ArrayVec::new();
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

        let mut moms = ArrayVec::new();
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

        let mut moms = ArrayVec::new();
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

        let (res, jac) = topo.deform::<float>(&moms, None, None, &mut cache);

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

        let (res, jac) = topo.deform::<f128::f128>(&moms, None, None, &mut cache);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..topo.n_loops].iter() {
            r.push((x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), x[3].to_f64().unwrap()));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }
});
