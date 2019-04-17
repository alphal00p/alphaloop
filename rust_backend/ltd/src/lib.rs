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
extern crate cuba;
extern crate f128;
extern crate nalgebra as na;
extern crate num_traits;

use num_traits::{FromPrimitive, ToPrimitive, Zero};

#[allow(non_camel_case_types)]
#[cfg(feature = "use_f128")]
pub type float = f128::f128;
#[allow(non_camel_case_types)]
#[cfg(not(feature = "use_f128"))]
pub type float = f64;

pub mod cts;
pub mod integrand;
pub mod ltd;
pub mod topologies;
pub mod utils;

use arrayvec::ArrayVec;
use num::Complex;
use serde::Deserialize;
use std::fmt;
use std::fs::File;
use vector::LorentzVector;

#[derive(Debug, Copy, Default, Clone, PartialEq, Deserialize)]
pub struct Scaling {
    pub lambda: f64,
    pub softmin_sigma: f64,
    pub expansion_check: bool,
    pub expansion_threshold: f64,
    pub positive_cut_check: bool,
    pub cut_propagator_check: bool,
    pub non_cut_propagator_check: bool,
    pub skip_hyperboloids: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum DeformationStrategy {
    #[serde(rename = "additive")]
    Additive,
    #[serde(rename = "multiplicative")]
    Multiplicative,
    #[serde(rename = "cutgroups")]
    CutGroups,
    #[serde(rename = "duals")]
    Duals,
    #[serde(rename = "constant")]
    Constant,
    #[serde(rename = "none")]
    None,
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
            "multiplicative" => DeformationStrategy::Multiplicative,
            "cutgroups" => DeformationStrategy::CutGroups,
            "duals" => DeformationStrategy::Duals,
            "constant" => DeformationStrategy::Constant,
            "none" => DeformationStrategy::None,
            _ => panic!("Unknown deformation strategy {}", s),
        }
    }
}

impl fmt::Display for DeformationStrategy {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DeformationStrategy::Additive => write!(f, "additive"),
            DeformationStrategy::Multiplicative => write!(f, "multiplicative"),
            DeformationStrategy::CutGroups => write!(f, "cutgroups"),
            DeformationStrategy::Duals => write!(f, "duals"),
            DeformationStrategy::Constant => write!(f, "constant"),
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
pub struct DeformationMultiplicativeSettings {
    #[serde(rename = "M_ij")]
    pub m_ij: f64,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct DeformationAdditiveSettings {
    pub mode: AdditiveMode,
    pub a_ij: f64,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct DeformationCutGroupsSettings {
    pub mode: AdditiveMode,
    #[serde(rename = "M_ij")]
    pub m_ij: f64,
    pub sigma: f64,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct DeformationSettings {
    pub scaling: Scaling,
    pub overall_scaling: OverallDeformationScaling,
    pub overall_scaling_constant: f64,
    pub multiplicative: DeformationMultiplicativeSettings,
    pub additive: DeformationAdditiveSettings,
    pub cutgroups: DeformationCutGroupsSettings,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
pub enum ParameterizationMode {
    #[serde(rename = "log")]
    Log,
    #[serde(rename = "linear")]
    Linear,
    #[serde(rename = "spherical")]
    Spherical,
}

impl Default for ParameterizationMode {
    fn default() -> ParameterizationMode {
        ParameterizationMode::Spherical
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct ParameterizationSettings {
    pub mode: ParameterizationMode,
    pub shifts: Vec<(f64, f64, f64, f64)>,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct GeneralSettings {
    pub log_file_prefix: String,
    pub screen_log_core: Option<usize>,
    pub log_points_to_screen: bool,
    pub log_stats_to_screen: bool,
    pub deformation_strategy: DeformationStrategy,
    pub cut_filter: Vec<usize>,
    pub topology: String,
    pub unstable_point_warning_percentage: f64,
    pub numerical_threshold: f64,
    pub relative_precision: f64,
    pub absolute_precision: f64,
    pub numerical_instability_check: bool,
    pub return_unstable_point: bool,
    pub integration_statistics: bool,
    pub statistics_interval: usize,
    pub debug: usize,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct IntegratorSettings {
    pub integrator: Integrator,
    pub n_increase: usize,
    pub n_max: usize,
    pub n_start: usize,
    pub n_new: usize,
    pub n_min: usize,
    pub flatness: f64,
    pub seed: i32,
    pub integrated_phase: IntegratedPhase,
}

#[derive(Debug, Clone, Deserialize)]
pub enum Integrator {
    #[serde(rename = "vegas")]
    Vegas,
    #[serde(rename = "suave")]
    Suave,
    #[serde(rename = "cuhre")]
    Cuhre,
}

impl Default for IntegratorSettings {
    fn default() -> IntegratorSettings {
        IntegratorSettings {
            integrator: Integrator::Vegas,
            n_increase: 0,
            n_start: 10000,
            n_max: 10000000,
            n_new: 1000,
            n_min: 2,
            flatness: 50.,
            seed: 1,
            integrated_phase: IntegratedPhase::Real,
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

// add bindings to the generated python module
py_module_initializer!(ltd, initltd, PyInit_ltd, |py, m| {
    m.add(py, "__doc__", "LTD")?;
    m.add_class::<LTD>(py)?;
    Ok(())
});

py_class!(class LTD |py| {
    data topo: RefCell<topologies::Topology>;
    data cache: RefCell<topologies::LTDCache<float>>;

    def __new__(_cls, topology_file: &str, name: &str, settings_file: &str)
    -> PyResult<LTD> {
        let settings = Settings::from_file(settings_file);
        let mut topologies = topologies::Topology::from_file(topology_file, &settings);
        let mut topo = topologies.remove(name).expect("Unknown topology");
        topo.process();
        let cache = topologies::LTDCache::<float>::new(&topo);

        LTD::create_instance(py, RefCell::new(topo), RefCell::new(cache))
    }

    def evaluate(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, res) = self.topo(py).borrow().evaluate::<float>(&x,
            &mut self.cache(py).borrow_mut());
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def parameterize(&self, x: Vec<f64>, loop_index: usize) -> PyResult<(f64, f64, f64, f64)> {
        let (x, jac) = self.topo(py).borrow().parameterize::<float>(&x, loop_index);
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
        topo.compute_complex_cut_energies(&moms, &mut cache);

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        match topo.evaluate_cut::<float>(&mut moms, cut, mat, &mut cache) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.))
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
        topo.compute_complex_cut_energies(&moms, &mut cache);

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        topo.set_loop_momentum_energies::<float>(&mut moms, cut, mat, &cache);

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

        let (res, jac) = topo.deform::<float>(&moms, None, &mut cache);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..topo.n_loops].iter() {
            r.push((x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), x[3].to_f64().unwrap()));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }
});
