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

use serde::Deserialize;
use std::fmt;
use std::fs::File;

pub type Complex = num::Complex<float>;
use arrayvec::ArrayVec;
use vector::LorentzVector;

#[derive(Debug, Copy, Clone, PartialEq, Deserialize)]
pub enum DeformationStrategy {
    #[serde(rename = "additive")]
    Additive,
    #[serde(rename = "multiplicative")]
    Multiplicative,
    #[serde(rename = "none")]
    None,
}

impl From<&str> for DeformationStrategy {
    fn from(s: &str) -> Self {
        match s {
            "additive" => DeformationStrategy::Additive,
            "multiplicative" => DeformationStrategy::Multiplicative,
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
pub struct DeformationSettings {
    pub lambda: f64,
    pub expansion_threshold: f64,
    pub multiplicative: DeformationMultiplicativeSettings,
    pub additive: DeformationAdditiveSettings,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct ParameterizationSettings {
    pub rescaling: f64,
    pub shift: (f64, f64, f64),
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct GeneralSettings {
    pub log_file_prefix: String,
    pub log_to_screen: bool,
    pub deformation_strategy: DeformationStrategy,
    pub topology: String,
    pub numerical_threshold: f64,
    pub relative_precision: f64,
    pub numerical_instability_check: bool,
    pub integration_statistics: bool,
    pub statistics_interval: usize,
    pub debug: usize,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct IntegratorSettings {
    pub integrator: String,
    pub n_increase: usize,
    pub n_max: usize,
    pub n_start: usize,
    pub integrated_phase: IntegratedPhase,
}

impl Default for IntegratorSettings {
    fn default() -> IntegratorSettings {
        IntegratorSettings {
            integrator: "vegas".to_owned(),
            n_increase: 0,
            n_start: 10000,
            n_max: 10000000,
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
        let f = File::open(filename).unwrap();
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

    def __new__(_cls, topology_file: &str, name: &str, settings_file: &str)
    -> PyResult<LTD> {
        let settings = Settings::from_file(settings_file);
        let mut topologies = topologies::Topology::from_file(topology_file, &settings);
        let topo = topologies.remove(name).expect("Unknown topology");

        LTD::create_instance(py, RefCell::new(topo))
    }

    def evaluate(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let (_, _k_def_rot, _jac_para_rot, _jac_def_rot, res) = self.topo(py).borrow_mut().evaluate(&x);
        Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap()))
    }

    def evaluate_cut(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();

        let mut moms = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(float::from_f64(l[0].0).unwrap(), float::from_f64(l[0].1).unwrap()),
                Complex::new(float::from_f64(l[1].0).unwrap(), float::from_f64(l[1].1).unwrap()),
                Complex::new(float::from_f64(l[2].0).unwrap(), float::from_f64(l[2].1).unwrap())));
        }

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        match topo.evaluate_cut(&mut moms, cut, mat) {
            Ok(res) => Ok((res.re.to_f64().unwrap(), res.im.to_f64().unwrap())),
            Err(_) => Ok((0., 0.))
        }
    }

    def get_loop_momentum_energies(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<Vec<(f64, f64)>> {
        let topo = self.topo(py).borrow();

        let mut moms = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(float::zero(), float::zero()),
                Complex::new(float::from_f64(l[0].0).unwrap(), float::from_f64(l[0].1).unwrap()),
                Complex::new(float::from_f64(l[1].0).unwrap(), float::from_f64(l[1].1).unwrap()),
                Complex::new(float::from_f64(l[2].0).unwrap(), float::from_f64(l[2].1).unwrap())));
        }

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        topo.set_loop_momentum_energies(&mut moms, cut, mat);

        let mut res = Vec::with_capacity(moms.len());
        for l in moms {
            res.push((l.t.re.to_f64().unwrap(), l.t.im.to_f64().unwrap()));
        }

        Ok(res)
    }

    def deform(&self, loop_momenta: Vec<Vec<f64>>) -> PyResult<(Vec<(f64, f64, f64)>, f64, f64)> {
        let topo = self.topo(py).borrow();

        let mut moms = Vec::with_capacity(loop_momenta.len());
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(float::zero(),
                float::from_f64(l[0]).unwrap(),
                float::from_f64(l[1]).unwrap(),
                float::from_f64(l[2]).unwrap()));
        }

        let (res, jac) = self.topo(py).borrow().deform(&moms);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..topo.n_loops].iter() {
            r.push((x[1].to_f64().unwrap(), x[2].to_f64().unwrap(), x[3].to_f64().unwrap()));
        }

        Ok((r, jac.re.to_f64().unwrap(), jac.im.to_f64().unwrap()))
    }
});
