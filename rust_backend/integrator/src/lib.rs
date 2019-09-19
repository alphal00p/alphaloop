//#[macro_use]
//extern crate cpython;
extern crate cuba;
extern crate deformation;
extern crate integrand;
extern crate num;
extern crate vector;

#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_yaml;

//use cpython::{exc, PyErr, PyList, PyResult};
//use std::cell::RefCell;

pub const REGION_ALL: usize = 0;
pub const REGION_EXT: usize = 1;
pub const REGION_INT: usize = 2;

type Complex = num::Complex<f64>;

pub mod aggregator;
pub mod evaluator;

/*

use evaluator::Evaluator;
use aggregator::Aggregator;


// add bindings to the generated python module
py_module_initializer!(integrator, initintegrator, PyInit_integrator, |py, m| {
    m.add(py, "__doc__", "Integrator / aggregrator")?;
    m.add_class::<Integrator>(py)?;
    Ok(())
});

py_class!(class Integrator |py| {
    data integrator: RefCell<Aggregator>;

    def __new__(_cls, name: &str, do_regions: bool, do_multichanneling: bool, e_cm_sq: f64, mu_sq: f64, ext_py: PyList) -> PyResult<Integrator> {
        let mut ext: Vec<vector::LorentzVector<f64>> = Vec::with_capacity(4);

        for mom in ext_py.iter(py) {
            let mut m: Vec<f64> = Vec::with_capacity(4);
            for x in mom.extract::<PyList>(py)?.iter(py) {
                m.push(x.extract(py)?);
            }

            if m.len() != 4 {
                return Err(PyErr::new::<exc::ValueError, _>(py, "Momentum does not have 4 components"));
            }

            ext.push(vector::LorentzVector::from_vec(m));
        }

        let i = Evaluator::new("box1L_direct_integration", e_cm_sq, mu_sq, external_momenta);

        let settings = IntegrationSettings {
        cores,
        samples,
        param_mode: "linear".to_owned(),
    };

    let i = match matches.value_of("topology").unwrap() {
        "box" => box_integrator(loops),
        "triangle" => double_triangle_integrator(),
        "trianglebox" => triangle_box_integrator(),
        _ => unreachable!(),
    };

    let mut aggregator = Aggregator::new(loops, regions, multichanneling, i, settings);
    let r = aggregator.aggregate();

        let int = integrator::Integrator::new(name, do_regions, do_multichanneling, e_cm_sq, mu_sq, ext)
                    .map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok(Integrator::create_instance(py, RefCell::new(int))?)
    }

    def evaluate(&self, momentum: PyList) -> PyResult<(f64, f64)> {
        let mut m: Vec<f64> = Vec::with_capacity(4);
        for x in momentum.iter(py) {
            m.push(x.extract(py)?);
        }

        if m.len() != 4 {
            return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
        }

        let r = self.integrator(py).borrow_mut().evaluate(&vector::LorentzVector::from_vec(m));
        Ok((r.re, r.im))
    }
});
*/
