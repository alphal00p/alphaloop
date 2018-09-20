#[macro_use] extern crate cpython;
extern crate num;
extern crate vector;
use cpython::{PyResult, PyList, ObjectProtocol, PyErr, exc};
use std::cell::RefCell;

mod integrands;

type Complex = num::Complex<f64>;

// add bindings to the generated python module
py_module_initializer!(integrand, initintegrand, PyInit_integrand, |py, m| {
    m.add(py, "__doc__", "Integrand evaluator")?;
    m.add_class::<Integrand>(py)?;
    Ok(())
});

py_class!(class Integrand |py| {
    data integrand: RefCell<integrands::Integrand>;

    def __new__(_cls, id: &str, channel: usize, region: usize, mu_sq: f64) -> PyResult<Integrand> {
        let int = integrands::Integrand::new(id, channel, region, mu_sq).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Integrand::create_instance(py, RefCell::new(int))
    }

    def set_qs(&self, qs_py: PyList) -> PyResult<bool> {
        // TODO: support Python LorentzVector?
        let mut qs: Vec<vector::LorentzVector<f64>> = Vec::with_capacity(4);

        for mom in qs_py.iter(py) {
            let mut m: Vec<f64> = Vec::with_capacity(4); 
            for x in mom.extract::<PyList>(py)?.iter(py) {
                m.push(x.extract(py)?);
            }

            if m.len() != 4 {
                return Err(PyErr::new::<exc::ValueError, _>(py, "Momentum does not have 4 components"));
            }

            qs.push(vector::LorentzVector::from_vec(m));
        }

        self.integrand(py).borrow_mut().qs = qs;
        Ok(true)
    }

    def set_channel(&self, channel: usize) -> PyResult<bool> {
        self.integrand(py).borrow_mut().channel = channel;
        Ok(true)
    }

    def evaluate(&self, momentum: PyList) -> PyResult<(f64, f64)> {
        let mut m: Vec<Complex> = Vec::with_capacity(4); 
        for x in momentum.iter(py) {
            m.push(Complex::new(x.getattr(py, "real")?.extract(py)?, x.getattr(py, "imag")?.extract(py)?));
        }

        if m.len() != 4 {
            return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
        }

        let res = self.integrand(py).borrow().evaluate(&vector::LorentzVector::from_vec(m)).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok((res.re, res.im))
    }
});