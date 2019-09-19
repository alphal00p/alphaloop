#[macro_use]
extern crate cpython;
extern crate num;
extern crate vector;
use cpython::{exc, ObjectProtocol, PyErr, PyList, PyResult};
use std::cell::RefCell;

pub mod integrands;

type Complex = num::Complex<f64>;

// add bindings to the generated python module
py_module_initializer!(integrand, initintegrand, PyInit_integrand, |py, m| {
    m.add(py, "__doc__", "Integrand evaluator")?;
    m.add_class::<Integrand>(py)?;
    Ok(())
});

py_class!(class Integrand |py| {
    data integrand: RefCell<integrands::Integrand>;

    def __new__(_cls, id: &str, channel: usize, region: usize, mu_sq: f64, tau: f64) -> PyResult<Integrand> {
        let int = integrands::Integrand::new(id, channel, region, mu_sq, tau).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Integrand::create_instance(py, RefCell::new(int))
    }

    def set_externals(&self, ext_py: PyList) -> PyResult<bool> {
        // TODO: support Python LorentzVector?
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

        self.integrand(py).borrow_mut().set_externals(ext);
        Ok(true)
    }

    def set_channel(&self, channel: usize) -> PyResult<bool> {
        self.integrand(py).borrow_mut().set_channel(channel);
        Ok(true)
    }

    def set_on_shell_flag(&self, onshell_flag: usize) -> PyResult<bool> {
        self.integrand(py).borrow_mut().set_on_shell_flag(onshell_flag);
        Ok(true)
    }

    def set_region(&self, region: usize) -> PyResult<bool> {
        self.integrand(py).borrow_mut().set_region(region);
        Ok(true)
    }

    def evaluate(&self, loop_momenta: PyList) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
        let mut moms = Vec::with_capacity(2);

        for m_py in loop_momenta.iter(py) {
            let mut m: Vec<Complex> = Vec::with_capacity(4);
            for x in m_py.extract::<PyList>(py)?.iter(py) {
                m.push(Complex::new(x.getattr(py, "real")?.extract(py)?, x.getattr(py, "imag")?.extract(py)?));
            }

            if m.len() != 4 {
                return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
            }

            moms.push(vector::LorentzVector::from_vec(m));
        }

        let res = self.integrand(py).borrow().evaluate(&moms).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok((res.0.re, res.0.im, res.1.re, res.1.im, res.2.re, res.2.im))
    }
});
