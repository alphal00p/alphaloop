#[macro_use] extern crate cpython;
extern crate num;
extern crate vector;
use cpython::{PyResult, PyList, ObjectProtocol, PyErr, exc};
use std::cell::RefCell;

mod deformation;
mod parameterization;

type Complex = num::Complex<f64>;

// add bindings to the generated python module
py_module_initializer!(deformation, initdeformation, PyInit_deformation, |py, m| {
    m.add(py, "__doc__", "Contour deformation")?;
    m.add_class::<Deformation>(py)?;
    Ok(())
});

py_class!(class Deformation |py| {
    data deformer: RefCell<deformation::Deformer>;

    def __new__(_cls, e_cm_sq: f64, qs_py: PyList) -> PyResult<Deformation> {
        let int = deformation::Deformer::new(e_cm_sq).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        let r = Deformation::create_instance(py, RefCell::new(int))?;
        r.set_qs(py, qs_py)?;
        Ok(r)
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

        self.deformer(py).borrow_mut().set_qs(qs);
        Ok(true)
    }

    def deform(&self, momentum: PyList) -> PyResult<Vec<(f64, f64)>> {
        let mut m: Vec<f64> = Vec::with_capacity(4);
        for x in momentum.iter(py) {
            m.push(x.extract(py)?);
        }

        if m.len() != 4 {
            return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
        }

        let res = self.deformer(py).borrow().deform(&vector::LorentzVector::from_vec(m)).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok(vec![(res.t.re, res.t.im), (res.x.re, res.x.im), (res.y.re, res.y.im), (res.z.re, res.z.im)])
    }
});