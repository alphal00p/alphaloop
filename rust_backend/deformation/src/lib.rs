#[macro_use]
extern crate cpython;
extern crate num;
extern crate vector;
use cpython::{exc, PyErr, PyList, PyResult};
use std::cell::RefCell;

pub mod deformation;
pub mod parameterization;
mod utils;

pub const REGION_ALL: usize = 0;
pub const REGION_EXT: usize = 1;
pub const REGION_INT: usize = 2;

type Complex = num::Complex<f64>;

// add bindings to the generated python module
py_module_initializer!(deformation, initdeformation, PyInit_deformation, |py, m| {
    m.add(py, "__doc__", "Contour deformation")?;
    m.add_class::<Deformation>(py)?;
    m.add_class::<Parameterization>(py)?;
    Ok(())
});

py_class!(class Deformation |py| {
    data deformer: RefCell<deformation::Deformer>;

    def __new__(_cls, e_cm_sq: f64, mu_sq: f64, region: usize, qs_py: PyList, masses: Vec<f64>) -> PyResult<Deformation> {
        let int = deformation::Deformer::new(e_cm_sq, mu_sq, region, masses)
                    .map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
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

    def set_region(&self, region: usize) -> PyResult<bool> {
        self.deformer(py).borrow_mut().set_region(region);
        Ok(true)
    }

    def deform(&self, momentum: PyList) -> PyResult<(Vec<(f64, f64)>, f64, f64)> {
        let mut m: Vec<f64> = Vec::with_capacity(4);
        for x in momentum.iter(py) {
            m.push(x.extract(py)?);
        }

        if m.len() != 4 {
            return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
        }

        let (res, jac) = self.deformer(py).borrow().deform(&vector::LorentzVector::from_vec(m)).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok((vec![(res.t.re, res.t.im), (res.x.re, res.x.im), (res.y.re, res.y.im), (res.z.re, res.z.im)], jac.re, jac.im))
    }

    // return analytical and numerical jacobian
    def numerical_jacobian(&self, momentum: PyList) -> PyResult<(f64, f64, f64, f64)> {
        let mut m: Vec<f64> = Vec::with_capacity(4);
        for x in momentum.iter(py) {
            m.push(x.extract(py)?);
        }

        if m.len() != 4 {
            return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
        }

        let k = vector::LorentzVector::from_vec(m);
        let (center, jac_ana) = self.deformer(py).borrow().deform(&k).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        let jac = self.deformer(py).borrow().numerical_jacobian(&k, &center);
        Ok((jac_ana.re, jac_ana.im, jac.re, jac.im))
    }

    def set_external_momenta(&self, ext_py: PyList) -> PyResult<bool> {
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

        self.deformer(py).borrow_mut().set_external_momenta(ext);
        Ok(true)
    }

    def deform_two_loops(&self, id: usize, momenta: PyList) -> PyResult<(Vec<(f64, f64)>, f64, f64, f64, f64)> {
        let mut qs: Vec<vector::LorentzVector<f64>> = Vec::with_capacity(4);

        for mom in momenta.iter(py) {
            let mut m: Vec<f64> = Vec::with_capacity(4);
            for x in mom.extract::<PyList>(py)?.iter(py) {
                m.push(x.extract(py)?);
            }

            if m.len() != 4 {
                return Err(PyErr::new::<exc::ValueError, _>(py, "Momentum does not have 4 components"));
            }

            qs.push(vector::LorentzVector::from_vec(m));
        }

        let (k, l) = self.deformer(py).borrow_mut().deform_two_loops(id, &qs[0], &qs[1]).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;

        // compute jacobian
        let jac = self.deformer(py).borrow_mut().numerical_jacobian_two_loops(id, &qs[0], &qs[1], (&k, &l));

        //let jac_central = self.deformer(py).borrow_mut().numerical_jacobian_center_two_loops(id, &qs[0], &qs[1], (&k, &l)).0;
        let jac_central = Complex::new(0., 0.);

        Ok((vec![(k.t.re, k.t.im), (k.x.re, k.x.im), (k.y.re, k.y.im), (k.z.re, k.z.im), (l.t.re, l.t.im), (l.x.re, l.x.im), (l.y.re, l.y.im), (l.z.re, l.z.im)], jac.re, jac.im,
            jac_central.re, jac_central.im))
    }

});

py_class!(class Parameterization |py| {
    data parameterizer: RefCell<parameterization::Parameterizer>;

    def __new__(_cls, e_cm_sq: f64, region: usize, channel: usize, qs_py: PyList) -> PyResult<Parameterization> {
        let int = parameterization::Parameterizer::new(e_cm_sq, region, channel).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        let r = Parameterization::create_instance(py, RefCell::new(int))?;
        r.set_qs(py, qs_py)?;
        Ok(r)
    }

    def set_mode(&self, mode: &str) -> PyResult<bool> {
        self.parameterizer(py).borrow_mut().set_mode(mode).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok(true)
    }

    def set_region(&self, region: usize) -> PyResult<bool> {
        self.parameterizer(py).borrow_mut().set_region(region);
        Ok(true)
    }

    def set_channel(&self, channel: usize) -> PyResult<bool> {
        self.parameterizer(py).borrow_mut().set_channel(channel);
        Ok(true)
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

        self.parameterizer(py).borrow_mut().set_qs(qs);
        Ok(true)
    }

    def map(&self, momentum: PyList) -> PyResult<(Vec<f64>, f64)> {
        let mut m: Vec<f64> = Vec::with_capacity(4);
        for x in momentum.iter(py) {
            m.push(x.extract(py)?);
        }

        if m.len() != 4 {
            return Err(PyErr::new::<exc::ValueError, _>(py, "Loop momentum does not have 4 components"));
        }

        let (k, jac) = self.parameterizer(py).borrow().map(&vector::LorentzVector::from_vec(m)).map_err(|m| PyErr::new::<exc::ValueError, _>(py, m))?;
        Ok((vec![k.t, k.x, k.y, k.z], jac))
    }
});
