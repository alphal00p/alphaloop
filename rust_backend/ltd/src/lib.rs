#[macro_use]
extern crate cpython;
extern crate arrayvec;
extern crate dual_num;
extern crate itertools;
extern crate num;
extern crate serde;
extern crate serde_yaml;
extern crate vector;
use cpython::PyResult;
use std::cell::RefCell;
extern crate cuba;
extern crate nalgebra as na;
extern crate num_traits;

pub mod cts;
pub mod ltd;
pub mod topologies;
pub mod utils;

type Complex = num::Complex<f64>;
use arrayvec::ArrayVec;
use vector::LorentzVector;

// add bindings to the generated python module
py_module_initializer!(ltd, initltd, PyInit_ltd, |py, m| {
    m.add(py, "__doc__", "LTD")?;
    m.add_class::<LTD>(py)?;
    Ok(())
});

py_class!(class LTD |py| {
    data topo: RefCell<topologies::Topology>;

    def __new__(_cls, config_file: &str, name: &str, deformation: &str)
    -> PyResult<LTD> {
        let mut topologies = topologies::Topology::from_file(config_file, deformation);
        let topo = topologies.remove(name).expect("Unknown topology");

        LTD::create_instance(py, RefCell::new(topo))
    }

    def evaluate(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.topo(py).borrow().evaluate(&x);
        Ok((res.re, res.im))
    }

    def evaluate_cut(&self, loop_momenta: Vec<Vec<(f64,f64)>>, cut_structure_index: usize, cut_index: usize) -> PyResult<(f64, f64)> {
        let topo = self.topo(py).borrow();

        let mut moms = ArrayVec::new();
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(
                Complex::new(0., 0.),
                Complex::new(l[0].0, l[0].1),
                Complex::new(l[1].0, l[1].1),
                Complex::new(l[2].0, l[2].1)));
        }

        let mat = &topo.cb_to_lmb_mat[cut_structure_index];
        let cut = &topo.ltd_cut_options[cut_structure_index][cut_index];

        let res = topo.evaluate_cut(&mut moms, cut, mat);
        Ok((res.re, res.im))
    }

    def deform(&self, loop_momenta: Vec<Vec<f64>>) -> PyResult<(Vec<(f64, f64, f64)>, f64, f64)> {
        let topo = self.topo(py).borrow();

        let mut moms = Vec::with_capacity(loop_momenta.len());
        for l in loop_momenta {
            moms.push(LorentzVector::from_args(0., l[0], l[1], l[2]));
        }

        let (res, jac) = self.topo(py).borrow().deform(&moms);

        let mut r = Vec::with_capacity(moms.len());
        for x in res[..topo.n_loops].iter() {
            r.push((x[1], x[2], x[3]));
        }

        Ok((r, jac.re, jac.im))
    }
});
