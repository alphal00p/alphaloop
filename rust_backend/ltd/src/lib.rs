#[macro_use]
extern crate cpython;
extern crate arrayvec;
extern crate dual_num;
extern crate num;
extern crate vector;
use cpython::{exc, PyErr, PyResult};
use std::cell::RefCell;
extern crate cuba;
extern crate nalgebra as na;
extern crate num_traits;

mod ltd;
mod topologies;
mod utils;

use na::Vector3;

// add bindings to the generated python module
py_module_initializer!(ltd, initltd, PyInit_ltd, |py, m| {
    m.add(py, "__doc__", "LTD")?;
    m.add_class::<LTD>(py)?;
    Ok(())
});

py_class!(class LTD |py| {
    data ltd: RefCell<ltd::LTD>;

    def __new__(_cls, name: &str)
    -> PyResult<LTD> {
        let (_loops, e_cm, loop_lines) = topologies::create_topology(name);
        let ltd = ltd::LTD::new(e_cm, loop_lines.clone());

        LTD::create_instance(py, RefCell::new(ltd))
    }

    def evaluate(&self, x: Vec<f64>) -> PyResult<(f64, f64)> {
        let res = self.ltd(py).borrow().evaluate(&x, true);
        Ok((res.re, res.im))
    }

    def evaluate_cut(&self, x: Vec<f64>, cut_structure: Vec<i32>, cut_indices: Vec<i32>) -> PyResult<(f64, f64)> {
        // set the loop line directions for this configuration
        let mut ll = self.ltd(py).borrow().loop_lines.clone();
        for (i, &loopline_cut) in cut_structure.iter().enumerate() {
            if loopline_cut < 0 {
                for (_, sign) in &mut ll[i].loop_momenta {
                    *sign = !*sign;
                }
            }
        }

        let res = self.ltd(py).borrow().evaluate_cut(&x, &ll, &cut_indices);
        Ok((res.re, res.im))
    }

    def deform(&self, loop_momentum: Vec<f64>) -> PyResult<((f64, f64, f64), f64, f64)> {
        let ltd = self.ltd(py).borrow();
        if ltd.loop_lines.len() > 1 {
            return Err(PyErr::new::<exc::TypeError, _>(py, "Only works for the one-loop case"));
        }

        let v = Vector3::from_row_slice(&loop_momentum);
        let (res, jac) = ltd.loop_lines[0].deform(ltd.e_cm, &[v]);
        Ok(((res[0], res[1], res[2]), jac.re, jac.im))
    }
});
