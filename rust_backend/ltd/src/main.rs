extern crate arrayvec;
extern crate cuba;
extern crate dual_num;
extern crate nalgebra as na;
extern crate num;
extern crate num_traits;
extern crate vector;

use cuba::{CubaIntegrator, CubaVerbosity};

mod ltd;
mod topologies;
mod utils;

#[derive(Debug)]
struct UserData {
    ltd: Vec<ltd::LTD>,
    sample_count: usize,
}

#[inline(always)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    _nvec: usize,
    core: i32,
) -> Result<(), &'static str> {
    let res = user_data.ltd[(core + 1) as usize].evaluate(x, true);
    user_data.sample_count += 1;

    if user_data.sample_count % 100000 == 0 {
        println!("Sample: {:?} {:e}", x, res.re);
    }

    if res.re.is_finite() {
        f[0] = res.re;
    } else {
        println!("Bad point: {:?}", x);
        f[0] = 0.;
    }

    Ok(())
}

fn main() {
    let mut ci = CubaIntegrator::new(integrand);
    let cores = 1;
    ci.set_mineval(10)
        .set_nstart(100000)
        .set_nincrease(0)
        .set_maxeval(100000000)
        .set_epsabs(0.)
        .set_epsrel(1e-15)
        .set_seed(1)
        .set_cores(cores, 1000);

    let (loops, e_cm, loop_lines) = topologies::create_topology("triangle-box");

    let r = ci.vegas(
        3 * loops,
        1,
        CubaVerbosity::Progress,
        0,
        UserData {
            ltd: vec![ltd::LTD::new(e_cm, loop_lines.clone()); cores + 1],
            sample_count: 0,
        },
    );
    println!("{:#?}", r);
}
