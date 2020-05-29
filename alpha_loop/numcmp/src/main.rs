
#[macro_use]
extern crate dlopen_derive;

use dlopen::wrapper::{Container, WrapperApi};
use libc::{c_char, c_double, c_int};
use num::Complex;
use itertools::Itertools;

#[derive(WrapperApi)]
pub struct CNumeratorAPI {
    evaluate: unsafe extern "C" fn(p: *const c_double, diag: c_int, conf: c_int, out: *mut c_double),
}

#[derive(WrapperApi)]
pub struct MGNumeratorAPI {
    c_get_numerator: unsafe extern "C" fn(
        p: *const c_double,
        proc_id: *const c_int,
        selected_diagam_left: *const c_int,
        selected_diagram_right: *const c_int,
        ans_re: *mut c_double,
        ans_im: *mut c_double,
    ),
    c_initialise: unsafe extern "C" fn(p: *const c_char),
}

pub fn initialise(api_container: &mut Container<MGNumeratorAPI>, param_card_filename: &str) {
    unsafe {
        let mut a = ['\0' as i8; 512];
        for (xa, c) in a.iter_mut().zip(param_card_filename.chars()) {
            *xa = c as i8;
        }
        api_container.c_initialise(a.as_ptr());
    }
}

pub fn get_form_numerator(
    api_container: &mut Container<CNumeratorAPI>,
    p: &[f64],
    diag: usize,
    conf: usize,
    poly: &mut[f64]
) {
    unsafe {
        api_container.evaluate(&p[0] as *const f64, diag as i32, conf as i32, &mut poly[0] as *mut f64);
    }
}

pub fn get_numerator(
    api_container: &mut Container<MGNumeratorAPI>,
    p: &[f64],
    proc_id: usize,
    left_diagram_id: usize,
    right_diagram_id: usize,
) -> Complex<f64> {
    let proc_id = proc_id as i32;
    let selected_diagram_left = left_diagram_id as i32;
    let selected_diagram_right = right_diagram_id as i32;

    unsafe {
        let mut ans = Complex::default();
        api_container.c_get_numerator(
            &p[0] as *const f64,
            &proc_id as *const i32,
            &selected_diagram_left as *const i32,
            &selected_diagram_right as *const i32,
            &mut ans.re as *mut f64,
            &mut ans.im as *mut f64,
        );
        ans
    }
}

pub fn load() -> Container<MGNumeratorAPI> {
    let mut path = std::env::var("MG_NUMERATOR_PATH")
        .expect("MG_NUMERATOR_PATH needs to be set in the mg_numerator mode.");

    let mut container: Container<MGNumeratorAPI> =
        unsafe { Container::load(&(path.clone() + "lib/libMGnumerators_dynamic.so")) }
            .expect("Could not open library or load symbols");

    path += "/Cards/param_card.dat";
    initialise(&mut container, &path);

    container
}

pub fn load_form() -> Container<CNumeratorAPI> {
    let mut path = std::env::var("MG_NUMERATOR_PATH")
        .expect("MG_NUMERATOR_PATH needs to be set in the mg_numerator mode.");

    let mut container: Container<CNumeratorAPI> =
        unsafe { Container::load(&(path.clone() + "lib/libFORM_numerators.so")) }
            .expect("Could not open library or load symbols");

    container
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
}

fn spatial_dot(a: &[f64], b: &[f64]) -> f64 {
    a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}

fn eval(
    lin: &[&[f64; 4]],
    ext: usize,
    mg_numerator: &mut Container<MGNumeratorAPI>,
    c_numerator: &mut Container<CNumeratorAPI>,
) {
    let mut form_input = vec![];

    for (ii, p) in lin[..ext].iter().enumerate() {
        form_input.push(p[0]);
        form_input.push(0.);
        for j in &lin[ii..ext] {
            form_input.push(dot(*p, *j));
            form_input.push(0.);
        }
    }

    for (ii, p) in lin[ext..].iter().enumerate() {
        form_input.push(p[0]);
        form_input.push(0.);
        for j in &lin[..ext] {
            form_input.push(dot(*p, *j));
            form_input.push(0.);
            form_input.push(spatial_dot(*p, *j));
            form_input.push(0.);
        }

        for (ii, p2) in lin[ext + ii..].iter().enumerate() {
            form_input.push(dot(*p, *p2));
            form_input.push(0.);
            form_input.push(spatial_dot(*p, *p2));
            form_input.push(0.);
        }
    }

    let mut mg_input = vec![];
    for k in lin {
        for kk in *k {
            mg_input.push(*kk);
            mg_input.push(0.);
        }
    }

    let exponent_map: Vec<Vec<usize>> = (0..=8).map(|rank| (0..lin.len() - ext).combinations_with_replacement(rank)).flatten().collect();
    let mut buffer = vec![0.; lin.len() * lin.len() * 4];

    let t = std::time::Instant::now();

    for i in 0..1_000_000 {
        get_form_numerator(c_numerator, &form_input, 0, 0, &mut buffer);
    }
    let t2 = std::time::Instant::now().duration_since(t);
    println!("time FORM for 1M={:#?}",  t2);

    get_form_numerator(c_numerator, &form_input, 0, 0, &mut buffer);

    let mut r = Complex::<f64>::default();
    for (expm, num) in exponent_map.iter().zip(buffer.chunks(2)) {
        //println!("expm={:?}, {} {}", expm, num[0], num[1]);
        r += Complex::new(num[0], num[1]) * expm.iter().fold(1., |prod, i| prod * lin[ext + *i][0]);
    }
    
    println!("=========================");
    println!("FORM={}", r);
    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 1);
    println!("MG #1,1 = {}", r1);
    println!("r={}", r.im / r1.im);
    println!("=========================");
    /*let r = get_form_numerator(c_numerator, &form_input, 1) / 16.;
    println!("FORM #1 ={}", r);
    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 2);
    println!("MG #1,2= {}", r1);
    println!("r={}", r.im / r1.im);
//    let r = get_form_numerator(c_numerator, &form_input, 2) / 16.;
//    println!("FORM #2 ={}", r);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 2, 2);
//    println!("MG #2,2= {}", r1);
//    println!("r={}", r.im / r1.im);
    println!("=========================");

    let t = std::time::Instant::now();

    // FIXME: uncommenting this gives a different result for the mg numerator
    for i in 0..1_000 {
        let r2 = get_numerator(mg_numerator, &mg_input, 0, 1, 1);
//        if r1 != r2 {
//            println!("{} {}", r1, r2);
//        }
    }
    println!("time MG={:#?}", std::time::Instant::now().duration_since(t));

    //println!("r={}", r.im / r1.im);*/
}

fn main() {
    let mut c_numerator = load_form();
    let mut mg_numerator = load();

    let mut p1 = [1., 0., 0., 1.];
    let mut p2 = [1., 0., 0., -1.];
    let mut k1 = [0.1, 0.0, 0.0, 0.0];
    let mut k2 = [0.1, 0.0, 0.0, 0.0];
    let mut k3 = [0.1, 0.0, 0.0, 0.0];
    let mut k4 = [0.1, 0.0, 0.0, 0.0];

    let lin = [&p1, &p2, &k1, &k2];//, &k3, &k4];

    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);

    let mut p1 = [2., 0., 0., 2.];
    let mut p2 = [2., 0., 0., -2.];
    let mut k1 = [0.1, 0.2, 0.3, 0.4];
    let mut k2 = [0.4, 0.2, 0.4, 0.2];
    let mut k3 = [0.4, 0.3, 0.1, 0.5];
    let mut k4 = [0.7, 0.3, 0.5, 0.2];

    let lin = [&p1, &p2, &k1, &k2];//, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);

    let mut p1 = [6., 0., 0., 6.];
    let mut p2 = [6., 0., 0., -6.];
    let mut k1 = [0.3, 0.5, 0.3, 0.];
    let mut k2 = [0.4, 0.2, 0.4, 0.7];
    let mut k3 = [0.11, 0.52, 0.4, 0.3];
    let mut k4 = [0.43, 0.15, 0.2, 0.7];

    let lin = [&p1, &p2, &k1, &k2];//, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);
}
