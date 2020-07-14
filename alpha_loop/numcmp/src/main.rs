
#[macro_use]
extern crate dlopen_derive;

use dlopen::wrapper::{Container, WrapperApi};
use libc::{c_char, c_double, c_int};
use num::Complex;
use itertools::Itertools;


#[derive(WrapperApi)]
pub struct FORMIntegrandAPI {
    evaluate:
        unsafe extern "C" fn(p: *const c_double, diag: c_int, conf: c_int) -> Complex<c_double>,
}

pub fn get_integrand(
    api_container: &mut Container<FORMIntegrandAPI>,
    p: &[f64],
    diag: usize,
    conf: usize,
) -> Complex<f64> {
    unsafe { api_container.evaluate(&p[0] as *const f64, diag as i32, conf as i32) }
}

pub fn load_integrand() -> Container<FORMIntegrandAPI> {
    let path = std::env::var("MG_NUMERATOR_PATH")
        .expect("MG_NUMERATOR_PATH needs to be set in the mg_numerator mode.");

    let container: Container<FORMIntegrandAPI> =
        unsafe { Container::load(&(path.clone() + "lib/libFORM_integrands.so")) }
            .expect("Could not open library or load symbols");

    container
}


#[derive(WrapperApi)]
pub struct CNumeratorAPI {
    evaluate: unsafe extern "C" fn(p: *const c_double, diag: c_int, conf: c_int, out: *mut c_double) -> c_int,
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
) -> usize {
    unsafe {
        api_container.evaluate(&p[0] as *const f64, diag as i32, conf as i32, &mut poly[0] as *mut f64) as usize
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
    mgleft: usize,
    mgright: usize,
    formid: usize,
    lin: &[&[f64; 4]],
    ext: usize,
    mg_numerator: &mut Container<MGNumeratorAPI>,
    c_numerator: &mut Container<CNumeratorAPI>,
    form_integrand: &mut Container<FORMIntegrandAPI>,
) {
    let mut form_input = vec![];

    for (ii, p) in lin[..ext].iter().enumerate() {
        form_input.push(p[0]);
        form_input.push(0.);
        for j in &lin[ii..ext] {
            form_input.push(dot(*p, *j));
            form_input.push(0.);
            form_input.push(spatial_dot(*p, *j));
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

    for i in 0..100 {
        form_input.push(0.);
    }

    get_integrand(form_integrand, &form_input, formid, 0);

    let exponent_map_k: Vec<Vec<usize>> = (0..=8).map(|rank| (0..1).combinations_with_replacement(rank)).flatten().collect();
    let exponent_map_l: Vec<Vec<usize>> = (0..=8).map(|rank| (1..2).combinations_with_replacement(rank)).flatten().collect();
    let mut buffer = vec![0.; lin.len() * lin.len() * 4];

    let t = std::time::Instant::now();

    /*for i in 0..1_000_000 {
        get_form_numerator(c_numerator, &form_input, 0, 0, &mut buffer);
    }
    let t2 = std::time::Instant::now().duration_since(t);*/
    //println!("time FORM for 1M={:#?}",  t2);

//    get_form_numerator(c_numerator, &form_input, 0, 0, &mut buffer);
//     get_form_numerator(c_numerator, &form_input, 77, 0, &mut buffer);
//    get_form_numerator(c_numerator, &form_input, 35, 0, &mut buffer);

//    get_form_numerator(c_numerator, &form_input, 94, 0, &mut buffer);

//    get_form_numerator(c_numerator, &form_input, 6, 0, &mut buffer);
//    get_form_numerator(c_numerator, &form_input, 10, 0, &mut buffer);

/*
        1,49 -> 39 SIGN
        3,40 -> 51 ZERO IMPRRECISE MG, so prolly OK
        3,50 -> 54 SIGN
        7,49 -> 64 WRONG
        7,50 -> 65 WRONG
        8,50 -> 71 WRONG
        37,49 -> 74 SIGN 
        38,49 -> 79 SIGN
        38,50 -> 80 SIGN
        49,49 -> 81 WRONG
        49,50 -> 82 WRONG
        1,50 -> 40 SIGN
*/

//    get_form_numerator(c_numerator, &form_input, 39, 0, &mut buffer);
//    get_form_numerator(c_numerator, &form_input, 81, 0, &mut buffer);
//    get_form_numerator(c_numerator, &form_input, 51, 0, &mut buffer);

    get_form_numerator(c_numerator, &form_input, formid, 0, &mut buffer);
    println!("FORM={}", Complex::new(buffer[0], buffer[1]));
    println!("=========================");

    let entries = get_form_numerator(c_numerator, &form_input, formid, 4, &mut buffer); // 1 is with k2

    let mut r = Complex::<f64>::default();
    for (expm, num) in exponent_map_k.iter().zip(buffer[..entries * 2].chunks(2)) {
        //println!("expm={:?}, {} {}", expm, num[0], num[1]);
        r += Complex::new(num[0], num[1]) * expm.iter().fold(1., |prod, i| prod * lin[ext + *i][0]);
    }
    
    println!("FORM k energy={}", r);
    println!("=========================");

    let entries = get_form_numerator(c_numerator, &form_input, formid, 1, &mut buffer);
    let mut r = Complex::<f64>::default();
    for (expm, num) in exponent_map_l.iter().zip(buffer[..entries * 2].chunks(2)) {
        //println!("expm={:?}, {} {}", expm, num[0], num[1]);
        r += Complex::new(num[0], num[1]) * expm.iter().fold(1., |prod, i| prod * lin[ext + *i][0]);
    }

    println!("FORM l energy={}", r);
    println!("=========================");

    let r1 = get_numerator(mg_numerator, &mg_input, 0, mgleft, mgright);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 9);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 9, 9);

//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 49);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 49, 49);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 3, 40);
    
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 1);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 38, 42);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 42);

//    let r1 = get_numerator(mg_numerator, &mg_input, 1, 1, 1);

//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 1, 7);
//    let r1 = get_numerator(mg_numerator, &mg_input, 0, 7, 8);
    println!("MG #1,1 = {}", r1);
    println!("r={}", r.norm() / r1.norm());
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
    let mut form_integrand = load_integrand();

    let mut p1 = [1., 0., 0., 1.];
    let mut p2 = [1., 0., 0., -1.];
    let mut k1 = [0.3, 0.0, 0.0, 0.];
    let mut k2 = [0.1, 0.0, 0.0, 0.];
    let mut k3 = [0.2, 0.0, 0.0, 1.0];
    let mut k4 = [0.3, 0.0, 0.0, 0.0];

    let mut p1 = [500., 0., 0., 500.];
    let mut p2 = [500., 0., 0., -500.];
    let mut k1 = [0.1871526174824760E+03,  -0.1257780476763981E+02,   0.3111387284817541E+02,  -0.6301450606135616E+02];
    let mut k2 = [0.2890726872654491E+03,  -0.6872525125230017E+02,  -0.2176061267288446E+03,   0.3947698029543776E+02];
    //let mut k3 = [0.1054303841594690E+03,  -0.7477954838939333E+02,  -0.6876840480238191E+02,   0.2818672644397240E+02];
    //let mut k4 = [0.3423857129250602E+03,   0.1919284080684538E+03,   0.2751076866100379E+03,  -0.6860920754231151E+02];

    

    let lin = [&p1, &p2, &k1, &k2];//, &k3, &k4];//, &k3, &k4];
    eval(1,2,1,&lin, 2, &mut mg_numerator, &mut c_numerator, &mut form_integrand);
    //eval(1,9,1,&lin, 2, &mut mg_numerator, &mut c_numerator);
    //eval(9,9,2,&lin, 2, &mut mg_numerator, &mut c_numerator);

/*
    let mut p1 = [500., 0., 0., 500.];
    let mut p2 = [500., 0., 0., -500.];
    let mut k1 = [0.8635406814378139E+02,  -0.1521338932026180E+02,   0.3763355129491628E+02,  -0.7621872268218542E+02];
    let mut k2 = [0.2801181818093764E+03,  -0.8312611165058223E+02,  -0.2632038567586505E+03,   0.4774908511602658E+02];
    let mut k3 = [100.0*0.11, 100.0*0.52, 100.0*0.4, 100.0*0.3];
    let mut k4 = [100.0*0.43, 100.0*0.15, 100.0*0.2, 100.0*0.7];

    let lin = [&p1, &p2, &k1, &k2, &k3, &k4];//, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);
*/
/*


    let mut p1 = [500., 0., 0., 500.];
    let mut p2 = [500., 0., 0., -500.];
    let mut k1 = [-0.8635406814378139E+02,  0.1521338932026180E+02,  -0.3763355129491628E+02,  0.7621872268218542E+02];
    let mut k2 = [-0.2801181818093764E+03,  0.8312611165058223E+02,  0.2632038567586505E+03,  -0.4774908511602658E+02];
    let mut k3 = [-0.1275225295696605E+03,  0.9044904129599348E+02,  0.8317830770307893E+02,  -0.3409304333925805E+02];
    let mut k4 = [-0.4141300683745435E+03,  -0.2321455649459386E+03,  -0.3327544367808187E+03,  0.8298575185244263E+02];

    let lin = [&p1, &p2, &k1, &k2, &k3, &k4];//, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);

    let mut p1 = [2., 0., 0., 2.];
    let mut p2 = [2., 0., 0., -2.];
    let mut k1 = [0.1, 0.2, 0.3, 0.4];
    let mut k2 = [0.4, 0.2, 0.4, 0.2];
    let mut k3 = [0.4, 0.3, 0.1, 0.5];
    let mut k4 = [0.7, 0.3, 0.5, 0.2];

    let lin = [&p1, &p2, &k1, &k2, &k3, &k4];//, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);

    let mut p1 = [6., 0., 0., 6.];
    let mut p2 = [6., 0., 0., -6.];
    let mut k1 = [0.3, 0.5, 0.3, 0.];
    let mut k2 = [0.4, 0.2, 0.4, 0.7];
    let mut k3 = [0.11, 0.52, 0.4, 0.3];
    let mut k4 = [0.43, 0.15, 0.2, 0.7];

    let mut p1 = [6., 0., 0., 6.];
    let mut p2 = [6., 0., 0., -6.];
    let mut k1 = [0.3, 0.5, 0.3, 0.1];
    let mut k2 = [0.4, 0.2, 0.4, 0.7];
    let mut k3 = [0.11, 0.52, 0.4, 0.3];
    let mut k4 = [0.43, 0.15, 0.2, 0.7];


    //let lin = [&p1, &p2, &k1, &k2, &k3, &k4];
    let lin = [&p1, &p2, &k1, &k2, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);

    let mut p1 = [100.0*6., 0., 0., 100.0*6.];
    let mut p2 = [100.0*6., 0., 0., 100.0*-6.];
    let mut k1 = [100.0*0.3, 100.0*0.5, 100.0*0.3, 100.0*0.1];
    let mut k2 = [100.0*0.4, 100.0*0.2, 100.0*0.4, 100.0*0.7];
    let mut k3 = [100.0*0.11, 100.0*0.52, 100.0*0.4, 100.0*0.3];
    let mut k4 = [100.0*0.43, 100.0*0.15, 100.0*0.2, 100.0*0.7];


    //let lin = [&p1, &p2, &k1, &k2, &k3, &k4];
    let lin = [&p1, &p2, &k1, &k2, &k3, &k4];
    eval(&lin, 2, &mut mg_numerator, &mut c_numerator);

*/
}
