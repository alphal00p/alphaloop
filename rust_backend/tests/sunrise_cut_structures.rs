extern crate arrayvec;
extern crate f128;
extern crate ltd;
extern crate num;
extern crate num_traits;
extern crate rand;
extern crate vector;

use arrayvec::ArrayVec;
use ltd::topologies::{LTDCache, Topology};
use ltd::utils::test_utils::{get_test_topology, numeriacal_eq};
use ltd::MAX_LOOP;
use num::Complex;
use num_traits::{NumCast, One, ToPrimitive, Zero};
use rand::prelude::*;
use lorentz_vector::LorentzVector;

type TestFloat = f128::f128;
//type TestFloat = f64;

fn cross_check<'a>(
    topo: &mut Topology,
    x: &'a [f64],
    cache: &mut LTDCache<TestFloat>,
) -> (Complex<TestFloat>, Complex<TestFloat>, Complex<TestFloat>) {
    // Parameterize
    let mut k = [LorentzVector::default(); MAX_LOOP];
    let mut jac_para = TestFloat::one();
    for i in 0..topo.n_loops {
        // set the loop index to i + 1 so that we can also shift k
        let (l_space, jac) = Topology::parameterize::<TestFloat>(
            &x[i * 3..(i + 1) * 3],
            topo.e_cm_squared,
            i,
            &topo.settings,
        );

        // there could be some rounding here
        let rot = topo.rotation_matrix;
        k[i] = LorentzVector::from_args(
            TestFloat::zero(),
            <TestFloat as NumCast>::from(rot[0][0]).unwrap() * l_space[0]
                + <TestFloat as NumCast>::from(rot[0][1]).unwrap() * l_space[1]
                + <TestFloat as NumCast>::from(rot[0][2]).unwrap() * l_space[2],
            <TestFloat as NumCast>::from(rot[1][0]).unwrap() * l_space[0]
                + <TestFloat as NumCast>::from(rot[1][1]).unwrap() * l_space[1]
                + <TestFloat as NumCast>::from(rot[1][2]).unwrap() * l_space[2],
            <TestFloat as NumCast>::from(rot[2][0]).unwrap() * l_space[0]
                + <TestFloat as NumCast>::from(rot[2][1]).unwrap() * l_space[1]
                + <TestFloat as NumCast>::from(rot[2][2]).unwrap() * l_space[2],
        );
        jac_para *= jac;
    }
    let mut k_def: ArrayVec<[LorentzVector<Complex<TestFloat>>; MAX_LOOP]> = (0..topo.n_loops)
        .map(|i| k[i].map(|x| Complex::new(x, TestFloat::zero())))
        .collect();
    topo.populate_ltd_cache(&k_def[..topo.n_loops], cache)
        .unwrap();

    // Compute using standard LTD
    topo.settings.general.partial_fractioning_threshold = -1.0;
    let topo2 = topo.clone();
    let mut start = std::time::Instant::now();
    let result_ltd = topo2
        .evaluate_all_dual_integrands(&mut k_def[..topo.n_loops], cache)
        .unwrap();
    let dt_ltd = start.elapsed();
    println!("\t#LTD in {:?}", dt_ltd);

    // Compute using partial fractioning
    topo.settings.general.partial_fractioning_threshold = 1e-99;
    let topo2 = topo.clone();
    start = std::time::Instant::now();
    topo2
        .numerator
        .evaluate_reduced_in_lb(&k_def, 0, cache, 0, true, true);

    let result_pf1 = topo2
        .evaluate_all_dual_integrands(&mut k_def[..topo.n_loops], cache)
        .unwrap();
    let dt_pf1 = start.elapsed();
    println!("\t#PF1 in {:?}", dt_pf1);

    // Compute using partial fractioning MultiLoop
    let start = std::time::Instant::now();
    topo.numerator
        .evaluate_reduced_in_lb(&k_def, 0, cache, 0, true, true);
    for prop_pow in cache.propagator_powers.iter_mut() {
        *prop_pow = 0;
    }
    for ll in topo.loop_lines.iter() {
        for p in ll.propagators.iter() {
            cache.propagator_powers[p.id] = p.power;
        }
    }

    let result_pf2 = topo.partial_fractioning_multiloops.evaluate(
        &topo.loop_lines,
        &topo.propagator_id_to_ll_id,
        cache,
    );
    let dt_pf2 = start.elapsed();
    println!(
        "\t#PF2 in {:?} :: PF2/LTD = {:?}",
        dt_pf2,
        dt_pf2.as_nanos() as f64 / dt_ltd.as_nanos() as f64
    );
    (result_ltd, result_pf1, result_pf2)
}

/*    BEGIN TESTS    */
mod cut_structures {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[allow(non_snake_case)]
    fn test_multiloops(name: &str, pt: &mut [f64], coeffs: &[(f64, f64)]) {
        let mut topo = get_test_topology("./tests/topologies/sunrise_cut_structures.yaml", name);
        topo.settings.general.partial_fractioning_threshold = 1e-99;
        topo.numerator_tensor_coefficients.clear();
        // Overwrite coefficients
        for coeff in coeffs.iter() {
            topo.numerator_tensor_coefficients.push(*coeff);
        }
        //topo.numerator_tensor_coefficients.push((1.0, 0.0));
        topo.process(true);

        // map the point back from momentum-space to the unit hypercube
        for (i, x) in pt.chunks_exact_mut(3).enumerate() {
            let r = Topology::inv_parametrize::<TestFloat>(
                &LorentzVector::from_args(0., x[0], x[1], x[2]).cast(),
                topo.e_cm_squared,
                i,
                &topo.settings,
            );
            x[0] = TestFloat::to_f64(&r.0[0]).unwrap();
            x[1] = TestFloat::to_f64(&r.0[1]).unwrap();
            x[2] = TestFloat::to_f64(&r.0[2]).unwrap();
        }
        // Compute values using e-fractioning the standard ltd sum over cuts
        let mut cache = LTDCache::<TestFloat>::new(&topo);
        println!("\nReport values for {}", name);
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);
        println!("\tresult_ltd: {:?}", result_ltd);
        println!("\tresult_pf1: {:?}", result_pf1);
        println!("\tresult_pf2: {:?}", result_pf2);

        if topo.n_loops > 0 {
            let check_1 = !numeriacal_eq(result_ltd, result_pf2);
            let check_2 = !numeriacal_eq(result_pf1, result_pf2);
            if check_1 {
                assert_eq!(result_ltd, result_pf2);
            }
            if check_2 {
                assert_eq!(result_pf1, result_pf2);
            }
        } else {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
        }
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_0() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_0", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_1() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_1", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_2() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_2", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_3() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_3", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_4() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_4", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_5() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_5", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_6() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_6", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_7() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_7", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_8() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_8", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_9() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_9", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_10() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_10", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn sunrise_CS_11() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let mut pt: Vec<f64> = (0..3 * 2)
            .map(|_| (1.0 - 1e-8) * (1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        let coeffs = vec![(1., 0.)];
        test_multiloops("sunrise_CS_11", &mut pt, &coeffs);
    }
}
