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
mod dot_propagators_1loop {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    fn test_hexagon_linear(
        pows: [usize; 6],
        pt: &mut [f64; 3],
        coeffs: &[(f64, f64); 5],
        negative_cut: bool,
    ) {
        let mut topo = get_test_topology("./tests/topologies/bubble_numerator.yaml", "1l");
        topo.settings.general.partial_fractioning_threshold = 1e-99;
        // Overwrite powers
        for i in 0..6 {
            topo.loop_lines[0].propagators[i].power = pows[i];
        }
        // Overwrite coefficients
        for i in 0..5 {
            topo.numerator_tensor_coefficients[i] = coeffs[i];
        }
        if negative_cut {
            topo.ltd_cut_structure[0][0] = -1;
        } else {
            topo.ltd_cut_structure[0][0] = 1;
        }
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
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);

        println!(
            "ltd = {:?}\npf1 = {:?}\npf2 = {:?}",
            result_ltd, result_pf1, result_pf2
        );
        if topo.n_loops == 1 {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
            if !numeriacal_eq(result_pf1, result_pf2) {
                assert_eq!(result_pf1, result_pf2);
            }
        } else {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
        }
    }

    fn test_hexagon_rank5(
        pows: [usize; 6],
        pt: &mut [f64; 3],
        coeffs: &[(f64, f64)],
        negative_cut: bool,
    ) {
        // Check if there are the right number of coefficients
        assert_eq!(coeffs.len(), 126);
        // Call topology
        let mut topo = get_test_topology("./tests/topologies/bubble_numerator.yaml", "1l");
        topo.settings.general.partial_fractioning_threshold = 1e-99;
        // Overwrite powers
        for i in 0..6 {
            topo.loop_lines[0].propagators[i].power = pows[i];
        }
        // Overwrite coefficients
        for (i, coeff) in coeffs.iter().enumerate() {
            topo.numerator_tensor_coefficients[i] = *coeff;
        }
        if negative_cut {
            topo.ltd_cut_structure[0][0] = -1;
        } else {
            topo.ltd_cut_structure[0][0] = 1;
        }
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
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);

        println!(
            "ltd = {:?}\npf1 = {:?}\npf2 = {:?}",
            result_ltd, result_pf1, result_pf2
        );
        if topo.n_loops == 1 {
            if !numeriacal_eq(result_pf1, result_pf2) {
                assert_eq!(result_pf1, result_pf2);
            }
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
        } else {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
        }
    }

    #[test]
    fn bubble_linear_0_1_energy() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [0, 1, 0, 0, 0, 0];
        let coeffs = [(0.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_2_0_energy() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [2, 0, 0, 0, 0, 0];
        let coeffs = [(0.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_1_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1, 0, 0, 0, 0];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_1_1_no_energy() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1, 0, 0, 0, 0];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (9.0, 3.0), (2.0, 12.0), (7.0, -1.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_1_2_energy() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 2, 0, 0, 0, 0];
        let coeffs = [
            (5.0, 7.0),
            (10.0, 0.0),
            (9.0, 3.0),
            (2.0, 12.0),
            (7.0, -1.0),
        ];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_1_1_energy() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1, 0, 0, 0, 0];
        let coeffs = [(0.0, 0.0), (10.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_3_4_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [3, 4, 0, 0, 0, 0];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_3_4_no_energy() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [3, 4, 0, 0, 0, 0];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (9.0, 3.0), (2.0, 12.0), (7.0, -1.0)];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_1_0_0_0_0_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [1, 0, 0, 0, 0, 1];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_1_0_1_0_0_0_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [1, 0, 1, 0, 0, 0];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_1_0_0_1_0_0_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [1, 0, 0, 1, 0, 0];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_0_1_0_0_1_0_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [0, 1, 0, 0, 1, 0];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_0_0_1_0_0_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [0, 0, 1, 0, 0, 1];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_1_1_1_1_1_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [1, 1, 1, 1, 1, 1];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_linear_2_1_2_2_1_2_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [2, 1, 2, 2, 1, 2];
        test_hexagon_linear(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn bubble_linear_1_1_energy_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (10.0, 5.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [1, 1, 0, 0, 0, 0];
        test_hexagon_linear(pows, &mut pt, &coeffs, true);
    }

    #[test]
    fn bubble_linear_2_1_energy_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (10.0, 5.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [2, 1, 0, 0, 0, 0];
        test_hexagon_linear(pows, &mut pt, &coeffs, true);
    }

    #[test]
    fn bubble_linear_2_2_energy_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let coeffs = [(5.0, 7.0), (10.0, 5.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let pows = [2, 2, 0, 0, 0, 0];
        test_hexagon_linear(pows, &mut pt, &coeffs, true);
    }

    #[test]
    fn box_rank5_1_1_1_1_dense() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1, 1, 1, 0, 0];
        // Create the 126 random coefficient needed for the rank 5 test at one loops
        let mut rng = rand::thread_rng();
        let coeffs: Vec<(f64, f64)> = (0..126)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        test_hexagon_rank5(pows, &mut pt, &coeffs, false);
    }

    #[test]
    fn hexagon_rank5_1_2_1_1_2_1_dense() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1, 1, 1, 1, 1];
        // Create the 126 random coefficient needed for the rank 5 test at one loops
        let mut rng = rand::thread_rng();
        let coeffs: Vec<(f64, f64)> = (0..126)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        test_hexagon_rank5(pows, &mut pt, &coeffs, false);
    }
}

mod dot_propagators_2loops {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[allow(non_snake_case)]
    fn test_2L_2PROPxLL_linear(pows: [usize; 6], pt: &mut [f64; 6], coeffs: &[(f64, f64); 9]) {
        let mut topo = get_test_topology(
            "./tests/topologies/2L_2PROPxLL_numerator.yaml",
            "2L_2PROPxLL",
        );
        topo.settings.general.partial_fractioning_threshold = 1e-99;
        // Overwrite powers
        topo.loop_lines[0].propagators[0].power = pows[0];
        topo.loop_lines[0].propagators[1].power = pows[1];
        topo.loop_lines[1].propagators[0].power = pows[2];
        topo.loop_lines[1].propagators[1].power = pows[3];
        topo.loop_lines[2].propagators[0].power = pows[4];
        topo.loop_lines[2].propagators[1].power = pows[5];
        // Overwrite coefficients
        for (i, coeff) in coeffs.iter().enumerate() {
            topo.numerator_tensor_coefficients[i] = *coeff;
        }
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
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);
        println!("Energies: {:?}", cache.complex_cut_energies);
        println!("ltd = {:?}\npf2 = {:?}", result_ltd, result_pf2);

        if topo.n_loops == 1 {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
            if !numeriacal_eq(result_pf1, result_pf2) {
                assert_eq!(result_pf1, result_pf2);
            }
        } else {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
        }
    }

    #[allow(non_snake_case)]
    fn test_2L_2PROPxLL_rank5(pows: [usize; 6], pt: &mut [f64; 6], coeffs: &[(f64, f64)]) {
        // Check if we have all the coefficients for the rank 5 polynomial
        assert_eq!(coeffs.len(), 1287);

        let mut topo = get_test_topology(
            "./tests/topologies/2L_2PROPxLL_numerator.yaml",
            "2L_2PROPxLL",
        );
        topo.settings.general.partial_fractioning_threshold = 1e-99;
        // Overwrite powers
        topo.loop_lines[0].propagators[0].power = pows[0];
        topo.loop_lines[0].propagators[1].power = pows[1];
        topo.loop_lines[1].propagators[0].power = pows[2];
        topo.loop_lines[1].propagators[1].power = pows[3];
        topo.loop_lines[2].propagators[0].power = pows[4];
        topo.loop_lines[2].propagators[1].power = pows[5];
        // Overwrite coefficients
        for (i, coeff) in coeffs.iter().enumerate() {
            topo.numerator_tensor_coefficients[i] = *coeff;
        }
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
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);
        println!("Energies: {:?}", cache.complex_cut_energies);
        println!("ltd = {:?}\npf2 = {:?}", result_ltd, result_pf2);

        if topo.n_loops == 1 {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
            if !numeriacal_eq(result_pf1, result_pf2) {
                assert_eq!(result_pf1, result_pf2);
            }
        } else {
            if !numeriacal_eq(result_ltd, result_pf2) {
                assert_eq!(result_ltd, result_pf2);
            }
        }
    }

    #[test]
    fn sunrise_linear_1_1_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 1, 0, 1, 0];
        let coeffs = [
            (5.0, 1.0),
            (0.0, 0.0),
            (2.0, 7.0),
            (3.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
            (4.0, -2.0),
            (6.0, -5.0),
            (7.0, 11.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_0_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 0, 0, 1, 0];
        let coeffs = [
            (5.0, 1.0),
            (0.0, 0.0),
            (2.0, 7.0),
            (3.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
            (4.0, -2.0),
            (6.0, -5.0),
            (7.0, 11.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_3_0_5_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [3, 0, 0, 0, 5, 0];
        let coeffs = [
            (5.0, 1.0),
            (0.0, 0.0),
            (2.0, 7.0),
            (3.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
            (4.0, -2.0),
            (6.0, -5.0),
            (7.0, 11.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_1_0_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 1, 0, 0, 0];
        let coeffs = [
            (5.0, 1.0),
            (0.0, 0.0),
            (2.0, 7.0),
            (3.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
            (4.0, -2.0),
            (6.0, -5.0),
            (7.0, 11.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_4_3_0_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [4, 0, 3, 0, 0, 0];
        let coeffs = [
            (5.0, 1.0),
            (0.0, 0.0),
            (2.0, 7.0),
            (3.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
            (4.0, -2.0),
            (6.0, -5.0),
            (7.0, 11.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_4_3_2_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [4, 0, 3, 0, 2, 0];
        let coeffs = [
            (5.0, 1.0),
            (0.0, 0.0),
            (2.0, 7.0),
            (3.0, 1.0),
            (4.0, 3.0),
            (0.0, 0.0),
            (4.0, -2.0),
            (6.0, -5.0),
            (7.0, 11.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_0_1_energy_k0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (-8.0, 7.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_0_1_energy_l0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_0_1_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_0_1_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
        let pows = [2, 0, 0, 0, 1, 0];
        let coeffs = [
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_0_1_energy_k0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0];
        let pows = [2, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_0_1_energy_l0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_0_1_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (10.0, -3.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_0_2_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 0, 0, 2, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (10.0, -3.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_3_0_1_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [3, 0, 0, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (10.0, -3.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_3_0_2_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [3, 0, 0, 0, 2, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (10.0, -3.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_1_0_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 1, 0, 0, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 6.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (3.0, -23.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_1_0_const() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 1, 0, 0, 0];
        let coeffs = [
            (-2.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_1_0_energy_k0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 1, 0, 0, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_2_1_0_energy_l0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 1, 0, 0, 0];
        let coeffs = [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_3_3_0_dense() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 0, 1, 0, 0, 0];
        let coeffs = [
            (0.00803423874513387, 0.39917307594008244),
            (0.5941221615536947, 0.5927276753863615),
            (0.889457150507375, 0.7453389005854838),
            (0.9707167826454934, 0.935066923049863),
            (0.10555829555344776, 0.4337427188651112),
            (0.7648616295229604, 0.9687629427493974),
            (0.8847179070075125, 0.6190056436276886),
            (0.5042209277500286, 0.9718517512723205),
            (0.5724303614849603, 0.9368408776767316),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_1_1_energy_k0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 1, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_1_1_energy_l0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 1, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_1_1_1_energies() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [1, 0, 1, 0, 1, 0];
        let coeffs = [
            (0.0, 0.0),
            (1.0, 6.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (3.0, -23.0),
            (1.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn sunrise_linear_4_3_2_dense() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [4, 0, 3, 0, 2, 0];
        let coeffs = [
            (0.00803423874513387, 0.39917307594008244),
            (0.5941221615536947, 0.5927276753863615),
            (0.889457150507375, 0.7453389005854838),
            (0.9707167826454934, 0.935066923049863),
            (0.10555829555344776, 0.4337427188651112),
            (0.7648616295229604, 0.9687629427493974),
            (0.8847179070075125, 0.6190056436276886),
            (0.5042209277500286, 0.9718517512723205),
            (0.5724303614849603, 0.9368408776767316),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[test]
    fn double_triangle_linear_23_2_12_dense() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 3, 2, 0, 1, 2];
        let coeffs = [
            (0.00803423874513387, 0.39917307594008244),
            (0.5941221615536947, 0.5927276753863615),
            (0.889457150507375, 0.7453389005854838),
            (0.9707167826454934, 0.935066923049863),
            (0.10555829555344776, 0.4337427188651112),
            (0.7648616295229604, 0.9687629427493974),
            (0.8847179070075125, 0.6190056436276886),
            (0.5042209277500286, 0.9718517512723205),
            (0.5724303614849603, 0.9368408776767316),
        ];
        test_2L_2PROPxLL_linear(pows, &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn topo_2L_2xLL_rank5_21_12_11_dense() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        let pows = [2, 1, 1, 2, 1, 1];
        // Create the 1287 random coefficient needed for the rank 5 test at two loops
        let mut rng = rand::thread_rng();
        let coeffs: Vec<(f64, f64)> = (0..1287)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_2L_2PROPxLL_rank5(pows, &mut pt, &coeffs);
    }
}
mod dot_propagators_multiloops {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    /// MultiLoop.yaml contains different topologies for each loop number (1 to 6)
    /// that have different at least by the ltd_cut_structure
    ///
    /// |ALIASES
    /// |:--
    /// |1L_v1 -> Box_3E
    /// |1L_v2 -> manual_Decagon_P2_no_ellipse_massless_ala_weinzierl
    /// |2L_v1 -> DoubleBox_physical
    /// |2L_v2 -> DoubleTriangle_massive_physical
    /// |2L_v3 -> FISHNET_1x2
    /// |2L_v4 -> TriangleBox_physical
    /// |2L_v5 -> manual_6pt_Weinzierl_a
    /// |2L_v6 -> manual_AltDoubleTriangle
    /// |3L_v1 -> PRL_Mercedes
    /// |3L_v2 -> TriangleBoxTriangle
    /// |3L_v3 -> manual_Mercedes
    /// |3L_v4 -> manual_Mercedesv2
    /// |3L_v5 -> manual_TriangleBoxBox
    /// |3L_v6 -> manual_TriangleBoxBox_alt
    /// |4L_v1 -> FISHNET_1x4
    /// |4L_v2 -> FISHNET_2x2
    /// |4L_v3 -> PRL_BoxBoxBoxBox
    /// |4L_v4 -> PRL_non_planar_four_loop
    /// |4L_v5 -> manual_TriangleBoxBoxTriangle
    /// |5L_v1 -> FISHNET_1x5
    /// |6L_v1 -> FISHNET_1x6
    /// |6L_v2 -> FISHNET_2x3
    ///
    /// |#Loop | Version | #LoopLines | #propagators per loop line
    /// |:--:|:--:|:--:|:--
    /// |  1   |    1    |      1     |   4
    /// |      |    2    |      1     |   10
    /// |  2   |    1    |      3     |   3\|1\|3
    /// |      |    2    |      3     |   3\|1\|3
    /// |      |    3    |      3     |   3\|3\|1
    /// |      |    4    |      3     |   2\|1\|3
    /// |      |    5    |      3     |   2\|6\|1
    /// |      |    6    |      3     |   2\|2\|1
    /// |  3   |    1    |      6     |   2\|1\|1\|1\|1\|1
    /// |      |    2    |      5     |   2\|2\|2\|1\|1
    /// |      |    3    |      6     |   1\|1\|2\|1\|1\|1
    /// |      |    4    |      6     |   2\|1\|1\|1\|1\|1
    /// |      |    5    |      6     |   2\|1\|3\|1\|1\|1
    /// |      |    6    |      6     |   2\|1\|3\|1\|1\|1
    /// |  4   |    1    |      7     |   3\|2\|2\|3\|1\|1\|1
    /// |      |    2    |      8     |   2\|2\|1\|1\|1\|2\|1\|1
    /// |      |    3    |      7     |   3\|1\|2\|1\|2\|1\|3
    /// |      |    4    |      9     |   1\|1\|1\|1\|1\|1\|1\|1\|1
    /// |      |    5    |      9     |   2\|1\|1\|1\|2\|1\|1\|1\|1
    /// |  5   |    1    |      9     |   3\|2\|2\|2\|3\|1\|1\|1\|1
    /// |  6   |    1    |      11    |   3\|2\|2\|2\|2\|3\|1\|1\|1\|1\|1
    /// |      |    2    |      13    |   2\|1\|2\|1\|1\|1\|1\|1\|2\|1\|1\|2\|1
    #[allow(non_snake_case)]
    fn test_multiloops(name: &str, pt: &mut [f64], coeffs: &[(f64, f64)]) {
        let mut topo = get_test_topology("./tests/topologies/MultiLoops.yaml", name);
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
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);

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
    fn MultiLoop_1L_v1_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 1).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..126)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("1L_v1", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_1L_v2_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 1).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..126)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("1L_v2", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_2L_v1_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 2).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..1287)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("2L_v1", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_2L_v2_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 2).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..1287)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("2L_v2", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_2L_v3_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 2).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..1287)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("2L_v3", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_2L_v4_rank4_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 2).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 4 test
        //let coeffs: Vec<(f64, f64)> = (0..1287)
        let coeffs: Vec<(f64, f64)> = (0..495)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("2L_v4", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_2L_v5_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 2).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..1287)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("2L_v5", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_2L_v6_rank4_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 2).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 4 test
        let coeffs: Vec<(f64, f64)> = (0..495)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("2L_v6", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_3L_v1_rank4_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 3).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 4 test
        let coeffs: Vec<(f64, f64)> = (0..1820)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("3L_v1", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_3L_v2_rank4_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 3).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 4 test
        let coeffs: Vec<(f64, f64)> = (0..1820)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("3L_v2", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_3L_v3_linear_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 3).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..6188)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("3L_v3", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_3L_v4_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 3).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..6188)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("3L_v4", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_3L_v5_rank5_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 3).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 5 test
        let coeffs: Vec<(f64, f64)> = (0..6188)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("3L_v5", &mut pt, &coeffs);
    }

    #[allow(non_snake_case)]
    #[test]
    fn MultiLoop_3L_v6_rank4_dense() {
        let mut rng = rand::thread_rng();
        // Select point in momentum space: 3*n_loops
        let mut pt: Vec<f64> = (0..3 * 3).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
        // Create the random coefficients needed for the rank 4 test
        let coeffs: Vec<(f64, f64)> = (0..1820)
            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
            .collect();
        // Test a two loop topology with two propagators per loop line
        test_multiloops("3L_v6", &mut pt, &coeffs);
    }

    // Uncomment to test 4-6 loops topologies
    // WARNING: they need large RAM

    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_4L_v1_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 4).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..4845)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("4L_v1", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_4L_v2_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 4).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..4845)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("4L_v2", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_4L_v3_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 4).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..4845)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("4L_v3", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_4L_v4_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 4).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..4845)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("4L_v4", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_4L_v5_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 4).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..4845)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("4L_v5", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_5L_v1_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 5).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..10626)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("5L_v1", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_6L_v1_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 6).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..20475)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("6L_v1", &mut pt, &coeffs);
    //    }
    //
    //    #[allow(non_snake_case)]
    //    #[test]
    //    fn MultiLoop_6L_v2_rank4_dense() {
    //        let mut rng = rand::thread_rng();
    //        // Select point in momentum space: 3*n_loops
    //        let mut pt: Vec<f64> = (0..3 * 6).map(|_| 1.0 - 2.0 * rng.gen::<f64>()).collect();
    //        // Create the random coefficients needed for the rank 4 test
    //        let coeffs: Vec<(f64, f64)> = (0..20475)
    //            .map(|_| (1.0 - 2.0 * rng.gen::<f64>(), 1.0 - 2.0 * rng.gen::<f64>()))
    //            .collect();
    //        // Test a two loop topology with two propagators per loop line
    //        test_multiloops("6L_v2", &mut pt, &coeffs);
    //    }
}
