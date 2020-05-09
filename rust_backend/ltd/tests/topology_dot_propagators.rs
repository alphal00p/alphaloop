extern crate arrayvec;
extern crate f128;
extern crate ltd;
extern crate num;
extern crate num_traits;
extern crate vector;

use arrayvec::ArrayVec;
use ltd::topologies::{LTDCache, Topology};
use ltd::{DeformationStrategy, Settings, MAX_LOOP};
use num::Complex;
use num_traits::{Float, NumCast, One, ToPrimitive, Zero};
use vector::LorentzVector;

fn get_test_topology(topology_file: &str, topology_name: &str) -> Topology {
    // Import settings and disable the deformation
    let mut settings = Settings::from_file("../../LTD/hyperparameters.yaml");
    settings.general.deformation_strategy = DeformationStrategy::None;
    settings.general.multi_channeling = false;
    settings.integrator.dashboard = false;

    // Import topology from folder
    let mut topologies = Topology::from_file(topology_file, &settings);
    settings.general.topology = topology_name.to_string();
    topologies
        .remove(&settings.general.topology)
        .expect("Unknown topology")
}

fn cross_check<'a>(
    topo: &mut Topology,
    x: &'a [f64],
    cache: &mut LTDCache<f128::f128>,
) -> (
    Complex<f128::f128>,
    Complex<f128::f128>,
    Complex<f128::f128>,
) {
    // Parameterize
    let mut k = [LorentzVector::default(); MAX_LOOP];
    let mut jac_para = f128::f128::one();
    for i in 0..topo.n_loops {
        // set the loop index to i + 1 so that we can also shift k
        let (l_space, jac) =
            Topology::parameterize(&x[i * 3..(i + 1) * 3], topo.e_cm_squared, i, &topo.settings);

        // there could be some rounding here
        let rot = topo.rotation_matrix;
        k[i] = LorentzVector::from_args(
            f128::f128::zero(),
            <f128::f128 as NumCast>::from(rot[0][0]).unwrap() * l_space[0]
                + <f128::f128 as NumCast>::from(rot[0][1]).unwrap() * l_space[1]
                + <f128::f128 as NumCast>::from(rot[0][2]).unwrap() * l_space[2],
            <f128::f128 as NumCast>::from(rot[1][0]).unwrap() * l_space[0]
                + <f128::f128 as NumCast>::from(rot[1][1]).unwrap() * l_space[1]
                + <f128::f128 as NumCast>::from(rot[1][2]).unwrap() * l_space[2],
            <f128::f128 as NumCast>::from(rot[2][0]).unwrap() * l_space[0]
                + <f128::f128 as NumCast>::from(rot[2][1]).unwrap() * l_space[1]
                + <f128::f128 as NumCast>::from(rot[2][2]).unwrap() * l_space[2],
        );
        jac_para *= jac;
    }
    let mut k_def: ArrayVec<[LorentzVector<Complex<f128::f128>>; MAX_LOOP]> = (0..topo.n_loops)
        .map(|i| k[i].map(|x| Complex::new(x, f128::f128::zero())))
        .collect();
    topo.compute_complex_cut_energies(&k_def[..topo.n_loops], cache)
        .unwrap();

    // Compute using standard LTD
    topo.settings.general.partial_fractioning_threshold = -1.0;
    let result_ltd = topo
        .clone()
        .evaluate_all_dual_integrands(&mut k_def[..topo.n_loops], cache)
        .unwrap();

    // Compute using partial fractioning
    topo.settings.general.partial_fractioning_threshold = 1e-99;
    let result_pf1 = if topo.n_loops == 1 {
        topo.clone()
            .evaluate_all_dual_integrands(&mut k_def[..topo.n_loops], cache)
            .unwrap()
    } else {
        Complex::zero()
    };

    // Compute using partial fractioning MultiLoop
    for prop_pow in cache.propagator_powers.iter_mut() {
        *prop_pow = 0;
    }
    for ll in topo.loop_lines.iter() {
        for p in ll.propagators.iter() {
            cache.propagator_powers[p.id] = p.power;
        }
    }

    let result_pf2 = topo.partial_fractioning_multiloops.evaluate(
        &topo.numerator,
        &topo.loop_lines,
        &topo.propagator_id_to_ll_id,
        cache,
    );
    (result_ltd, result_pf1, result_pf2)
    //                         //println!("{:?}", self.partial_fractioning_multiloops.partial_fractioning_element);
    //                         for block in self
    //                             .partial_fractioning_multiloops
    //                             .partial_fractioning_element
    //                             .iter()
    //                         {
    //                             println!("factor: {:?}", block.factor);
    //                             print!("\t[");
    //                             for den in block.dens.iter() {
    //                                 let d1 = den.energies_and_shifts[0];
    //                                 let d2 = den.energies_and_shifts[1];
    //                                 let e1 = den.indices[0];
    //                                 let e2 = den.indices[1];
    //                                 if d1.1 < 0.0 {
    //                                     //println!("\t [{:?}, {:?}] @ [{},{}]", d1, d2, e1, e2);
    //                                     print!("({},{}), ", e1, e2);
    //                                 } else {
    //                                     print!("({},{}), ", e2, e1);
    //                                     //println!("\t [{:?}, {:?}] @ [{},{}]", d2, d1, e2, e1);
    //                                 }
    //                                 //println!("\t {:?} @ {:?}", &den.energies_and_shifts[..den.size], &den.indices[..den.size]);
    //                             }
    //                             println!("]");
    //                         }
    //                        //println!("{:?}", self.partial_fractioning.partial_fractioning_element);
    //                        for block in self.partial_fractioning.partial_fractioning_element.iter() {
    //                            if block.numerator.len() > 1 {
    //                                continue;
    //                            }
    //                            println!("numerator: {:?}", block.numerator);
    //                            println!("\t {:?}", block.ellipsoids_product);
    //                        }
    // if ((result.re - result_pf1.re) / result.re).powf(Into::<T>::into(2.0))
    //     > Into::<T>::into(1e-5)
    //     || ((result.im - result_pf1.im) / result.im).powf(Into::<T>::into(2.0))
    //         > Into::<T>::into(1e-5)
    // {
    //     println!("{} vs {}", result, result_pf1);
    //     println!("k = {:?}", k_def);
    // }
}

fn numeriacal_eq(x: Complex<f128::f128>, y: Complex<f128::f128>) -> bool {
    let real_check = if x.re == y.re {
        true
    } else {
        ((x.re - y.re) / x.re).abs().to_f64().unwrap() < 1e-29
    };
    let imag_check = if x.im == y.im {
        true
    } else {
        ((x.im - y.im) / x.im).abs().to_f64().unwrap() < 1e-29
    };

    println!("{:?} : {}", ((x.re - y.re) / x.re).abs(), real_check);
    println!("{:?} : {}", ((x.im - y.im) / x.im).abs(), imag_check);
    real_check && imag_check
}

/*    BEGIN TESTS    */
mod dot_propagators_1loop {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    fn test_bubble(pows: [usize; 2], pt: &mut [f64; 3], negative_cut: bool) {
        let mut topo = get_test_topology("./tests/topologies/bubble.yaml", "1l");
        topo.loop_lines[0].propagators[0].power = pows[0];
        topo.loop_lines[0].propagators[1].power = pows[1];
        if negative_cut {
            topo.ltd_cut_structure[0][0] = -1;
        } else {
            topo.ltd_cut_structure[0][0] = 1;
        }
        topo.process(true);

        // map the point back from momentum-space to the unit hypercube
        for (i, x) in pt.chunks_exact_mut(3).enumerate() {
            let r = Topology::inv_parametrize::<f128::f128>(
                &LorentzVector::from_args(0., x[0], x[1], x[2]).cast(),
                topo.e_cm_squared,
                i,
                &topo.settings,
            );
            x[0] = f128::f128::to_f64(&r.0[0]).unwrap();
            x[1] = f128::f128::to_f64(&r.0[1]).unwrap();
            x[2] = f128::f128::to_f64(&r.0[2]).unwrap();
        }
        // Compute values using e-fractioning the standard ltd sum over cuts
        let mut cache = LTDCache::<f128::f128>::new(&topo);
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);

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
    fn bubble_1_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1];

        test_bubble(pows, &mut pt, false);
    }

    #[test]
    fn bubble_1_2() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 2];

        test_bubble(pows, &mut pt, false);
    }

    #[test]
    fn bubble_2_2() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [2, 2];

        test_bubble(pows, &mut pt, false);
    }

    #[test]
    fn bubble_3_4() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [3, 4];

        test_bubble(pows, &mut pt, false);
    }

    #[test]
    fn bubble_2_0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [2, 0];

        test_bubble(pows, &mut pt, false);
    }

    #[test]
    fn bubble_3_0_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [3, 0];

        test_bubble(pows, &mut pt, true);
    }


    #[test]
    fn bubble_1_1_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [1, 1];

        test_bubble(pows, &mut pt, true);
    }

    #[test]
    fn bubble_2_1_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [2, 1];

        test_bubble(pows, &mut pt, true);
    }

    #[test]
    fn bubble_2_2_negative_cut() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [1., 2., 3.];
        let pows = [2, 2];

        test_bubble(pows, &mut pt, true);
    }
}

mod dot_propagators_2loops {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    fn test_2L_PROPxLL(pows: [usize; 6], pt: &mut [f64; 6]) {
        let mut topo = get_test_topology("./tests/topologies/2L_2PROPxLL.yaml", "2L_2PROPxLL");
        topo.loop_lines[0].propagators[0].power = pows[0];
        topo.loop_lines[0].propagators[1].power = pows[1];
        topo.loop_lines[1].propagators[0].power = pows[2];
        topo.loop_lines[1].propagators[1].power = pows[3];
        topo.loop_lines[2].propagators[0].power = pows[4];
        topo.loop_lines[2].propagators[1].power = pows[5];
        topo.process(true);

        // map the point back from momentum-space to the unit hypercube
        for (i, x) in pt.chunks_exact_mut(3).enumerate() {
            let r = Topology::inv_parametrize::<f128::f128>(
                &LorentzVector::from_args(0., x[0], x[1], x[2]).cast(),
                topo.e_cm_squared,
                i,
                &topo.settings,
            );
            x[0] = f128::f128::to_f64(&r.0[0]).unwrap();
            x[1] = f128::f128::to_f64(&r.0[1]).unwrap();
            x[2] = f128::f128::to_f64(&r.0[2]).unwrap();
        }
        // Compute values using e-fractioning the standard ltd sum over cuts
        let mut cache = LTDCache::<f128::f128>::new(&topo);
        let (result_ltd, result_pf1, result_pf2) = cross_check(&mut topo, pt, &mut cache);

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
    fn sunrise_1_1_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 1, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_1_0_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 0, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_3_0_5() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [3, 0, 0, 0, 5, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_1_1_0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 1, 0, 0, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_4_3_0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [4, 0, 3, 0, 0, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_2_1_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 0, 1, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_1_2_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 2, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_1_1_2() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 1, 0, 2, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_3_1_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [3, 0, 1, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }
    
    #[test]
    fn sunrise_1_3_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 3, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }


    #[test]
    fn sunrise_2_2_1() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 0, 2, 0, 1, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_2_1_2() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 0, 1, 0, 2, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_1_2_2() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 2, 0, 2, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn sunrise_2_2_2() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 0, 2, 0, 2, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }
    #[test]
    fn double_triangle_11_1_11() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 1, 1, 0, 1, 1];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn double_triangle_23_0_22() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 3, 0, 0, 2, 2];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn double_triangle_21_11_0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 1, 1, 1, 0, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn double_triangle_22_11_0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 2, 1, 1, 0, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn double_triangle_10_21_0() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [1, 0, 2, 1, 0, 0];

        test_2L_PROPxLL(pows, &mut pt);
    }

    #[test]
    fn double_triangle_21_2_32() {
        // Select point in momentum space: 3*n_loops
        let mut pt = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
        let pows = [2, 1, 2, 0, 3, 2];

        test_2L_PROPxLL(pows, &mut pt);
    }
}
