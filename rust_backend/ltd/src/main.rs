extern crate arrayvec;
extern crate clap;
extern crate cuba;
extern crate dual_num;
#[macro_use]
extern crate itertools;
extern crate colored;
extern crate f128;
extern crate ltd;
extern crate nalgebra as na;
extern crate num;
extern crate num_traits;
extern crate rand;
extern crate serde;
extern crate serde_yaml;
use serde::{Deserialize, Serialize};

use arrayvec::ArrayVec;
use clap::{App, Arg, ArgMatches, SubCommand};
use ltd::topologies::{Cut, CutList};
use ltd::LorentzVector;
use num_traits::real::Real;
use num_traits::ToPrimitive;
use num_traits::{NumCast, One, Zero};
use rand::prelude::*;
use std::str::FromStr;
use std::time::Instant;

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};

use cuba::{CubaIntegrator, CubaResult, CubaVerbosity};

use ltd::integrand::Integrand;
use ltd::topologies::{LTDCache, Surface, Topology};
use ltd::utils::Signum;
use ltd::{float, IntegratedPhase, Integrator, Settings};

use colored::*;

#[derive(Serialize, Deserialize)]
struct CubaResultDef {
    pub neval: i64,
    pub fail: i32,
    pub result: Vec<f64>,
    pub error: Vec<f64>,
    pub prob: Vec<f64>,
}

impl CubaResultDef {
    fn new(o: &CubaResult) -> CubaResultDef {
        CubaResultDef {
            neval: o.neval,
            fail: o.fail,
            result: o.result.clone(),
            error: o.error.clone(),
            prob: o.prob.clone(),
        }
    }
}

struct UserData {
    integrand: Vec<Integrand>,
    integrated_phase: IntegratedPhase,
}

#[inline(always)]
fn integrand(
    x: &[f64],
    f: &mut [f64],
    user_data: &mut UserData,
    _nvec: usize,
    core: i32,
) -> Result<(), &'static str> {
    let res = user_data.integrand[(core + 1) as usize].evaluate(x);

    if res.is_finite() {
        match user_data.integrated_phase {
            IntegratedPhase::Real => {
                f[0] = res.re.to_f64().unwrap();
            }
            IntegratedPhase::Imag => {
                f[0] = res.im.to_f64().unwrap();
            }
            IntegratedPhase::Both => {
                f[0] = res.re.to_f64().unwrap();
                f[1] = res.im.to_f64().unwrap();
            }
        }
    } else {
        f[0] = 0.;
    }

    Ok(())
}

fn bench(topo: &Topology, settings: &Settings) {
    let mut x = vec![0.; 3 * topo.n_loops];
    let mut rng = rand::thread_rng();

    let mut cache = LTDCache::<float>::new(&topo);

    let now = Instant::now();
    for _ in 0..settings.integrator.n_max {
        for xi in x.iter_mut() {
            *xi = rng.gen();
        }

        let _r = topo.evaluate(&x, &mut cache, &None);
    }

    println!("{:#?}", now.elapsed());
}

fn cb_to_lm(
    topo: &Topology,
    surf: &Surface,
    mat: &[i8],
    cut_momenta: &[LorentzVector<f128::f128>],
    loop_momenta: &mut [LorentzVector<f128::f128>],
) {
    // transform from cut momentum basis to loop momentum basis
    for (i, l) in loop_momenta.iter_mut().enumerate() {
        *l = LorentzVector::default();
        for (j, (c, e)) in mat[i * topo.n_loops..(i + 1) * topo.n_loops]
            .iter()
            .zip(&cut_momenta[..topo.n_loops])
            .enumerate()
        {
            *l += e.multiply_sign(*c);

            // subtract the shifts
            let mut index = 0;
            for (&cut, ll) in surf.cut.iter().zip(topo.loop_lines.iter()) {
                if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = cut {
                    if j == index {
                        *l -= ll.propagators[i].q.cast().multiply_sign(*c);
                    }
                    index += 1;
                }
            }
        }
    }
}

fn point_generator<'a>(
    topo: &Topology,
    surf: &Surface,
    rescaling: f64,
    mat: &[i8],
    pos_mom: &mut [LorentzVector<f128::f128>],
    neg_mom: &mut [LorentzVector<f128::f128>],
) -> bool {
    let mut rng = rand::thread_rng();

    let mut neg_surface_signs_count = surf.signs.iter().filter(|x| **x == -1).count();
    let mut pos_surface_signs_count = surf.signs.iter().filter(|x| **x == 1).count();
    if surf.delta_sign > 0 {
        pos_surface_signs_count += 1;
    } else {
        neg_surface_signs_count += 1;
    }

    let mut loop_momenta: ArrayVec<[LorentzVector<f128::f128>; ltd::MAX_LOOP]> = (0..topo.n_loops)
        .map(|_| LorentzVector::default())
        .collect();

    let mut cut_momenta: ArrayVec<[LorentzVector<f128::f128>; ltd::MAX_LOOP]> = (0..topo.n_loops)
        .map(|_| LorentzVector::default())
        .collect();

    let cut_masses: ArrayVec<[f128::f128; ltd::MAX_LOOP]> = surf
        .cut
        .iter()
        .zip(topo.loop_lines.iter())
        .filter_map(|(cut, ll)| {
            if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = cut {
                Some(ll.propagators[*i].m_squared.into())
            } else {
                None
            }
        })
        .collect();

    let surf_mass: f128::f128 = topo.loop_lines[surf.onshell_ll_index].propagators
        [surf.onshell_prop_index]
        .m_squared
        .into();

    // sample a random point
    for cm in cut_momenta.iter_mut() {
        for index in 1..4 {
            cm[index] =
                ((rng.gen::<f64>() * 2.0 - 1.0) * topo.e_cm_squared.sqrt() * rescaling).into();
        }
    }

    if pos_surface_signs_count == 1 && neg_surface_signs_count == 1 {
        // find the relevant cut
        let index = surf.signs.iter().position(|x| *x != 0).unwrap();

        // the one-loop case can be solved analytically
        let mut k = cut_momenta[index].spatial_squared().sqrt();
        let p = surf.shift.spatial_squared().sqrt();

        let mut costheta = f128::f128::zero();
        while k < f128::f128::INFINITY {
            costheta = (cut_masses[index] * cut_masses[index] - surf_mass * surf_mass - p * p
                + Into::<f128::f128>::into(2.)
                    * (k * k + cut_masses[index]).sqrt()
                    * surf.shift.t.multiply_sign(-surf.delta_sign)
                + surf.shift.t * surf.shift.t)
                / (Into::<f128::f128>::into(2.) * k * p);
            if costheta >= -f128::f128::one() && costheta <= f128::f128::one() {
                break;
            }

            k *= Into::<f128::f128>::into(2.0);
        }

        let pv: LorentzVector<f128::f128> = surf.shift.cast();
        let perp = if pv.z.is_zero() && (pv.x - pv.y).is_zero() {
            LorentzVector::from_args(f128::f128::zero(), pv.y - pv.z, pv.x, pv.x)
        } else {
            LorentzVector::from_args(f128::f128::zero(), pv.z, pv.z, -pv.x - pv.y)
        };
        let k_perp = (k * k - (k * k * costheta * costheta)).sqrt() / perp.spatial_squared().sqrt();
        cut_momenta[index] = pv * (k * costheta / p) + perp * k_perp;
    }

    // transform from cut momentum basis to loop momentum basis
    cb_to_lm(topo, surf, mat, &cut_momenta, &mut loop_momenta);
    let res = evaluate_surface(topo, surf, &loop_momenta);

    if res.abs() < Into::<f128::f128>::into(1e-15) {
        for i in 0..topo.n_loops {
            pos_mom[i] = loop_momenta[i];
            neg_mom[i] = loop_momenta[i];
        }

        return true;
    }

    if pos_surface_signs_count == 1 && neg_surface_signs_count == 1 {
        println!(
            "{} {}",
            "One-loop hyperboloid not correctly sampled:".red(),
            res
        );
        return false;
    }

    let need_positive = res < f128::f128::zero();

    if need_positive {
        for i in 0..topo.n_loops {
            neg_mom[i] = loop_momenta[i];
        }
    } else {
        for i in 0..topo.n_loops {
            pos_mom[i] = loop_momenta[i];
        }
    }

    let mut branch; // for debugging
    if need_positive && pos_surface_signs_count > 1 || !need_positive && neg_surface_signs_count > 1
    {
        let mut indices = [100; 2];
        let mut index = 0;

        for (i, &s) in surf.signs.iter().enumerate() {
            if s > 0 && need_positive || s < 0 && !need_positive {
                indices[index] = i;
                index += 1;

                if index == 2 {
                    break;
                }
            }
        }

        if index == 1 {
            // the second index would be the surface term!
            // in this case we simply scale
            let q1 = (cut_momenta[indices[0]].spatial_squared() + cut_masses[indices[0]]).sqrt();

            // evaluate the contributions without the index we want and without the surface term
            let mut res1 = surf.shift.t;
            for (i, (s, cm, mass)) in izip!(&surf.signs, &cut_momenta, &cut_masses).enumerate() {
                if i != indices[0] {
                    res1 += (cm.spatial_squared() + mass).sqrt().multiply_sign(*s);
                }
            }

            cut_momenta[indices[0]] *= res1 / q1;
            branch = 1;
        } else {
            // let q1 and q2 both have a + or -
            // q1 -> q1 - s * q2 *lambda
            // q2 -> q2 + q2 * lambda
            // then the negative surface term stays the same
            // s is the relative surface sign of q1 and q2
            let s = surf.sig_ll_in_cb[indices[0]] * surf.sig_ll_in_cb[indices[1]];

            let q1 = (cut_momenta[indices[0]].spatial_squared() + cut_masses[indices[0]]).sqrt();
            let q2 = (cut_momenta[indices[1]].spatial_squared() + cut_masses[indices[1]]).sqrt();
            let lambda = (-res
                + q1.multiply_sign(surf.signs[indices[0]])
                + q2.multiply_sign(surf.signs[indices[1]]))
                / q2
                - f128::f128::one();

            cut_momenta[indices[0]] -= cut_momenta[indices[1]] * lambda * s.into();
            cut_momenta[indices[1]] += cut_momenta[indices[1]] * lambda;
            branch = 2;
        }
    } else {
        // we are in the case where we only have 1 term (or 0) with the sign we need
        // set it equal to minus the shift and set the rest to their mass
        // this is the minimal solution

        branch = 3; // branch 3: sign is on the surface term
        for (cm, &s, &mom_sign) in izip!(
            cut_momenta.iter_mut(),
            surf.signs.iter(),
            surf.sig_ll_in_cb.iter()
        ) {
            if s > 0 && need_positive || s < 0 && !need_positive {
                *cm = -surf.shift.cast().multiply_sign(mom_sign);
                branch = 4;
            } else {
                *cm = LorentzVector::default();
            }
        }
    }

    cb_to_lm(topo, surf, mat, &cut_momenta, &mut loop_momenta);
    let res2 = evaluate_surface(topo, surf, &loop_momenta);

    if res.signum() == res2.signum() {
        println!(
            "{} {} vs {}, branch: {}",
            "FAILED to get a different sign:".red(),
            res,
            res2,
            branch
        );

        for (i, x) in cut_momenta.iter().enumerate() {
            println!(
                "q{}={}; |q{}| = {}",
                i + 1,
                x,
                i + 1,
                x.spatial_squared().sqrt()
            );
        }
    }

    if need_positive {
        for i in 0..topo.n_loops {
            pos_mom[i] = loop_momenta[i];
        }
    } else {
        for i in 0..topo.n_loops {
            neg_mom[i] = loop_momenta[i];
        }
    }

    res.signum() != res2.signum()
}

fn evaluate_surface(
    topo: &Topology,
    surf: &Surface,
    loop_momenta: &[LorentzVector<f128::f128>],
) -> f128::f128 {
    let mut res = f128::f128::zero();

    let mut cut_index = 0;
    for (cut, ll) in izip!(surf.cut.iter(), topo.loop_lines.iter()) {
        if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = *cut {
            if surf.signs[cut_index] == 0 {
                cut_index += 1;
                continue;
            }

            // construct the cut energy
            let mut mom = LorentzVector::<f128::f128>::default();
            for (&cut_sign, lm) in ll.signature.iter().zip(loop_momenta.iter()) {
                mom += lm * cut_sign.into();
            }
            // compute the postive cut energy
            let q: LorentzVector<f128::f128> = ll.propagators[i].q.cast();
            let energy = ((mom + q).spatial_squared()
                + <f128::f128 as NumCast>::from(ll.propagators[i].m_squared).unwrap())
            .sqrt();

            res += energy.multiply_sign(surf.signs[cut_index]);

            cut_index += 1;
        }
    }

    // now for the surface term
    let mut mom = LorentzVector::<f128::f128>::default();
    let onshell_ll = &topo.loop_lines[surf.onshell_ll_index];
    let onshell_prop = &onshell_ll.propagators[surf.onshell_prop_index];
    for (&surf_sign, lm) in onshell_ll.signature.iter().zip(loop_momenta.iter()) {
        mom += lm * surf_sign.into();
    }
    let energy = ((mom + onshell_prop.q.cast()).spatial_squared()
        + <f128::f128 as NumCast>::from(onshell_prop.m_squared).unwrap())
    .sqrt();

    res += energy.multiply_sign(surf.delta_sign);
    res += <f128::f128 as NumCast>::from(surf.shift.t).unwrap();
    res
}

// TODO: move to diagnostics
fn surface_prober<'a>(topo: &Topology, settings: &Settings, matches: &ArgMatches<'a>) {
    let mut loop_momenta = vec![LorentzVector::<f128::f128>::default(); topo.n_loops];

    let mut k_def: ArrayVec<[LorentzVector<num::Complex<f128::f128>>; ltd::MAX_LOOP]>;
    let mut cache = LTDCache::new(topo);

    let mut positive_lm = vec![LorentzVector::<f128::f128>::default(); topo.n_loops];
    let mut negative_lm = vec![LorentzVector::<f128::f128>::default(); topo.n_loops];

    let ids: Vec<_> = match matches.values_of("ids") {
        Some(x) => x.map(|x| usize::from_str(x).unwrap()).collect(),
        None => vec![],
    };

    let samples = usize::from_str(matches.value_of("samples").unwrap()).unwrap();
    let rescaling = f64::from_str(matches.value_of("rescaling").unwrap()).unwrap();

    for (surf_index, surf) in topo.surfaces.iter().enumerate() {
        if !ids.is_empty() && !ids.contains(&surf_index) {
            continue;
        }

        println!(
            "-> group={}, ellipsoid={}, prop={:?} cut={}, mom_map={:?}, signs={:?}, marker={}, shift={}",
            surf.group, surf.ellipsoid, (surf.onshell_ll_index, surf.onshell_prop_index), CutList(&surf.cut), surf.sig_ll_in_cb,
            surf.signs, surf.delta_sign, surf.shift
        );

        let onshell_ll = &topo.loop_lines[surf.onshell_ll_index];
        let onshell_prop = &onshell_ll.propagators[surf.onshell_prop_index];

        for _ in 0..samples {
            let mut did_break = false;

            if point_generator(
                topo,
                surf,
                rescaling,
                &topo.cb_to_lmb_mat[surf.cut_structure_index],
                &mut positive_lm,
                &mut negative_lm,
            ) {
                // try to bisect
                for _ in 0..1000 {
                    for (lm, pl, nl) in izip!(
                        loop_momenta.iter_mut(),
                        positive_lm.iter(),
                        negative_lm.iter()
                    ) {
                        *lm = (pl + nl) * Into::<f128::f128>::into(0.5);
                    }

                    let res = evaluate_surface(topo, surf, &loop_momenta);

                    // update the bounds
                    if res < f128::f128::zero() {
                        for (nl, ll) in negative_lm.iter_mut().zip(loop_momenta.iter()) {
                            *nl = ll.clone();
                        }
                    }
                    if res > f128::f128::zero() {
                        for (pl, ll) in positive_lm.iter_mut().zip(loop_momenta.iter()) {
                            *pl = ll.clone();
                        }
                    }

                    if res.abs() < 1e-14.into() {
                        if settings.general.debug > 0 {
                            println!("Found point {}", res);
                        }

                        // check the pole for ellipsoids
                        if surf.ellipsoid {
                            // set the loop momenta
                            let (kappas, _) = topo.deform(&loop_momenta, None, &mut cache);
                            k_def = (0..topo.n_loops)
                                .map(|i| {
                                    loop_momenta[i]
                                        .map(|x| num::Complex::new(x, f128::f128::zero()))
                                        + kappas[i]
                                            .map(|x| num::Complex::new(f128::f128::zero(), x))
                                })
                                .collect();

                            // do a full evaluation
                            let mut result = num::Complex::<f128::f128>::default();
                            if topo
                                .compute_complex_cut_energies(&k_def, &mut cache)
                                .is_ok()
                            {
                                for (cuts, mat) in
                                    topo.ltd_cut_options.iter().zip(topo.cb_to_lmb_mat.iter())
                                {
                                    for cut in cuts.iter() {
                                        result += topo
                                            .evaluate_cut(&mut k_def, cut, mat, &mut cache)
                                            .unwrap();

                                        // check the pole of the on-shell propagator for ellipsoids
                                        if surf.ellipsoid && *cut == surf.cut {
                                            let mut e: num::Complex<f128::f128> =
                                                num::Complex::default();
                                            for (l, &c) in
                                                k_def.iter().zip(onshell_ll.signature.iter())
                                            {
                                                e += l.t.multiply_sign(c);
                                            }

                                            let r = ltd::utils::powi(
                                                e + Into::<f128::f128>::into(onshell_prop.q.t),
                                                2,
                                            ) - cache.complex_prop_spatial[onshell_prop.id];

                                            if r.im <= f128::f128::zero() {
                                                println!(
                                                    "{} for cut={}, sig={:?}, prop_q={}: {}",
                                                    "Bad pole detected".red(),
                                                    CutList(cut),
                                                    onshell_ll.signature,
                                                    onshell_prop.q,
                                                    r.im
                                                );
                                            }
                                        }
                                    }
                                }
                            }

                            //println!("Result: {} at distance {}", result, res);
                            if surf.ellipsoid && result.norm() > Into::<f128::f128>::into(1e8) {
                                println!(
                                    "{} {:e} at distance {} for {:?}",
                                    "Large result".red(),
                                    result,
                                    res,
                                    loop_momenta
                                );
                            }
                        } else {
                            // check the dual cancelations by probing points close to the dual canceling surface
                            let mut probes = [
                                num::Complex::<f128::f128>::default(),
                                num::Complex::<f128::f128>::default(),
                                num::Complex::<f128::f128>::default(),
                            ];
                            if !surf.ellipsoid {
                                for (probe, lambda) in probes.iter_mut().zip(&[
                                    1.000000001,
                                    1.0000000089999991,
                                    1.00000008999991,
                                ]) {
                                    if settings.general.debug > 4 {
                                        println!("Testing lambda {}", lambda);
                                    }

                                    for lm in &mut loop_momenta {
                                        *lm *= Into::<f128::f128>::into(*lambda);
                                    }

                                    // set the loop momenta
                                    let (kappas, _) = topo.deform(&loop_momenta, None, &mut cache);
                                    k_def = (0..topo.n_loops)
                                        .map(|i| {
                                            loop_momenta[i]
                                                .map(|x| num::Complex::new(x, f128::f128::zero()))
                                                + kappas[i].map(|x| {
                                                    num::Complex::new(f128::f128::zero(), x)
                                                })
                                        })
                                        .collect();

                                    // do a full evaluation
                                    if topo
                                        .compute_complex_cut_energies(&k_def, &mut cache)
                                        .is_ok()
                                    {
                                        for (cuts, mat) in topo
                                            .ltd_cut_options
                                            .iter()
                                            .zip(topo.cb_to_lmb_mat.iter())
                                        {
                                            for cut in cuts.iter() {
                                                *probe += topo
                                                    .evaluate_cut(&mut k_def, cut, mat, &mut cache)
                                                    .unwrap();
                                            }
                                        }
                                    }
                                }

                                let mut a: Vec<_> = probes.iter().map(|x| x.norm()).collect();
                                a.sort_by(|a, b| a.partial_cmp(b).unwrap());

                                let pv: f128::f128 =
                                    (a.last().unwrap() - a.first().unwrap()) / (a[a.len() / 2]);

                                if pv > Into::<f128::f128>::into(1e-3) {
                                    println!(
                                        "{}: pv={:e}, probes={:?}",
                                        "Dual cancellation breakdown detected".red(),
                                        pv,
                                        probes
                                    );
                                }
                            }
                        }
                        did_break = true;
                        break;
                    }
                }

                if !did_break {
                    println!(
                        "Could not bisect for surface with cut={} and os={:?}",
                        CutList(&surf.cut),
                        (surf.onshell_ll_index, surf.onshell_prop_index)
                    );
                }
            }
        }
    }
}

fn inspect<'a>(topo: &Topology, settings: &mut Settings, matches: &ArgMatches<'a>) {
    let pt: Vec<_> = matches
        .values_of("point")
        .unwrap()
        .map(|x| f64::from_str(x).unwrap())
        .collect();
    if pt.len() != 3 * topo.n_loops {
        panic!(
            "Dimension of the input point is incorrect. It should be {} but is {}.",
            topo.n_loops * 3,
            pt.len()
        );
    }

    if matches.is_present("full_integrand") {
        settings.general.screen_log_core = Some(1);
        settings.general.log_points_to_screen = true;
        let mut i = Integrand::new(&topo, settings.clone(), 1);
        i.evaluate(&pt);
        return;
    }

    // TODO: prevent code repetition
    if matches.is_present("use_f128") {
        let mut cache = LTDCache::<f128::f128>::new(&topo);
        let (x, k_def, jac_para, jac_def, result) = topo.clone().evaluate(&pt, &mut cache, &None);
        match topo.n_loops {
            1 => {
                println!(
                    "result={:e}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                    result, x, k_def[0], jac_para, jac_def
                );
            }
            2 => {
                println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result, x, k_def[0], k_def[1], jac_para, jac_def
                        );
            }
            3 => {
                println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    );
            }
            4 => {
                println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}n  | n={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], k_def[3], jac_para, jac_def
                    );
            }
            _ => {}
        }
    } else {
        let mut cache = LTDCache::<float>::new(&topo);
        let (x, k_def, jac_para, jac_def, result) = topo.clone().evaluate(&pt, &mut cache, &None);
        match topo.n_loops {
            1 => {
                println!(
                    "result={:e}\n  | x={:?}\n  | k={:e}\n  | jac_para={:e}, jac_def={:e}",
                    result, x, k_def[0], jac_para, jac_def
                );
            }
            2 => {
                println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result, x, k_def[0], k_def[1], jac_para, jac_def
                        );
            }
            3 => {
                println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], jac_para, jac_def
                    );
            }
            4 => {
                println!(
                        "result={:e}\n  | x={:?}\n  | k={:e}\n  | l={:e}\n  | m={:e}n  | n={:e}\n  | jac_para={:e}, jac_def={:e}",
                        result,x, k_def[0], k_def[1], k_def[2], k_def[3], jac_para, jac_def
                    );
            }
            _ => {}
        }
    }
}

fn main() {
    let matches = App::new("Feynman diagram integrator")
        .version("0.1")
        .about("Numerically integrate your favourite integrals with LTD")
        .arg(
            Arg::with_name("cores")
                .short("c")
                .long("cores")
                .value_name("NUMCORES")
                .help("Set the number of cores"),
        )
        .arg(
            Arg::with_name("samples")
                .short("s")
                .long("samples")
                .value_name("SAMPLES")
                .help("Number of samples per integration"),
        )
        .arg(
            Arg::with_name("topologies")
                .short("l")
                .long("topologies")
                .value_name("TOPOLOGY_FILE")
                .default_value("../LTD/topologies.yaml")
                .help("Set the topology file"),
        )
        .arg(
            Arg::with_name("config")
                .short("f")
                .long("config")
                .value_name("CONFIG_FILE")
                .default_value("../LTD/hyperparameters.yaml")
                .help("Set the configuration file"),
        )
        .arg(
            Arg::with_name("deformation")
                .short("d")
                .long("deformation")
                .value_name("DEFORMATION")
                .default_value("none")
                .help("Set the deformation"),
        )
        .arg(
            Arg::with_name("topology")
                .short("t")
                .long("topology")
                .value_name("TOPOLOGY")
                .help("Set the active topology"),
        )
        .subcommand(SubCommand::with_name("bench").about("Run a benchmark"))
        .subcommand(
            SubCommand::with_name("probe")
                .about("Sample points on hyperboloids and ellipsoids")
                .arg(
                    Arg::with_name("ids")
                        .long("ids")
                        .min_values(1)
                        .help("Only sample these surface ids"),
                )
                .arg(
                    Arg::with_name("samples")
                        .short("s")
                        .default_value("100")
                        .help("Number of samples per surface"),
                )
                .arg(
                    Arg::with_name("rescaling")
                        .short("r")
                        .default_value("1.")
                        .help("Rescale the sampling range by this factor"),
                ),
        )
        .subcommand(
            SubCommand::with_name("inspect")
                .about("Inspect a single input point")
                .arg(Arg::with_name("point").required(true).min_values(3))
                .arg(
                    Arg::with_name("use_f128")
                        .short("f128")
                        .long("use_f128")
                        .help("Use f128 evaluation"),
                )
                .arg(
                    Arg::with_name("full_integrand")
                        .long("full_integrand")
                        .help("Evaluate the integrand and possibly its rotated vesion"),
                ),
        )
        .get_matches();

    let mut settings = Settings::from_file(matches.value_of("config").unwrap());

    let mut cores = 1;
    if let Some(x) = matches.value_of("cores") {
        cores = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("samples") {
        settings.integrator.n_max = usize::from_str(x).unwrap();
    }

    if let Some(x) = matches.value_of("topology") {
        settings.general.topology = x.to_owned();
    }

    let topology_file = matches.value_of("topologies").unwrap();

    if let Some(x) = matches.value_of("deformation_strategy") {
        settings.general.deformation_strategy = x.into();
    }
    let mut ci = CubaIntegrator::new(integrand);

    ci.set_mineval(10)
        .set_nstart(settings.integrator.n_start as i64)
        .set_nincrease(settings.integrator.n_increase as i64)
        .set_maxeval(settings.integrator.n_max as i64)
        .set_epsabs(0.)
        .set_epsrel(1e-15)
        .set_seed(settings.integrator.seed)
        .set_cores(cores, 1000);

    // load the example file
    let mut topologies = Topology::from_file(topology_file, &settings);
    let topo = topologies
        .get_mut(&settings.general.topology)
        .expect("Unknown topology");
    topo.process();

    if let Some(_) = matches.subcommand_matches("bench") {
        bench(&topo, &settings);
        return;
    }

    if let Some(matches) = matches.subcommand_matches("probe") {
        surface_prober(&topo, &settings, matches);
        return;
    }

    if let Some(matches) = matches.subcommand_matches("inspect") {
        inspect(&topo, &mut settings, matches);
        return;
    }

    println!(
        "Integrating {} with {} samples and deformation '{}'",
        settings.general.topology, settings.integrator.n_max, settings.general.deformation_strategy
    );

    match topo.analytical_result_real {
        Some(_) => println!(
            "Analytic result: {:e}",
            num::Complex::<f64>::new(
                topo.analytical_result_real.unwrap(),
                topo.analytical_result_imag.unwrap()
            )
        ),
        _ => println!("Analytic result not available."),
    }

    let cuba_result = match settings.integrator.integrator {
        Integrator::Vegas => ci.vegas(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            CubaVerbosity::Progress,
            0,
            UserData {
                integrand: (0..=cores)
                    .map(|i| Integrand::new(topo, settings.clone(), i))
                    .collect(),
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
        Integrator::Suave => ci.suave(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            settings.integrator.n_new,
            settings.integrator.n_min,
            settings.integrator.flatness,
            CubaVerbosity::Progress,
            UserData {
                integrand: (0..=cores)
                    .map(|i| Integrand::new(topo, settings.clone(), i))
                    .collect(),
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
        Integrator::Cuhre => ci.cuhre(
            3 * topo.n_loops,
            if settings.integrator.integrated_phase == IntegratedPhase::Both {
                2
            } else {
                1
            },
            CubaVerbosity::Progress,
            UserData {
                integrand: (0..=cores)
                    .map(|i| Integrand::new(topo, settings.clone(), i))
                    .collect(),
                integrated_phase: settings.integrator.integrated_phase,
            },
        ),
    };
    println!("{:#?}", cuba_result);
    match topo.analytical_result_real {
        Some(_) => println!(
            "Analytic result: {:e}",
            num::Complex::<f64>::new(
                topo.analytical_result_real.unwrap(),
                topo.analytical_result_imag.unwrap()
            )
        ),
        _ => println!("Analytic result not available."),
    }
    let f = OpenOptions::new()
        .create(true)
        .write(true)
        .open(settings.general.topology.clone() + "_res.dat")
        .expect("Unable to create result file");
    let mut result_file = BufWriter::new(f);

    // write the result to a file
    writeln!(
        &mut result_file,
        "{}",
        serde_yaml::to_string(&CubaResultDef::new(&cuba_result)).unwrap()
    )
    .unwrap();
    writeln!(&mut result_file, "...").unwrap(); // write end-marker, for easy streaming
}
