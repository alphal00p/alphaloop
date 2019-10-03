use arrayvec::ArrayVec;
use colored::Colorize;
use disjoint_sets::UnionFind;
use dual_num::{DimName, DualN};
use itertools::Itertools;
use num::Complex;
use num_traits::ops::inv::Inv;
use num_traits::{Float, FloatConst, FromPrimitive, NumCast, One, Signed, Zero};
use topologies::{CacheSelector, Cut, CutList, LTDCache, LoopLine, Surface, SurfaceType, Topology};
use utils::Signum;
use vector::LorentzVector;
use {float, PythonNumerator};
use {
    AdditiveMode, DeformationStrategy, FloatLike, OverallDeformationScaling,
    ParameterizationMapping, ParameterizationMode, MAX_LOOP,
};

use utils;

type Dual4<T> = DualN<T, dual_num::U4>;
type Dual7<T> = DualN<T, dual_num::U7>;
type Dual10<T> = DualN<T, dual_num::U10>;
type Dual13<T> = DualN<T, dual_num::U13>;
type Dual16<T> = DualN<T, dual_num::U16>;
type Dual19<T> = DualN<T, dual_num::U19>;

impl LoopLine {
    /// Return the inverse of the evaluated loop line
    fn evaluate<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<num::Complex<T>>],
        cut: &Cut,
        topo: &Topology,
        cache: &LTDCache<T>,
    ) -> Result<num::Complex<T>, &'static str> {
        // only set the values when we use them
        let (kinematics_scale, threshold) = if topo.settings.general.numerical_threshold > 0. {
            (
                Into::<T>::into(topo.e_cm_squared),
                Into::<T>::into(topo.settings.general.numerical_threshold),
            )
        } else {
            (T::zero(), T::zero())
        };

        // construct the loop energy
        // if the line is cut, we can get it from the cache
        let e = match cut {
            Cut::PositiveCut(i) => {
                cache.complex_cut_energies[self.propagators[*i].id]
                    - Into::<T>::into(self.propagators[*i].q.t)
            }
            Cut::NegativeCut(i) => {
                -cache.complex_cut_energies[self.propagators[*i].id]
                    - Into::<T>::into(self.propagators[*i].q.t)
            }
            Cut::NoCut => {
                let mut e: num::Complex<T> = Complex::default();
                for (l, &c) in loop_momenta.iter().zip_eq(self.signature.iter()) {
                    e += l.t * T::from_i8(c).unwrap();
                }
                e
            }
        };

        let mut res = num::Complex::one();
        if topo.settings.general.debug > 3 {
            println!(
                "Loop line evaluation for cut {}\n  | signature = {:?}",
                cut, self.signature
            );
        }
        for (i, p) in self.propagators.iter().enumerate() {
            match cut {
                Cut::PositiveCut(j) | Cut::NegativeCut(j) if i == *j => {
                    let r = cache.complex_cut_energies[p.id] * Into::<T>::into(2.);
                    res *= r;

                    if topo.settings.general.debug > 3 {
                        println!("  | prop x {}={}", i, r);
                    }
                }
                _ => {
                    // multiply dual propagator
                    let r = utils::powi(e + Into::<T>::into(p.q.t), 2)
                        - cache.complex_prop_spatial[p.id];

                    if topo.settings.general.debug > 3 {
                        println!("  | prop   {}={}", i, r);
                    }

                    if !r.is_finite()
                        || (topo.settings.general.numerical_threshold > 0.
                            && r.re * r.re < kinematics_scale * threshold)
                    {
                        return Err("numerical instability");
                    }

                    res *= r;
                }
            }
        }
        Ok(res)
    }
}

impl Topology {
    /// Compute LTD related quantities for this topology.
    pub fn process(&mut self) {
        if self.n_loops > MAX_LOOP {
            panic!("Please increase MAX_LOOP to {}", self.n_loops);
        }

        assert!(
            self.settings.general.partial_fractioning
                || self.settings.general.deformation_strategy != DeformationStrategy::Intersections
        );

        // set the identity rotation matrix
        self.rotation_matrix = [
            [float::one(), float::zero(), float::zero()],
            [float::zero(), float::one(), float::zero()],
            [float::zero(), float::zero(), float::one()],
        ];

        // copy the signature to the propagators
        let mut prop_id = 0;
        for l in &mut self.loop_lines {
            for p in &mut l.propagators {
                p.signature = l.signature.clone();
                p.id = prop_id;
                prop_id += 1;
            }
        }

        // determine the e_cm_squared and on-shell flag
        let mut v: LorentzVector<f64> = LorentzVector::default();
        for e in &self.external_kinematics {
            if e.t > 0. {
                v += *e;
            }
        }

        self.e_cm_squared = v.square_impr().abs();
        if self.e_cm_squared == 0. {
            eprintln!("e_cm is zero: taking the abs of the spatial part instead");
            self.e_cm_squared = self.external_kinematics[0].spatial_squared_impr();
        }

        // determine the on-shell flag
        for (i, e) in self.external_kinematics.iter().enumerate() {
            if e.square_impr().abs() < 1e-10 * self.e_cm_squared {
                self.on_shell_flag |= 2_usize.pow(i as u32);
            }
        }

        // compute the cartesian product of cut structures
        self.ltd_cut_options = self
            .ltd_cut_structure
            .iter()
            .map(|cs| {
                cs.iter()
                    .zip_eq(self.loop_lines.iter())
                    .map(|(ll_cut, ll)| match ll_cut {
                        0 => vec![Cut::NoCut],
                        -1 => (0..ll.propagators.len())
                            .map(|i| Cut::NegativeCut(i))
                            .collect::<Vec<_>>(),
                        1 => (0..ll.propagators.len())
                            .map(|i| Cut::PositiveCut(i))
                            .collect::<Vec<_>>(),
                        _ => unreachable!("Cut structure can only be -1,0,1"),
                    })
                    .multi_cartesian_product()
                    .collect::<Vec<_>>()
            })
            .collect();

        // solve the linear system for fixing the energy component of each loop line
        for (cut_index, (residue_sign, cut_options)) in self
            .ltd_cut_structure
            .iter()
            .zip_eq(self.ltd_cut_options.iter())
            .enumerate()
        {
            let mut cut_signatures_matrix = vec![];
            let mut cut_residue_sign = vec![];

            // note that negative cuts will be treated later
            for (c, ll) in residue_sign.iter().zip_eq(self.loop_lines.iter()) {
                if *c != 0 {
                    for mom in &ll.signature {
                        cut_signatures_matrix.push(*mom as f64);
                    }
                    cut_residue_sign.push(*c);
                }
            }

            let b = na::DMatrix::from_row_slice(self.n_loops, self.n_loops, &cut_signatures_matrix);
            let c = b.try_inverse().unwrap();
            let c_i8 = c.map(|x| x as i8);
            // note: this transpose is ONLY there because the iteration is column-wise instead of row-wise
            let mat = c.transpose().iter().map(|x| *x as i8).collect(); // TODO: check if safe
            self.cb_to_lmb_mat.push(mat);

            // find all ellipsoids per cut option by expressing each propagator in terms of
            // the cut momenta
            for (cut_option_index, cut_option) in cut_options.iter().enumerate() {
                let mut cut_shift: Vec<LorentzVector<float>> = vec![]; // qs of cut
                let mut cut_mass = vec![];
                for (cut, ll) in cut_option.iter().zip_eq(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut
                    {
                        cut_shift.push(ll.propagators[*cut_prop_index].q.cast());
                        cut_mass.push(
                            float::from_f64(ll.propagators[*cut_prop_index].m_squared).unwrap(),
                        );
                    }
                }

                for (ll_index, ll) in self.loop_lines.iter().enumerate() {
                    // map the signature from the loop line to the cut momentum basis
                    let sig_ll_in_cb = c_i8.transpose()
                        * na::DMatrix::from_row_slice(self.n_loops, 1, &ll.signature);

                    for (onshell_prop_index, onshell_prop) in ll.propagators.iter().enumerate() {
                        // do not consider the cut propagator to have an ellipsoid
                        if cut_option[ll_index] == Cut::PositiveCut(onshell_prop_index)
                            || cut_option[ll_index] == Cut::NegativeCut(onshell_prop_index)
                        {
                            continue;
                        }

                        let mut surface_shift: LorentzVector<float> = onshell_prop.q.cast();
                        for (&sign, q) in sig_ll_in_cb.iter().zip_eq(cut_shift.iter()) {
                            surface_shift -= q.multiply_sign(sign);
                        }

                        // now see if we have an ellipsoid
                        // 1. surface_shift != 0
                        // 2. surface_shift^2 - (sum_i m_i)^2 >= 0
                        // 3. all signs need to be the same (except 0)
                        // 4. surface_shift.t needs to have the opposite sign as in step 3.
                        if surface_shift.t != float::zero() {
                            // multiply the residue sign
                            let mut surface_signs = sig_ll_in_cb.clone();
                            for (ss, sign) in surface_signs.iter_mut().zip_eq(&cut_residue_sign) {
                                if *sign < 0 {
                                    *ss *= -1;
                                }
                            }

                            let mut cut_mass_sum = float::zero();
                            for (ss, mass) in surface_signs.iter().zip_eq(cut_mass.iter()) {
                                cut_mass_sum += mass.multiply_sign(*ss);
                            }

                            let surface_mass =
                                float::from_f64(onshell_prop.m_squared.sqrt()).unwrap();

                            let group = self.surfaces.len(); // every surface is in a different group at first
                            let surface_sign_sum = surface_signs.iter().sum::<i8>();
                            let surface_signs_abs_sum =
                                surface_signs.iter().map(|x| x.abs()).sum::<i8>();
                            // check the two branches Delta+ and Delta-
                            for &delta_sign in &[-1, 1] {
                                // if all non-zero coefficients are the same, we have an ellipse
                                // otherwise, we have a hyperboloid
                                if (surface_sign_sum + delta_sign).abs()
                                    == surface_signs_abs_sum + delta_sign.abs()
                                {
                                    if surface_shift.square()
                                        - (cut_mass_sum.abs() + surface_mass).powi(2)
                                        >= Into::<float>::into(-1e-13 * self.e_cm_squared)
                                        && surface_shift.t.multiply_sign(delta_sign) < float::zero()
                                    {
                                        let is_pinch = surface_shift.square().abs()
                                            < Into::<float>::into(1e-13 * self.e_cm_squared);

                                        self.surfaces.push(Surface {
                                            group,
                                            surface_type: if is_pinch {
                                                SurfaceType::Pinch
                                            } else {
                                                SurfaceType::Ellipsoid
                                            },
                                            cut_structure_index: cut_index,
                                            cut_option_index,
                                            cut: cut_option.clone(),
                                            onshell_ll_index: ll_index,
                                            onshell_prop_index,
                                            delta_sign,
                                            sig_ll_in_cb: sig_ll_in_cb.iter().cloned().collect(),
                                            signs: surface_signs.iter().cloned().collect(),
                                            shift: surface_shift.clone(),
                                            id: vec![],
                                        });
                                    }
                                } else {
                                    let mut neg_surface_signs_count =
                                        surface_signs.iter().filter(|x| **x == -1).count();
                                    let mut pos_surface_signs_count =
                                        surface_signs.iter().filter(|x| **x == 1).count();

                                    if delta_sign > 0 {
                                        pos_surface_signs_count += 1;
                                    } else {
                                        neg_surface_signs_count += 1;
                                    }

                                    // if we have at least two positive and two negative foci,
                                    // the hyperboloid always exists. If we have one negative focus,
                                    // we require the surface equation to be negative at qi=0.
                                    // For one positive focus, we require it to be positive.
                                    // For the case with one positive and one negative focus, we
                                    // use a known condition.
                                    if match (pos_surface_signs_count, neg_surface_signs_count) {
                                        (1, 1) => {
                                            surface_shift.square()
                                                - (surface_mass - cut_mass_sum.abs()).powi(2)
                                                <= Into::<float>::into(-1e-13 * self.e_cm_squared)
                                        }
                                        (1, _) | (_, 1) => {
                                            let mut eval = float::zero();
                                            for (&ss, mass) in
                                                surface_signs.iter().zip_eq(cut_mass.iter())
                                            {
                                                if pos_surface_signs_count == 1 && ss == 1
                                                    || neg_surface_signs_count == 1 && ss == -1
                                                {
                                                    eval += (surface_shift.spatial_squared()
                                                        + mass * mass)
                                                        .sqrt()
                                                        .multiply_sign(ss);
                                                } else {
                                                    eval += mass.multiply_sign(ss);
                                                }
                                            }

                                            if pos_surface_signs_count == 1 && delta_sign == 1
                                                || neg_surface_signs_count == 1 && delta_sign == -1
                                            {
                                                eval += (surface_shift.spatial_squared()
                                                    + Into::<float>::into(onshell_prop.m_squared))
                                                .sqrt()
                                                .multiply_sign(delta_sign);
                                            } else {
                                                eval += Into::<float>::into(onshell_prop.m_squared)
                                                    .sqrt()
                                                    .multiply_sign(delta_sign);
                                            }

                                            eval += surface_shift.t;

                                            pos_surface_signs_count == 1 && eval > float::zero()
                                                || neg_surface_signs_count == 1
                                                    && eval < float::zero()
                                        }
                                        _ => true,
                                    } {
                                        self.surfaces.push(Surface {
                                            group,
                                            surface_type: SurfaceType::Hyperboloid,
                                            cut_structure_index: cut_index,
                                            cut_option_index,
                                            cut: cut_option.clone(),
                                            onshell_ll_index: ll_index,
                                            onshell_prop_index,
                                            delta_sign,
                                            sig_ll_in_cb: sig_ll_in_cb.iter().cloned().collect(),
                                            signs: surface_signs.iter().cloned().collect(),
                                            shift: surface_shift.clone(),
                                            id: vec![],
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if self.settings.general.debug > 1 {
            println!("Surfaces for {}:", self.name);
        }

        // Identify similar surfaces and put them in the same group
        // If a surface is the first of a new group, the group id will be the index
        // in the surface list
        let mut unique_ellipsoids = 0;
        let mut group_representative: Vec<(Vec<((usize, usize), i8, i8)>, usize)> = vec![];
        for (surf_index, s) in &mut self.surfaces.iter_mut().enumerate() {
            // create a vector of tuples that identifies a surface:
            // [(x_1, a_1, b_1)] where (x,a,b) means a*E_x+b*p_x^0
            // the sorted list is a unique ID when set the a_1=1 by
            // an overall sign multiplication
            let mut s_cut_sorted = s
                .cut
                .iter()
                .enumerate()
                .filter(|(_, x)| **x != Cut::NoCut)
                .zip_eq(s.signs.iter())
                .zip_eq(s.sig_ll_in_cb.iter())
                .filter_map(|(((lli, c), s), ss)| match c {
                    Cut::PositiveCut(i) | Cut::NegativeCut(i) if *s != 0 => {
                        Some(((lli, *i), *s, -*ss))
                    }
                    _ => None,
                })
                .collect::<Vec<_>>();

            // now add the surface focus. it always as +p^0
            s_cut_sorted.push(((s.onshell_ll_index, s.onshell_prop_index), s.delta_sign, 1));
            s_cut_sorted.sort();

            // normalize such that the first focus has positive sign
            if s_cut_sorted[0].1 == -1 {
                for focus in s_cut_sorted.iter_mut() {
                    *focus = (focus.0, -focus.1, -focus.2);
                }
            }

            s.id = s_cut_sorted.clone();

            match group_representative.iter().find(|r| r.0 == s_cut_sorted) {
                Some(i) => s.group = i.1,
                None => {
                    if s.surface_type == SurfaceType::Ellipsoid {
                        unique_ellipsoids += 1;
                    }
                    s.group = surf_index;
                    group_representative.push((s_cut_sorted, surf_index));
                }
            }

            if self.settings.general.debug > 1 {
                println!(
                    "  | id={}, group={}, type={:?}, prop={:?} cut={}, id={:?}, shift={}",
                    surf_index,
                    s.group,
                    s.surface_type,
                    (s.onshell_ll_index, s.onshell_prop_index),
                    CutList(&s.cut),
                    s.id,
                    s.shift
                );
            }
        }

        if self.settings.general.debug > 1 {
            println!("Number of unique ellipsoids: {}", unique_ellipsoids);
            println!("Surfaces not appearing in cut:");
        }

        // make a list of all the ellipsoids that do not appear in a specific cut
        let mut ellipsoids_not_in_cuts = vec![];
        for cut_options in self.ltd_cut_options.iter() {
            let mut accum = vec![];
            for cut in cut_options {
                let mut ellipsoids_not_in_cut = vec![];

                for (surf_index, s) in self.surfaces.iter().enumerate() {
                    // only go through ellipsoid representatives
                    if s.surface_type != SurfaceType::Ellipsoid || s.group != surf_index {
                        continue;
                    }

                    let mut in_cut = false;
                    for s1 in &self.surfaces[surf_index..] {
                        if s1.group == surf_index {
                            if s1.cut == *cut {
                                in_cut = true;
                                break;
                            }
                        }
                    }
                    if !in_cut {
                        ellipsoids_not_in_cut.push(surf_index);
                    }
                }

                if self.settings.general.debug > 1 {
                    println!("  | {}: {:?}", CutList(cut), ellipsoids_not_in_cut);
                }
                accum.push(ellipsoids_not_in_cut);
            }
            ellipsoids_not_in_cuts.push(accum);
        }
        self.ellipsoids_not_in_cuts = ellipsoids_not_in_cuts;

        // now find all dual canceling groups
        let mut dual_groups_rep: Vec<(Vec<usize>, usize)> = vec![];
        let mut dual_groups = UnionFind::new(self.ltd_cut_options.iter().map(|x| x.len()).sum());
        for s in &self.surfaces {
            if s.surface_type != SurfaceType::Hyperboloid {
                continue;
            }

            let mut cs = s
                .cut
                .iter()
                .zip_eq(&self.loop_lines)
                .filter_map(|(c, ll)| {
                    if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = c {
                        Some(ll.propagators[*i].id)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();

            cs.push(self.loop_lines[s.onshell_ll_index].propagators[s.onshell_prop_index].id);
            cs.sort();

            // get the cut option index of the cut
            let cut_index = self.ltd_cut_options[..s.cut_structure_index]
                .iter()
                .map(|x| x.len())
                .sum::<usize>()
                + s.cut_option_index;

            let mut new = true;
            for x in &mut dual_groups_rep {
                if x.0 == cs {
                    new = false;
                    dual_groups.union(x.1, cut_index);
                    break;
                }
            }

            if new {
                dual_groups_rep.push((cs, cut_index));
            }
        }

        if self.settings.general.debug > 1 {
            let mut cut_index = 0;
            let mut cut_opt_index = 0;
            println!("Dual grouping:");
            for x in &dual_groups.to_vec() {
                println!(
                    "  | {}: {}",
                    CutList(&self.ltd_cut_options[cut_index][cut_opt_index]),
                    x
                );

                cut_opt_index += 1;
                if self.ltd_cut_options[cut_index].len() == cut_opt_index {
                    cut_index += 1;
                    cut_opt_index = 0;
                }
            }
        }

        self.all_excluded_surfaces = vec![false; self.surfaces.len()];
        for d in &mut self.fixed_deformation {
            let mut found = false;
            for surf_id in &d.excluded_surface_ids {
                for (i, s) in self.surfaces.iter().enumerate() {
                    if s.id == *surf_id && i == s.group {
                        d.excluded_surface_indices.push(i);
                        self.all_excluded_surfaces[i] = true;
                        found = true;
                    }
                }

                if !found {
                    panic!("Unkown surface id in fixed deformation: {:?}", surf_id);
                }
            }
        }

        if self.settings.general.deformation_strategy == DeformationStrategy::Fixed
            && self.fixed_deformation.is_empty()
        {
            panic!("Fixed deformation strategy selected but none was provided.");
        }
    }

    /// Map a vector in the unit hypercube to the infinite hypercube.
    /// Also compute the Jacobian.
    pub fn parameterize<T: FloatLike>(&self, x: &[f64], loop_index: usize) -> ([T; 3], T) {
        let e_cm = Into::<T>::into(self.e_cm_squared).sqrt()
            * Into::<T>::into(self.settings.parameterization.shifts[loop_index].0);
        let mut l_space = [T::zero(); 3];
        let mut jac = T::one();

        // rescale the input to the desired range
        let mut x_r = [0.; 3];
        for (xd, xi, &(lo, hi)) in izip!(
            &mut x_r,
            x,
            &self.settings.parameterization.input_rescaling[loop_index]
        ) {
            *xd = lo + xi * (hi - lo);
            jac *= Into::<T>::into(hi - lo);
        }

        match self.settings.parameterization.mode {
            ParameterizationMode::Cartesian => match self.settings.parameterization.mapping {
                ParameterizationMapping::Log => {
                    for i in 0..3 {
                        let x = Into::<T>::into(x_r[i]);
                        l_space[i] = e_cm * (x / (T::one() - x)).ln();
                        jac *= e_cm / (x - x * x);
                    }
                }
                ParameterizationMapping::Linear => {
                    for i in 0..3 {
                        let x = Into::<T>::into(x_r[i]);
                        l_space[i] = e_cm * (T::one() / (T::one() - x) - T::one() / x);
                        jac *= e_cm
                            * (T::one() / (x * x) + T::one() / ((T::one() - x) * (T::one() - x)));
                    }
                }
            },
            ParameterizationMode::Spherical => {
                let radius = match self.settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        // r = e_cm * ln(1 + b*x/(1-x))
                        let x = Into::<T>::into(x_r[0]);
                        let b = Into::<T>::into(self.settings.parameterization.b);
                        let radius = e_cm * (T::one() + b * x / (T::one() - x)).ln();
                        jac *= e_cm * b / (T::one() - x) / (T::one() + x * (b - T::one()));

                        radius
                    }
                    ParameterizationMapping::Linear => {
                        // r = e_cm * b * x/(1-x)
                        let b = Into::<T>::into(self.settings.parameterization.b);
                        let radius = e_cm * b * Into::<T>::into(x_r[0])
                            / (T::one() - Into::<T>::into(x_r[0]));
                        jac *= <T as num_traits::Float>::powi(e_cm * b + radius, 2) / e_cm / b;
                        radius
                    }
                };
                let phi = Into::<T>::into(2.) * <T as FloatConst>::PI() * Into::<T>::into(x_r[1]);
                jac *= Into::<T>::into(2.) * <T as FloatConst>::PI();

                let cos_theta = -T::one() + Into::<T>::into(2.) * Into::<T>::into(x_r[2]); // out of range
                jac *= Into::<T>::into(2.);
                let sin_theta = (T::one() - cos_theta * cos_theta).sqrt();

                l_space[0] = radius * sin_theta * phi.cos();
                l_space[1] = radius * sin_theta * phi.sin();
                l_space[2] = radius * cos_theta;

                jac *= radius * radius; // spherical coord
            }
        }

        // add a shift such that k=l is harder to be picked up by integrators such as cuhre
        l_space[0] += e_cm * Into::<T>::into(self.settings.parameterization.shifts[loop_index].1);
        l_space[1] += e_cm * Into::<T>::into(self.settings.parameterization.shifts[loop_index].2);
        l_space[2] += e_cm * Into::<T>::into(self.settings.parameterization.shifts[loop_index].3);

        (l_space, jac)
    }

    pub fn inv_parametrize<T: FloatLike>(
        &self,
        mom: &LorentzVector<f64>,
        loop_index: usize,
    ) -> ([T; 3], T) {
        if self.settings.parameterization.mode != ParameterizationMode::Spherical {
            panic!("Inverse mapping is only implemented for spherical coordinates");
        }

        let mut jac = T::one();
        let e_cm = Into::<T>::into(self.e_cm_squared).sqrt()
            * Into::<T>::into(self.settings.parameterization.shifts[loop_index].0);

        let x: T = Into::<T>::into(mom.x)
            - e_cm * Into::<T>::into(self.settings.parameterization.shifts[loop_index].1);
        let y: T = Into::<T>::into(mom.y)
            - e_cm * Into::<T>::into(self.settings.parameterization.shifts[loop_index].2);
        let z: T = Into::<T>::into(mom.z)
            - e_cm * Into::<T>::into(self.settings.parameterization.shifts[loop_index].3);

        let k_r_sq = x * x + y * y + z * z;
        let k_r = k_r_sq.sqrt();

        let x2 = if y < T::zero() {
            T::one() + Into::<T>::into(0.5) * T::FRAC_1_PI() * T::atan2(y, x)
        } else {
            Into::<T>::into(0.5) * T::FRAC_1_PI() * T::atan2(y, x)
        };

        // cover the degenerate case
        if k_r_sq.is_zero() {
            return ([T::zero(), x2, T::zero()], T::zero());
        }

        let x1 = match self.settings.parameterization.mapping {
            ParameterizationMapping::Log => {
                let b = Into::<T>::into(self.settings.parameterization.b);
                let x1 = T::one() - b / (-T::one() + b + (k_r / e_cm).exp());
                jac /= e_cm * b / (T::one() - x1) / (T::one() + x1 * (b - T::one()));
                x1
            }
            ParameterizationMapping::Linear => {
                let b = Into::<T>::into(self.settings.parameterization.b);
                jac /= <T as num_traits::Float>::powi(e_cm * b + k_r, 2) / e_cm / b;
                k_r / (e_cm * b + k_r)
            }
        };

        let x3 = Into::<T>::into(0.5) * (T::one() + z / k_r);

        jac /= Into::<T>::into(2.) * <T as FloatConst>::PI();
        jac /= Into::<T>::into(2.);
        jac /= k_r * k_r;

        let mut x = [x1, x2, x3];
        for (xi, &(lo, hi)) in x
            .iter_mut()
            .zip_eq(&self.settings.parameterization.input_rescaling[loop_index])
        {
            *xi = (*xi - Into::<T>::into(lo)) / Into::<T>::into(hi - lo);
            jac /= Into::<T>::into(hi - lo);
        }

        (x, jac)
    }

    #[inline]
    fn compute_lambda_factor<U: DimName, T: FloatLike>(
        x: DualN<T, U>,
        y: DualN<T, U>,
    ) -> DualN<T, U>
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
    {
        // FIXME: not smooth
        if x * Into::<T>::into(2.) < y {
            y * Into::<T>::into(0.25)
        } else if y < T::zero() {
            x - y * Into::<T>::into(0.5)
        } else {
            x - y * Into::<T>::into(0.25)
        }
    }

    /// Determine the global lambda scaling by going through each propagator
    /// for each cut and taking the minimum of their respective lambdas
    /// Additionally, check for the expansion condition and make sure that
    /// the real part of the cut propagator is positive
    fn determine_lambda<U: DimName, T: FloatLike>(
        &self,
        kappas: &[LorentzVector<DualN<T, U>>],
        lambda_max: f64,
        cache: &mut LTDCache<T>,
    ) -> [DualN<T, U>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut lambda_sq = DualN::from_real(<T as Float>::powi(Into::<T>::into(lambda_max), 2));

        let sigma = DualN::from_real(Into::<T>::into(
            self.settings.deformation.scaling.softmin_sigma,
        ));
        let mut smooth_min_num = lambda_sq * (-lambda_sq / sigma).exp();
        let mut smooth_min_den = (-lambda_sq / sigma).exp();

        // update the cut information with the kappa-dependent information
        let cache = cache.get_cache_mut();

        let mut lambda_i_sq = [lambda_sq; MAX_LOOP];

        for ll in &self.loop_lines {
            let mut kappa_cut = LorentzVector::default();
            for (kappa, &sign) in kappas.iter().zip_eq(ll.signature.iter()) {
                kappa_cut += kappa.multiply_sign(sign);
            }

            for p in &ll.propagators {
                cache.computed_cut_ll[p.id] = 0;

                let info = &mut cache.cut_info[p.id];
                info.kappa = kappa_cut;
                info.kappa_sq = kappa_cut.spatial_squared_impr();
                info.kappa_dot_mom = kappa_cut.spatial_dot_impr(&info.momentum);

                if self.settings.deformation.scaling.non_cut_propagator_check {
                    info.a = (info.kappa_sq * info.real_energy.powi(2)
                        - info.kappa_dot_mom.powi(2))
                        / info.real_energy.powi(3);
                    info.b = info.kappa_dot_mom / info.real_energy;
                    info.c = info.real_energy;
                }

                if info.kappa_sq.is_zero() {
                    continue;
                }

                // we have to make sure that our linear expansion for the deformation vectors is reasonable
                // for that we need lambda < c * (q_i^2^cut+m_i^2)/|kappa_i^cut * q_i^cut|
                if self.settings.deformation.scaling.expansion_check {
                    let lambda_exp = DualN::from_real(Into::<T>::into(
                        self.settings.deformation.scaling.expansion_threshold,
                    )) * info.spatial_and_mass_sq
                        / info.kappa_dot_mom.abs(); // note: not holomorphic

                    let lambda_exp_sq = lambda_exp * lambda_exp;

                    if sigma.is_zero() {
                        if lambda_exp_sq < lambda_sq {
                            lambda_sq = lambda_exp_sq;
                        }
                    } else {
                        let e = (-lambda_exp_sq / sigma).exp();
                        smooth_min_num += lambda_exp_sq * e;
                        smooth_min_den += e;
                    }
                }

                // prevent a discontinuity in the cut delta by making sure that the real part of the cut propagator is > 0
                // for that we need lambda_sq < (q_i^2^cut+m_i^2)/(kappa_i^cut^2)
                if self.settings.deformation.scaling.positive_cut_check {
                    let lambda_disc_sq = info.spatial_and_mass_sq / info.kappa_sq
                        * DualN::from_real(Into::<T>::into(0.95));

                    if sigma.is_zero() {
                        let loop_var_count = Into::<T>::into(
                            ll.signature.iter().filter(|i| **i != 0).count() as f64,
                        );
                        for (lambda_ii, loop_mom) in lambda_i_sq.iter_mut().zip(&ll.signature) {
                            if *loop_mom != 0 {
                                if lambda_disc_sq / loop_var_count / loop_var_count < *lambda_ii {
                                    *lambda_ii = lambda_disc_sq / loop_var_count / loop_var_count;
                                }
                            }
                        }
                    } else {
                        // FIXME: not supported
                        let e = (-lambda_disc_sq / sigma).exp();
                        smooth_min_num += lambda_disc_sq * e;
                        smooth_min_den += e;
                    }
                }

                // check the on-shell propagator for zeros
                if self.settings.deformation.scaling.cut_propagator_check {
                    let x = (info.kappa_dot_mom / info.kappa_sq).powi(2);
                    let y = info.spatial_and_mass_sq / info.kappa_sq;

                    let prop_lambda_sq = Topology::compute_lambda_factor(x, y);

                    if sigma.is_zero() {
                        if prop_lambda_sq < lambda_sq {
                            lambda_sq = prop_lambda_sq;
                        }
                    } else {
                        let e = (-prop_lambda_sq / sigma).exp();
                        smooth_min_num += prop_lambda_sq * e;
                        smooth_min_den += e;
                    }
                }
            }
        }

        if self.settings.deformation.scaling.non_cut_propagator_check {
            let mut cut_infos = [&cache.cut_info[0]; MAX_LOOP];

            for (cuts, cb_to_lmb_mat) in self
                .ltd_cut_options
                .iter()
                .zip_eq(self.cb_to_lmb_mat.iter())
            {
                for cut in cuts {
                    let mut index = 0;
                    for (c, ll) in cut.iter().zip_eq(self.loop_lines.iter()) {
                        if let Cut::NegativeCut(i) | Cut::PositiveCut(i) = c {
                            cut_infos[index] = &cache.cut_info[ll.propagators[*i].id];
                            index += 1;
                        }
                    }

                    for (ll_cut, onshell_ll) in cut.iter().zip_eq(self.loop_lines.iter()) {
                        // store a bit flag to see if we have checked a loop line with a cut already
                        // 1 => positive cut computed, 2 => negative cut computed
                        match *ll_cut {
                            Cut::PositiveCut(cut_prop) => {
                                if cache.computed_cut_ll[onshell_ll.propagators[cut_prop].id] & 1
                                    == 0
                                {
                                    cache.computed_cut_ll[onshell_ll.propagators[cut_prop].id] |= 1;
                                } else {
                                    continue;
                                }
                            }
                            Cut::NegativeCut(cut_prop) => {
                                if cache.computed_cut_ll[onshell_ll.propagators[cut_prop].id] & 2
                                    == 0
                                {
                                    cache.computed_cut_ll[onshell_ll.propagators[cut_prop].id] |= 2;
                                } else {
                                    continue;
                                }
                            }
                            Cut::NoCut => {}
                        }

                        // TODO: this should be precomputed
                        // determine the map from loop line momentum to cut momenta
                        let mut onshell_signs = [0; MAX_LOOP];
                        for ii in 0..self.n_loops {
                            for jj in 0..self.n_loops {
                                // note we are taking the transpose of cb_to_lmb_mat
                                onshell_signs[ii] += cb_to_lmb_mat[jj * self.n_loops + ii]
                                    * onshell_ll.signature[jj];
                            }
                        }

                        // multiply the signature by the sign of the cut
                        let mut index = 0;
                        for x in cut {
                            match x {
                                Cut::PositiveCut(_) => index += 1,
                                Cut::NegativeCut(_) => {
                                    onshell_signs[index] *= -1;
                                    index += 1;
                                }
                                Cut::NoCut => {}
                            }
                        }

                        let mut ellipse_flag = 0;
                        for x in &onshell_signs {
                            match x {
                                1 => ellipse_flag |= 1,
                                -1 => ellipse_flag |= 2,
                                0 => {}
                                _ => unreachable!(),
                            }
                        }

                        if ellipse_flag == 3 && self.settings.deformation.scaling.skip_hyperboloids
                        {
                            continue;
                        }

                        // treat the cut part of the surface equation
                        let mut a = DualN::from_real(T::zero());
                        let mut b = DualN::from_real(T::zero());
                        let mut c = DualN::from_real(T::zero());
                        for (cut_info, &sign) in cut_infos[..self.n_loops]
                            .iter()
                            .zip_eq(onshell_signs[..self.n_loops].iter())
                        {
                            if sign != 0 {
                                a += cut_info.a.multiply_sign(sign);
                                b += cut_info.b.multiply_sign(sign);
                                c += cut_info.c.multiply_sign(sign);
                            }
                        }

                        // TODO: we are checking surfaces more than once. The ones that depend on only
                        // one loop momentum for example
                        for (prop_index, onshell_prop) in onshell_ll.propagators.iter().enumerate()
                        {
                            let on_shell_info = &cache.cut_info[onshell_prop.id];

                            if *ll_cut == Cut::PositiveCut(prop_index)
                                || *ll_cut == Cut::NegativeCut(prop_index)
                            {
                                continue;
                            }

                            // construct the on-shell part of the propagator
                            for &sign in &[-T::one(), T::one()] {
                                if self.settings.deformation.scaling.skip_hyperboloids
                                    && (ellipse_flag == 1 && !sign.is_one()
                                        || ellipse_flag == 2 && sign.is_one())
                                {
                                    continue;
                                }

                                let a_tot = (a + on_shell_info.a * sign) * Into::<T>::into(-0.5);
                                let a_tot_inv = a_tot.inv();
                                let b_tot = b + on_shell_info.b * sign;
                                let c_tot = c + on_shell_info.c * sign + on_shell_info.shift.t;
                                let x = (b_tot * a_tot_inv).powi(2) * Into::<T>::into(0.25);
                                let y = -c_tot * a_tot_inv;

                                let prop_lambda_sq = Topology::compute_lambda_factor(x, y);

                                if a_tot.real().is_zero() {
                                    continue;
                                }

                                if sigma.is_zero() {
                                    if prop_lambda_sq < lambda_sq {
                                        lambda_sq = prop_lambda_sq;
                                    }
                                } else {
                                    let e = (-prop_lambda_sq / sigma).exp();
                                    smooth_min_num += prop_lambda_sq * e;
                                    smooth_min_den += e;
                                }
                            }
                        }
                    }
                }
            }
        }

        for lambda_ii in &mut lambda_i_sq[..self.n_loops] {
            if lambda_sq < *lambda_ii {
                *lambda_ii = lambda_sq;
            }

            *lambda_ii = lambda_ii.sqrt();
        }

        if sigma.is_zero() {
            lambda_i_sq
        } else {
            unimplemented!()
            //(smooth_min_num / smooth_min_den).sqrt()
        }
    }

    /// Compute the normal vector for each ellipsoid.
    /// Also evaluates the ellipsoid if `full_propagator` is false,
    /// otherwise it evaluates the full propagator.
    fn compute_ellipsoid_deformation_vector<U: DimName, T: FloatLike>(
        &self,
        full_propagator: bool,
        cache: &mut LTDCache<T>,
    ) -> usize
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut cut_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut deform_dirs = [LorentzVector::default(); MAX_LOOP];

        let cache = cache.get_cache_mut();
        let kappa_surf = &mut cache.deform_dirs;
        let inv_surf_prop = &mut cache.ellipsoid_eval;

        // surface equation:
        // m*|sum_i^L a_i q_i^cut + p_vec| + sum_i^L a_i*r_i * |q_i^cut| - p^0 = 0
        // where m=delta_sign a = sig_ll_in_cb, r = residue_sign
        let mut group_counter = 0;
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            // only deform the set of different ellipsoids
            if surf.surface_type != SurfaceType::Ellipsoid || surf.group != surf_index {
                continue;
            }

            // compute v_i = q_i^cut / sqrt(|q_i^cut|^2 + m^2) and the cut energy
            // construct the normalized 3-momenta that flow through the cut propagators
            // the 0th component is to be ignored
            let mut cut_counter = 0;
            let mut cut_energy = DualN::default();
            for (c, ll) in surf.cut.iter().zip_eq(self.loop_lines.iter()) {
                if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = c {
                    cut_dirs[cut_counter] = cache.cut_info[ll.propagators[*i].id].momentum
                        / cache.cut_info[ll.propagators[*i].id].real_energy;

                    if let Cut::PositiveCut(i) = c {
                        cut_energy += cache.cut_info[ll.propagators[*i].id]
                            .real_energy
                            .multiply_sign(surf.sig_ll_in_cb[cut_counter]);
                    } else {
                        cut_energy -= cache.cut_info[ll.propagators[*i].id]
                            .real_energy
                            .multiply_sign(surf.sig_ll_in_cb[cut_counter]);
                    }

                    cut_counter += 1;
                }
            }

            // compute w=(sum_i^L a_i q_i^cut + p_vec) / sqrt(|sum_i^L a_i q_i^cut + p_vec|^2 + m^2)
            // determine sum_i^L a_i q_i^cut using the loop-momentum basis
            let surface_prop =
                &self.loop_lines[surf.onshell_ll_index].propagators[surf.onshell_prop_index];
            let surface_dir_norm = cache.cut_info[surface_prop.id].momentum
                / cache.cut_info[surface_prop.id].real_energy;

            // compute n_i = (v_i + a_i * w) * abs(a_i)
            // n_i is the normal vector of the ith entry of the cut momentum basis
            // this n_i is 0 when the surface does not depend on the ith cut loop line
            // we update cut_dirs from v_i to n_i
            for (cut_dir, &sign) in cut_dirs[..self.n_loops]
                .iter_mut()
                .zip_eq(surf.sig_ll_in_cb.iter())
            {
                if sign != 0 {
                    *cut_dir += surface_dir_norm.multiply_sign(sign);
                }
            }

            // convert from cut momentum basis to loop momentum basis
            for (i, deform_dir) in deform_dirs[..self.n_loops].iter_mut().enumerate() {
                *deform_dir = LorentzVector::default();
                for (&sign, cut_dir) in self.cb_to_lmb_mat[surf.cut_structure_index]
                    [i * self.n_loops..(i + 1) * self.n_loops]
                    .iter()
                    .zip_eq(cut_dirs[..self.n_loops].iter())
                {
                    *deform_dir += cut_dir.multiply_sign(sign);
                }
            }

            // evaluate the inverse propagator of the surface
            if full_propagator {
                inv_surf_prop[group_counter] = (cut_energy + Into::<T>::into(surf.shift.t)).powi(2)
                    - cache.cut_info[surface_prop.id].spatial_and_mass_sq;
            } else {
                inv_surf_prop[group_counter] = cut_energy
                    + Into::<T>::into(surf.shift.t)
                    + cache.cut_info[surface_prop.id]
                        .real_energy
                        .multiply_sign(surf.delta_sign);
            }

            for (loop_index, dir) in deform_dirs[..self.n_loops].iter().enumerate() {
                // note the sign
                kappa_surf[group_counter * self.n_loops + loop_index] =
                    -dir / dir.spatial_squared_impr().sqrt();
            }

            group_counter += 1;
        }

        group_counter
    }

    /// Construct a deformation vector by going through all the ellipsoids
    fn deform_ellipsoids<U: DimName, T: FloatLike>(
        &self,
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let group_counter = self.compute_ellipsoid_deformation_vector(true, cache);

        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let cache = cache.get_cache_mut();
        let kappa_surf = &mut cache.deform_dirs;
        let inv_surf_prop = &mut cache.ellipsoid_eval;

        // now combine the kappas from the surface using the chosen strategy
        let mut aij: DualN<T, U> = NumCast::from(self.settings.deformation.additive.a_ij).unwrap();
        let mut lambda: DualN<T, U> = DualN::one();

        for (i, &inv) in inv_surf_prop[..group_counter].iter().enumerate() {
            if i < self.settings.deformation.additive.a_ijs.len() {
                aij = NumCast::from(self.settings.deformation.additive.a_ijs[i]).unwrap();
            }

            if i < self.settings.deformation.lambdas.len() {
                lambda = NumCast::from(self.settings.deformation.lambdas[i]).unwrap();
            }

            let dampening = match self.settings.deformation.additive.mode {
                AdditiveMode::Exponential => (-inv * inv
                    / (aij * <T as Float>::powi(Into::<T>::into(self.e_cm_squared), 2)))
                .exp(),
                AdditiveMode::Hyperbolic => {
                    let t = inv * inv / <T as Float>::powi(Into::<T>::into(self.e_cm_squared), 2);
                    t / (t + aij)
                }
                AdditiveMode::Unity => DualN::one(),
                AdditiveMode::SoftMin => unimplemented!(),
            };

            for (loop_index, kappa) in kappas[..self.n_loops].iter_mut().enumerate() {
                let dir = kappa_surf[i * self.n_loops + loop_index];
                if self.settings.general.debug > 2 {
                    println!(
                        "Deformation contribution for surface {}:\n  | dir={}\n  | exp={}",
                        group_counter,
                        dir,
                        dampening.real()
                    );
                }

                *kappa += dir * dampening * lambda;
            }
        }

        kappas
    }

    /// Compute the deformation for an ellipsoid and its (potential)
    /// intersections with other ellipsoids.
    fn deform_intersections<U: DimName, T: FloatLike>(
        &self,
        ellipse_id: usize,
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut kappas = [LorentzVector::default(); MAX_LOOP];

        // determine all vectors per ellipsoid once at the start
        if ellipse_id == 0 {
            self.compute_ellipsoid_deformation_vector(false, cache);
        }

        let cache = cache.get_cache_mut();

        // the vector of the surface is always active
        // it can get interpolated with vectors from other surfaces near intersections
        for (loop_index, kappa) in kappas[..self.n_loops].iter_mut().enumerate() {
            *kappa = cache.deform_dirs[ellipse_id * self.n_loops + loop_index];
        }

        // TODO: filter non-intersecting ellipsoids
        // TODO: allow for triple intersections
        let num_ellipsoids = self
            .surfaces
            .iter()
            .enumerate()
            .filter(|(i, x)| x.surface_type == SurfaceType::Ellipsoid && x.group == *i)
            .count();
        for j in 0..num_ellipsoids {
            if ellipse_id == j {
                continue;
            }

            // this factor will determine how far away from an intersection we are
            let j_eval = cache.ellipsoid_eval[j].powi(2);

            for (loop_index, kappa) in kappas[..self.n_loops].iter_mut().enumerate() {
                let dampening = DualN::from_real(T::one())
                    - (j_eval / (j_eval + Into::<T>::into(self.settings.deformation.fixed.m_ij)));
                *kappa += cache.deform_dirs[j * self.n_loops + loop_index] * dampening;
            }
        }

        kappas
    }

    fn get_deformation_for_cut<U: DimName, T: FloatLike>(
        &self,
        cut: &[Cut],
        cut_structure_index: usize,
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut deform_dirs = [LorentzVector::default(); MAX_LOOP];
        let cache = cache.get_cache();

        // TODO: precompute this map by mapping cut indices to cut propagator ids
        let mut index = 0;
        let mut cut_infos = [&cache.cut_info[0]; MAX_LOOP];
        for (ll_cut, ll) in cut.iter().zip_eq(self.loop_lines.iter()) {
            if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = *ll_cut {
                cut_infos[index] = &cache.cut_info[ll.propagators[i].id];
                index += 1;
            }
        }

        // transform the cut momenta with shift to loop momentum basis
        for i in 0..self.n_loops {
            let mut deform_dir = LorentzVector::default();
            for (&sign, &cut_info) in self.cb_to_lmb_mat[cut_structure_index]
                [i * self.n_loops..(i + 1) * self.n_loops]
                .iter()
                .zip_eq(cut_infos.iter())
            {
                let mom_normalized = cut_info.momentum / cut_info.real_energy;
                deform_dir += mom_normalized.multiply_sign(sign);
            }

            deform_dirs[i] = deform_dir;
        }

        deform_dirs
    }

    /// Evaluate a surface with complex momenta.
    fn evaluate_surface_complex<T: FloatLike>(
        &self,
        surf: &Surface,
        loop_momenta: &[LorentzVector<num::Complex<T>>],
    ) -> num::Complex<T> {
        // TODO: use a cache
        let mut res = num::Complex::zero();

        let mut cut_index = 0;
        for (cut, ll) in surf.cut.iter().zip_eq(self.loop_lines.iter()) {
            if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = *cut {
                if surf.signs[cut_index] == 0 {
                    cut_index += 1;
                    continue;
                }

                // construct the cut energy
                let mut mom = LorentzVector::<num::Complex<T>>::default();
                for (&cut_sign, lm) in ll.signature.iter().zip_eq(loop_momenta.iter()) {
                    mom += lm.multiply_sign(cut_sign);
                }
                // compute the postive cut energy
                let q: LorentzVector<T> = ll.propagators[i].q.cast();
                let energy = ((mom + q).spatial_squared()
                    + <T as NumCast>::from(ll.propagators[i].m_squared).unwrap())
                .sqrt();

                res += energy.multiply_sign(surf.signs[cut_index]);

                cut_index += 1;
            }
        }

        // now for the surface term
        let mut mom = LorentzVector::<num::Complex<T>>::default();
        let onshell_ll = &self.loop_lines[surf.onshell_ll_index];
        let onshell_prop = &onshell_ll.propagators[surf.onshell_prop_index];
        for (&surf_sign, lm) in onshell_ll.signature.iter().zip_eq(loop_momenta.iter()) {
            mom += lm.multiply_sign(surf_sign);
        }

        let q: LorentzVector<T> = onshell_prop.q.cast();
        let energy = ((mom + q).spatial_squared()
            + <T as NumCast>::from(onshell_prop.m_squared).unwrap())
        .sqrt();

        res += energy.multiply_sign(surf.delta_sign);
        res += <T as NumCast>::from(surf.shift.t).unwrap();
        res
    }

    /// Construct a constant deformation vector.
    /// If an ellipse id is specified, it will be a different constant per ellipse.
    fn deform_constant<U: DimName, T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        ellipsoid_id: Option<usize>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut kappas = [LorentzVector::default(); MAX_LOOP];

        for i in 0..self.n_loops {
            kappas[i] = loop_momenta[i] / loop_momenta[i].spatial_squared_impr().sqrt();

            if let Some(e) = ellipsoid_id {
                kappas[i] *= DualN::from_real(Into::<T>::into(e as f64 + 1.));
            }
        }

        kappas
    }

    fn deform_fixed<U: DimName, T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let cache = cache.get_cache_mut();

        // evaluate all excluded ellipsoids
        for (i, surf) in self.surfaces.iter().enumerate() {
            if surf.surface_type != SurfaceType::Ellipsoid
                || surf.group != i
                || !self.all_excluded_surfaces[i]
            {
                cache.ellipsoid_eval[i] = DualN::zero();
                continue;
            }

            let mut cut_counter = 0;
            let mut cut_energy = DualN::default();
            for (c, ll) in surf.cut.iter().zip_eq(self.loop_lines.iter()) {
                if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = c {
                    if let Cut::PositiveCut(i) = c {
                        cut_energy += cache.cut_info[ll.propagators[*i].id]
                            .real_energy
                            .multiply_sign(surf.sig_ll_in_cb[cut_counter]);
                    } else {
                        cut_energy -= cache.cut_info[ll.propagators[*i].id]
                            .real_energy
                            .multiply_sign(surf.sig_ll_in_cb[cut_counter]);
                    }

                    cut_counter += 1;
                }
            }

            let surface_prop =
                &self.loop_lines[surf.onshell_ll_index].propagators[surf.onshell_prop_index];
            cache.ellipsoid_eval[i] = cut_energy
                + Into::<T>::into(surf.shift.t)
                + cache.cut_info[surface_prop.id]
                    .real_energy
                    .multiply_sign(surf.delta_sign);

            // TODO: properly normalize
            // take the square so that the number if always positive
            cache.ellipsoid_eval[i] =
                cache.ellipsoid_eval[i].powi(2) / Into::<T>::into(self.e_cm_squared);
        }

        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let mut lambda = DualN::one();
        let mut mij = Into::<T>::into(self.settings.deformation.fixed.m_ij);

        for (i, d) in self.fixed_deformation.iter().enumerate() {
            if i < self.settings.deformation.lambdas.len() {
                lambda = NumCast::from(self.settings.deformation.lambdas[i]).unwrap();
            }

            let mut s = DualN::one();
            let mut softmin_num = DualN::zero();
            let mut softmin_den = DualN::zero();

            for &surf_index in &d.excluded_surface_indices {
                if surf_index < self.settings.deformation.fixed.m_ijs.len() {
                    mij = Into::<T>::into(self.settings.deformation.fixed.m_ijs[surf_index]);
                }

                // t is the weighing factor that is 0 if we are on both cut_i and cut_j
                // at the same time and goes to 1 otherwise
                debug_assert!(self.surfaces[surf_index].surface_type == SurfaceType::Ellipsoid);
                let t = cache.ellipsoid_eval[surf_index];

                let sup = if self.settings.deformation.fixed.mode == AdditiveMode::SoftMin {
                    if self.settings.deformation.fixed.sigma.is_zero() {
                        s = s.min(t);
                        s
                    } else {
                        let e = (-t / Into::<T>::into(self.settings.deformation.fixed.sigma)).exp();
                        softmin_num += t * e;
                        softmin_den += e;
                        e
                    }
                } else {
                    let sup = t / (t + mij);
                    s *= sup;
                    sup
                };

                if self.settings.general.debug > 2 {
                    println!(
                        "  | surf {}: t={:e}, suppression={:e}",
                        surf_index,
                        t.real(),
                        sup.real()
                    );
                }

                if self.settings.deformation.fixed.mode == AdditiveMode::SoftMin {
                    if !self.settings.deformation.fixed.sigma.is_zero() {
                        s = softmin_num / softmin_den;
                    }

                    if !self.settings.deformation.fixed.m_ij.is_zero() {
                        s = s / (s + mij);
                    }
                }
            }

            if self.settings.general.debug > 2 {
                println!(
                    "  | k={:e}\n  | dirs={:e}\n  | suppression={:e}\n  | contribution={:e}",
                    loop_momenta[0].real(),
                    &d.deformation_sources[0],
                    s.real(),
                    &d.deformation_sources[0].cast() * s.real(),
                );
            }

            for ii in 0..self.n_loops {
                let dir = loop_momenta[ii] - d.deformation_sources[ii].cast();
                kappas[ii] -= dir / dir.spatial_distance() * s * lambda;
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{}={:e}\n", ii + 1, kappas[ii].real());
            }
        }

        kappas
    }

    fn deform_generic<U: DimName, T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        cut: Option<(usize, usize)>,
        ellipsoid_id: Option<usize>,
        cache: &mut LTDCache<T>,
    ) -> ([LorentzVector<T>; MAX_LOOP], Complex<T>)
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        // compute all cut energies
        if self.settings.general.deformation_strategy != DeformationStrategy::Constant
            || self.settings.deformation.scaling.lambda > 0.
        {
            let cache = cache.get_cache_mut();
            for ll in &self.loop_lines {
                let mut mom = LorentzVector::default();
                for (l, &c) in loop_momenta.iter().zip_eq(ll.signature.iter()) {
                    mom += l.multiply_sign(c);
                }

                // update the cut info
                for p in &ll.propagators {
                    let q: LorentzVector<DualN<T, U>> = p.q.cast();
                    let mass = DualN::from_real(Into::<T>::into(p.m_squared));
                    let cm = (mom + q).spatial_squared() + mass;
                    let energy = cm.sqrt();
                    cache.cut_energies[p.id] = energy;

                    let info = &mut cache.cut_info[p.id];
                    info.id = p.id;
                    info.momentum = mom + q;
                    info.real_energy = energy;
                    info.spatial_and_mass_sq = cm;
                    // TODO: these two never change and can be set at the start!
                    info.shift = q;
                    info.mass = mass;
                }
            }
        }

        let mut kappas = match self.settings.general.deformation_strategy {
            DeformationStrategy::Duals => {
                let co = cut.unwrap();
                let cut = &self.ltd_cut_options[co.0][co.1];
                self.get_deformation_for_cut(cut, co.0, cache)
            }
            DeformationStrategy::Intersections => {
                self.deform_intersections(ellipsoid_id.unwrap(), cache)
            }
            DeformationStrategy::Additive => self.deform_ellipsoids(cache),
            DeformationStrategy::Constant => self.deform_constant(loop_momenta, ellipsoid_id),
            DeformationStrategy::Fixed => self.deform_fixed(loop_momenta, cache),
            DeformationStrategy::None => unreachable!(),
        };

        // make sure the kappa has the right dimension by multiplying in the scale
        let scale = DualN::from_real(
            T::from_f64(
                self.e_cm_squared.sqrt() * self.settings.deformation.overall_scaling_constant,
            )
            .unwrap(),
        );
        for (kappa, k) in &mut kappas[..self.n_loops]
            .iter_mut()
            .zip_eq(loop_momenta.iter())
        {
            match self.settings.deformation.overall_scaling {
                OverallDeformationScaling::Constant => *kappa *= scale,
                OverallDeformationScaling::Linear => {
                    *kappa *= k.spatial_squared_impr().sqrt()
                        * Into::<T>::into(self.settings.deformation.overall_scaling_constant)
                }
                OverallDeformationScaling::Sigmoid => {
                    let k_scale = k.spatial_squared_impr().sqrt();
                    *kappa *=
                        k_scale * Into::<T>::into(2.) / (DualN::one() + (k_scale / scale).exp());
                }
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{} scaled={:e}\n", ii + 1, kappas[ii].real());
            }
        }

        let lambda = if self.settings.deformation.scaling.lambda > 0. {
            self.determine_lambda(
                &kappas[..self.n_loops],
                self.settings.deformation.scaling.lambda,
                cache,
            )
        } else {
            [NumCast::from(self.settings.deformation.scaling.lambda.abs()).unwrap(); MAX_LOOP]
        };

        for (k, lambda_i) in kappas[..self.n_loops]
            .iter_mut()
            .zip_eq(&lambda[..self.n_loops])
        {
            if self.settings.general.debug > 2 {
                println!("  | lambda={:e}\n", lambda_i.real());
            }

            *k *= *lambda_i;
            k.t = DualN::zero(); // make sure we do not have a left-over deformation
        }

        let jac_mat = &mut cache.get_cache_mut().deformation_jacobian;
        for i in 0..3 * self.n_loops {
            for j in 0..3 * self.n_loops {
                // first index: loop momentum, second: xyz, third: dual
                jac_mat[i * 3 * self.n_loops + j] =
                    Complex::new(T::zero(), kappas[i / 3][i % 3 + 1][j + 1]);
            }
            jac_mat[i * 3 * self.n_loops + i] += Complex::one();
        }

        let jac = utils::determinant(jac_mat, 3 * self.n_loops);
        // take the real part
        let mut r = [LorentzVector::default(); MAX_LOOP];
        for (rr, k) in r[..self.n_loops].iter_mut().zip_eq(&kappas[..self.n_loops]) {
            *rr = k.real();
        }
        (r, jac)
    }

    pub fn deform<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<T>],
        cut: Option<(usize, usize)>,
        ellipsoid_id: Option<usize>,
        cache: &mut LTDCache<T>,
    ) -> ([LorentzVector<T>; MAX_LOOP], Complex<T>) {
        if DeformationStrategy::None == self.settings.general.deformation_strategy {
            let r = [LorentzVector::default(); MAX_LOOP];
            return (r, Complex::new(T::one(), T::zero()));
        }

        match self.n_loops {
            1 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];

                r[0] = loop_momenta[0].map(|x| Dual4::from_real(x));
                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                }

                return self.deform_generic(&r[..self.n_loops], cut, ellipsoid_id, cache);
            }
            2 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual7::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual7::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                }
                self.deform_generic(&r[..self.n_loops], cut, ellipsoid_id, cache)
            }
            3 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual10::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual10::from_real(x));
                r[2] = loop_momenta[2].map(|x| Dual10::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                    r[2][i + 1][i + 7] = T::one();
                }
                self.deform_generic(&r[..self.n_loops], cut, ellipsoid_id, cache)
            }
            4 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual13::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual13::from_real(x));
                r[2] = loop_momenta[2].map(|x| Dual13::from_real(x));
                r[3] = loop_momenta[3].map(|x| Dual13::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                    r[2][i + 1][i + 7] = T::one();
                    r[3][i + 1][i + 10] = T::one();
                }
                self.deform_generic(&r[..self.n_loops], cut, ellipsoid_id, cache)
            }
            5 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual16::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual16::from_real(x));
                r[2] = loop_momenta[2].map(|x| Dual16::from_real(x));
                r[3] = loop_momenta[3].map(|x| Dual16::from_real(x));
                r[4] = loop_momenta[4].map(|x| Dual16::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                    r[2][i + 1][i + 7] = T::one();
                    r[3][i + 1][i + 10] = T::one();
                    r[4][i + 1][i + 13] = T::one();
                }
                self.deform_generic(&r[..self.n_loops], cut, ellipsoid_id, cache)
            }
            6 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual19::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual19::from_real(x));
                r[2] = loop_momenta[2].map(|x| Dual19::from_real(x));
                r[3] = loop_momenta[3].map(|x| Dual19::from_real(x));
                r[4] = loop_momenta[4].map(|x| Dual19::from_real(x));
                r[5] = loop_momenta[5].map(|x| Dual19::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                    r[2][i + 1][i + 7] = T::one();
                    r[3][i + 1][i + 10] = T::one();
                    r[4][i + 1][i + 13] = T::one();
                    r[5][i + 1][i + 16] = T::one();
                }
                self.deform_generic(&r[..self.n_loops], cut, ellipsoid_id, cache)
            }
            n => panic!("Binding for deformation at {} loops is not implemented", n),
        }
    }

    /// Set the energy component of the loop momenta according to
    /// `cut`. It takes the cut energies from the cache.
    #[inline]
    pub fn set_loop_momentum_energies<T: FloatLike>(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &LTDCache<T>,
    ) {
        // compute the cut energy for each loop line
        let mut cut_energy = [Complex::default(); MAX_LOOP];
        let mut index = 0;
        for (&ll_cut, ll) in cut.iter().zip_eq(self.loop_lines.iter()) {
            if let Cut::PositiveCut(j) | Cut::NegativeCut(j) = ll_cut {
                let e = cache.complex_cut_energies[ll.propagators[j].id];
                if let Cut::PositiveCut(_) = ll_cut {
                    cut_energy[index] = e - Into::<T>::into(ll.propagators[j].q.t);
                } else {
                    cut_energy[index] = -e - Into::<T>::into(ll.propagators[j].q.t);
                }
                index += 1;
            }
        }

        // compute the energies for the loop momenta
        for (i, l) in k_def.iter_mut().enumerate() {
            l.t = Complex::default();
            for (c, e) in mat[i * self.n_loops..(i + 1) * self.n_loops]
                .iter()
                .zip_eq(&cut_energy[..self.n_loops])
            {
                l.t += e.multiply_sign(*c);
            }
        }
    }

    #[inline]
    pub fn evaluate_cut<T: FloatLike>(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &'static str> {
        self.set_loop_momentum_energies(k_def, cut, mat, cache);

        let mut r = Complex::new(T::one(), T::zero());
        for (i, (ll_cut, ll)) in cut.iter().zip_eq(self.loop_lines.iter()).enumerate() {
            // get the loop line result from the cache if possible
            r *= match ll_cut {
                Cut::PositiveCut(j) => cache.complex_loop_line_eval[i][*j][0],
                Cut::NegativeCut(j) => cache.complex_loop_line_eval[i][*j][1],
                _ => ll.evaluate(&k_def, ll_cut, &self, cache)?,
            };
        }
        r = r.inv(); // normal inverse may overflow but has better precision than finv, which overflows later

        if self.settings.general.debug > 1 {
            match self.n_loops {
                1 => {
                    println!(
                        "Cut {}:\n  | result={:e}\n  | k0={:e}",
                        CutList(cut),
                        r,
                        k_def[0].t,
                    );
                }
                2 => {
                    println!(
                        "Cut {}:\n  | result={:e}\n  | k0={:e}\n  | l0={:e}",
                        CutList(cut),
                        r,
                        k_def[0].t,
                        k_def[1].t,
                    );
                }
                3 => {
                    println!(
                        "Cut {}:\n  | result={:e}\n  | k0={:e}\n  | l0={:e}\n  | m0={:e}",
                        CutList(cut),
                        r,
                        k_def[0].t,
                        k_def[1].t,
                        k_def[2].t,
                    );
                }
                4 => {
                    println!(
                        "Cut {}:\n  | result={:e}\n  | k0={:e}\n  | l0={:e}\n  | m0={:e}\n  | n0={:e}",
                        CutList(cut),
                        r,
                        k_def[0].t,
                        k_def[1].t,
                        k_def[2].t,
                        k_def[3].t,
                    );
                }
                _ => {}
            }
        }

        Ok(r)
    }

    /// Compute the complex cut energies and evaluate the cut loop lines
    pub fn compute_complex_cut_energies<T: FloatLike>(
        &self,
        k_def: &ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cache: &mut LTDCache<T>,
    ) -> Result<(), &'static str> {
        // compute all complex cut energies
        for (ll_index, ll) in self.loop_lines.iter().enumerate() {
            let mut mom = LorentzVector::default();
            for (l, &c) in k_def.iter().zip_eq(ll.signature.iter()) {
                mom += l.multiply_sign(c);
            }

            for p in &ll.propagators {
                let q: LorentzVector<Complex<T>> = p.q.cast();
                let cm = (mom + q).spatial_squared() + Into::<T>::into(p.m_squared);

                if self.settings.deformation.scaling.positive_cut_check
                    && cm.re < T::zero()
                    && cm.im < T::zero()
                {
                    eprintln!(
                        "{} for prop {}, ll sig={:?}, ks={:?}: {}",
                        "Branch cut detected".red(),
                        p.id,
                        ll.signature,
                        k_def,
                        cm
                    );
                }

                cache.complex_prop_spatial[p.id] = cm;
                cache.complex_cut_energies[p.id] = cm.sqrt();
            }

            let has_positive_cut = self.ltd_cut_structure.iter().any(|x| x[ll_index] == 1);
            let has_negative_cut = self.ltd_cut_structure.iter().any(|x| x[ll_index] == -1);

            for i in 0..ll.propagators.len() {
                // compute the entire dual loop line
                // note that the deformed momenta are not used here
                if has_positive_cut {
                    cache.complex_loop_line_eval[ll_index][i][0] =
                        ll.evaluate(k_def, &Cut::PositiveCut(i), &self, cache)?;
                }
                if has_negative_cut {
                    cache.complex_loop_line_eval[ll_index][i][1] =
                        ll.evaluate(k_def, &Cut::NegativeCut(i), &self, cache)?;
                }
            }
        }
        Ok(())
    }

    pub fn evaluate_multi_channel<T: FloatLike>(
        &self,
        k: &[LorentzVector<T>],
        cache: &mut LTDCache<T>,
        python_numerator: &Option<PythonNumerator>,
    ) -> Result<(Complex<T>, ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>), &'static str> {
        let num_ellipsoids = self
            .surfaces
            .iter()
            .enumerate()
            .filter(|(i, x)| x.surface_type == SurfaceType::Ellipsoid && x.group == *i)
            .count();
        let mut complex_evals_sq = vec![Complex::default(); num_ellipsoids];
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]> = ArrayVec::default();

        let mut res = Complex::default();

        let mut sum_fac: Complex<T> = Complex::default();
        for i in 0..num_ellipsoids {
            // compute the deformation vector for this ellipsoid
            let (kappas, jac) = self.deform(k, None, Some(i), cache);
            k_def = (0..self.n_loops)
                .map(|i| {
                    k[i].map(|x| Complex::new(x, T::zero()))
                        + kappas[i].map(|x| Complex::new(T::zero(), x))
                })
                .collect();

            self.compute_complex_cut_energies(&k_def, cache)?;

            // compute multi-channeling factors by evaluating every surface
            let mut index = 0;
            let mut suppresion_factor = Complex::one();
            for (ei, e) in self.surfaces.iter().enumerate() {
                if e.surface_type == SurfaceType::Ellipsoid && e.group == ei {
                    complex_evals_sq[index] =
                        utils::powi(self.evaluate_surface_complex(e, &k_def), 2);
                    suppresion_factor *= complex_evals_sq[index];
                    index += 1;
                }
            }

            let mut normalization = Complex::zero();
            for j in 0..num_ellipsoids {
                normalization += suppresion_factor / complex_evals_sq[j];
            }

            let (r, kd) = self.evaluate_all_dual_integrands(k, k_def, cache, python_numerator)?;
            k_def = kd;

            res += r * jac * suppresion_factor / complex_evals_sq[i] / normalization;
            sum_fac += suppresion_factor / complex_evals_sq[i] / normalization;
        }

        Ok((res, k_def))
    }

    fn evaluate_all_dual_integrands<T: FloatLike>(
        &self,
        k: &[LorentzVector<T>],
        mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cache: &mut LTDCache<T>,
        python_numerator: &Option<PythonNumerator>,
    ) -> Result<(Complex<T>, ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>), &'static str> {
        // evaluate all dual integrands
        let mut result = Complex::default();
        let mut cut_counter = 0;
        for (cut_structure_index, (cuts, mat)) in self
            .ltd_cut_options
            .iter()
            .zip_eq(self.cb_to_lmb_mat.iter())
            .enumerate()
        {
            // for each cut coming from the same cut structure
            for (cut_option_index, cut) in cuts.iter().enumerate() {
                if self.settings.general.cut_filter.is_empty()
                    || self.settings.general.cut_filter.contains(&cut_counter)
                {
                    let mut dual_jac_def = Complex::one();
                    if self.settings.general.deformation_strategy == DeformationStrategy::Duals {
                        let (kappas, jac) = self.deform(
                            &k,
                            Some((cut_structure_index, cut_option_index)),
                            None,
                            cache,
                        );
                        k_def = (0..self.n_loops)
                            .map(|i| {
                                k[i].map(|x| Complex::new(x, T::zero()))
                                    + kappas[i].map(|x| Complex::new(T::zero(), x))
                            })
                            .collect();
                        dual_jac_def = jac;
                        self.compute_complex_cut_energies(&k_def, cache)?;
                    }

                    let v = self.evaluate_amplitude_cut(&mut k_def, cut, mat, cache)?;

                    // k_def has the correct energy component at this stage
                    if let Some(pn) = python_numerator {
                        result += v * pn.evaluate_numerator(&k_def[..self.n_loops]) * dual_jac_def
                    } else {
                        result += v * dual_jac_def
                    }
                }
                cut_counter += 1;
            }
        }
        Ok((result, k_def))
    }

    #[inline]
    pub fn evaluate<'a, T: FloatLike>(
        &self,
        x: &'a [f64],
        cache: &mut LTDCache<T>,
        python_numerator: &Option<PythonNumerator>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        T,
        Complex<T>,
        Complex<T>,
    ) {
        // parameterize
        let mut k = [LorentzVector::default(); MAX_LOOP];
        let mut jac_para = T::one();
        for i in 0..self.n_loops {
            // set the loop index to i + 1 so that we can also shift k
            let (l_space, jac) = self.parameterize(&x[i * 3..(i + 1) * 3], i);

            // there could be some rounding here
            let rot = self.rotation_matrix;
            k[i] = LorentzVector::from_args(
                T::zero(),
                <T as NumCast>::from(rot[0][0]).unwrap() * l_space[0]
                    + <T as NumCast>::from(rot[0][1]).unwrap() * l_space[1]
                    + <T as NumCast>::from(rot[0][2]).unwrap() * l_space[2],
                <T as NumCast>::from(rot[1][0]).unwrap() * l_space[0]
                    + <T as NumCast>::from(rot[1][1]).unwrap() * l_space[1]
                    + <T as NumCast>::from(rot[1][2]).unwrap() * l_space[2],
                <T as NumCast>::from(rot[2][0]).unwrap() * l_space[0]
                    + <T as NumCast>::from(rot[2][1]).unwrap() * l_space[1]
                    + <T as NumCast>::from(rot[2][2]).unwrap() * l_space[2],
            );
            jac_para *= jac;
        }

        // deform
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]> = ArrayVec::default();
        let mut jac_def = Complex::one();

        let r = if self.settings.general.partial_fractioning {
            self.evaluate_multi_channel(&k, cache, python_numerator)
        } else {
            if self.settings.general.deformation_strategy != DeformationStrategy::Duals {
                let (kappas, jac) = self.deform(&k, None, None, cache);
                k_def = (0..self.n_loops)
                    .map(|i| {
                        k[i].map(|x| Complex::new(x, T::zero()))
                            + kappas[i].map(|x| Complex::new(T::zero(), x))
                    })
                    .collect();
                jac_def = jac;

                if self.compute_complex_cut_energies(&k_def, cache).is_err() {
                    return (x, k_def, jac_para, jac_def, Complex::default());
                }
            }

            self.evaluate_all_dual_integrands(&k, k_def, cache, python_numerator)
        };

        let (mut result, k_def) = match r {
            Ok(v) => v,
            Err(_) => {
                return (
                    x,
                    ArrayVec::default(),
                    jac_para,
                    jac_def,
                    Complex::default(),
                )
            }
        };

        result *= utils::powi(
            num::Complex::new(T::zero(), Into::<T>::into(-2.) * <T as FloatConst>::PI()),
            self.n_loops,
        ); // factor of delta cut

        result *= utils::powi(
            num::Complex::new(
                Into::<T>::into(1.)
                    / <T as Float>::powi(Into::<T>::into(2.) * <T as FloatConst>::PI(), 4),
                T::zero(),
            ),
            self.n_loops,
        ); // loop momentum factor
        result *= jac_def * jac_para;

        if self.settings.general.use_amplitude {
            (
                x,
                k_def,
                jac_para,
                jac_def,
                result
                    + Complex::new(
                        Into::<T>::into(self.amplitude.int_ct[0].re),
                        Into::<T>::into(self.amplitude.int_ct[0].im),
                    ),
            )
        } else {
            (x, k_def, jac_para, jac_def, result)
        }
    }
}
