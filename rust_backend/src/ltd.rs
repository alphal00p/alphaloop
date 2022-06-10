use crate::topologies::{
    CacheSelector, Cut, CutList, LTDCache, SOCPProblem, Surface, SurfaceType, Topology,
};
use crate::utils::Signum;
use crate::{
    AdditiveMode, DeformationStrategy, ExpansionCheckStrategy, FloatLike, IRHandling,
    OverallDeformationScaling, ParameterizationMapping, ParameterizationMode, PoleCheckStrategy,
    Settings, MAX_LOOP,
};
use hyperdual::Hyperdual;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use num::Complex;
use num_traits::ops::inv::Inv;
use num_traits::{Float, FloatConst, NumCast, One, Pow, Signed, Zero};

use crate::utils;

impl Topology {
    /// Compute LTD related quantities for this topology.
    /// If `external_momenta_set` is true, the existence conditions
    /// for ellipsoids and e_cm are determined as well.
    pub fn process(&mut self, external_momenta_set: bool) {
        assert!(
            self.n_loops <= MAX_LOOP,
            "MAX_LOOP is too small: it should be at least {}",
            self.n_loops
        );

        // set the identity rotation matrix
        self.rotation_matrix = [
            [f64::one(), f64::zero(), f64::zero()],
            [f64::zero(), f64::one(), f64::zero()],
            [f64::zero(), f64::zero(), f64::one()],
        ];

        // copy the signature to the propagators
        let mut prop_id = 0;
        self.propagator_id_to_ll_id = vec![];
        for (ll_index, l) in self.loop_lines.iter_mut().enumerate() {
            for (p_index, p) in l.propagators.iter_mut().enumerate() {
                p.signature = l.signature.clone();
                p.id = prop_id;
                self.propagator_id_to_ll_id.push((ll_index, p_index));
                prop_id += 1;
            }
        }

        if external_momenta_set {
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
                self.e_cm_squared = self.external_kinematics[0].spatial_squared_impr().abs();
            }

            // determine the on-shell flag
            for (i, e) in self.external_kinematics.iter().enumerate() {
                if e.square_impr().abs() < 1e-10 * self.e_cm_squared {
                    self.on_shell_flag |= 2_usize.pow(i as u32);
                }
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

        // Construct all subspaces required for the deformation.
        // We store the cut structure index, so that we can use the basis transformation to go into
        // the subspace.
        self.subspaces.clear();
        for n in 0..self.n_loops {
            for subspace in (0..self.loop_lines.len()).combinations(n) {
                let cut_structure_index = match self
                    .ltd_cut_structure
                    .iter()
                    .enumerate()
                    .filter(|(_i, c)| subspace.iter().all(|s| c[*s] != 0))
                    .next()
                {
                    Some(x) => x.0,
                    None => {
                        continue;
                    }
                };

                if subspace.len() == 0 {
                    self.subspaces.push((cut_structure_index, 0, vec![]));
                }

                for subspace_opts in subspace
                    .iter()
                    .map(|&li| (0..self.loop_lines[li].propagators.len()).map(move |pi| (li, pi)))
                    .multi_cartesian_product()
                {
                    let cut_option_index = self.ltd_cut_options[cut_structure_index]
                        .iter()
                        .enumerate()
                        .filter(|(_i, c)| {
                            subspace_opts.iter().all(|&(sli, spi)| {
                                c[sli] == Cut::PositiveCut(spi) || c[sli] == Cut::NegativeCut(spi)
                            })
                        })
                        .map(|(i, _)| i)
                        .next()
                        .unwrap();

                    self.subspaces
                        .push((cut_structure_index, cut_option_index, subspace_opts));
                }
            }
        }

        self.surfaces.clear();
        for (cut_index, (residue_sign, cut_options)) in self
            .ltd_cut_structure
            .iter()
            .zip_eq(self.ltd_cut_options.iter())
            .enumerate()
        {
            let mut cut_signatures_matrix = vec![];
            let mut cut_residue_sign = vec![];

            // solve the linear system for fixing the energy component of each loop line
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
                let mut cut_shift: Vec<LorentzVector<f64>> = vec![]; // qs of cut
                let mut cut_mass = vec![];
                for (cut, ll) in cut_option.iter().zip_eq(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut
                    {
                        cut_shift.push(ll.propagators[*cut_prop_index].q.cast());
                        cut_mass.push(ll.propagators[*cut_prop_index].m_squared.sqrt());
                    }
                }

                for (ll_index, ll) in self.loop_lines.iter().enumerate() {
                    // if this loop line has no loop momentum, skip it
                    if ll.signature.iter().all(|s| *s == 0) {
                        continue;
                    }

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

                        let mut surface_shift: LorentzVector<f64> = onshell_prop.q.cast();
                        for (&sign, q) in sig_ll_in_cb.iter().zip_eq(cut_shift.iter()) {
                            surface_shift -= q.multiply_sign(sign);
                        }

                        // multiply the residue sign
                        let mut surface_signs = sig_ll_in_cb.clone();
                        for (ss, sign) in surface_signs.iter_mut().zip_eq(&cut_residue_sign) {
                            if *sign < 0 {
                                *ss *= -1;
                            }
                        }

                        let mut cut_mass_sum = 0.;
                        for (ss, mass) in surface_signs.iter().zip_eq(cut_mass.iter()) {
                            cut_mass_sum += mass.multiply_sign(*ss);
                        }

                        let surface_mass = onshell_prop.m_squared.sqrt();

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
                                // now see if the ellipsoid exists if the external momenta are set
                                // 1. surface_shift != 0
                                // 2. surface_shift^2 - (sum_i m_i)^2 >= 0
                                // 3. all signs need to be the same (except 0)
                                // 4. surface_shift.t needs to have the opposite sign as in step 3.
                                let mut exists = external_momenta_set;
                                let mut is_pinch = false;

                                if external_momenta_set {
                                    if surface_shift.square()
                                        - (cut_mass_sum.abs() + surface_mass).powi(2)
                                        >= Into::<f64>::into(-1e-13 * self.e_cm_squared)
                                        && surface_shift.t.multiply_sign(delta_sign) < 0.
                                    {
                                        if surface_shift.square()
                                            - (cut_mass_sum.abs() + surface_mass).powi(2)
                                            < Into::<f64>::into(1e-10 * self.e_cm_squared)
                                        {
                                            is_pinch = true;
                                        }
                                    } else {
                                        exists = false;
                                    }
                                }

                                self.surfaces.push(Surface {
                                    unique_ellipsoid_id: None,
                                    exists,
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
                                    massless: (cut_mass_sum + surface_mass).abs()
                                        < 2. * f64::EPSILON,
                                });
                            } else {
                                if !external_momenta_set {
                                    // if the external momenta are unknown, we do not
                                    // store information about hyperboloids
                                    continue;
                                }

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
                                            <= Into::<f64>::into(-1e-13 * self.e_cm_squared)
                                    }
                                    (1, _) | (_, 1) => {
                                        let mut eval = 0.;
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
                                                + Into::<f64>::into(onshell_prop.m_squared))
                                            .sqrt()
                                            .multiply_sign(delta_sign);
                                        } else {
                                            eval += Into::<f64>::into(onshell_prop.m_squared)
                                                .sqrt()
                                                .multiply_sign(delta_sign);
                                        }

                                        eval += surface_shift.t;

                                        pos_surface_signs_count == 1 && eval > 0.
                                            || neg_surface_signs_count == 1 && eval < 0.
                                    }
                                    _ => true,
                                } {
                                    self.surfaces.push(Surface {
                                        unique_ellipsoid_id: None,
                                        exists: true, // we only store existing hyperboloids
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
                                        massless: false, // unused
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }

        if external_momenta_set && self.settings.general.debug > 1 {
            println!("Surfaces for {}:", self.name);
        }

        // Identify similar surfaces and put them in the same group
        // If a surface is the first of a new group, the group id will be the index
        // in the surface list
        let mut unique_ellipsoids = 0;
        let mut total_unique_ellipsoids = 0;
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
                        total_unique_ellipsoids += 1;
                        if s.exists {
                            s.unique_ellipsoid_id = Some(unique_ellipsoids);
                            unique_ellipsoids += 1;
                        }
                    }
                    s.group = surf_index;
                    group_representative.push((s_cut_sorted, surf_index));
                }
            }

            if self.settings.general.debug > 1 && s.exists {
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

        if external_momenta_set && self.settings.general.debug > 1 {
            println!("Number of unique ellipsoids: {}", unique_ellipsoids);
        }

        self.all_excluded_surfaces = vec![false; self.surfaces.len()];
        for d_lim in &mut self.fixed_deformation {
            for d in &mut d_lim.deformation_per_overlap {
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
        }

        if self.settings.general.deformation_strategy == DeformationStrategy::Fixed {
            if total_unique_ellipsoids > 0 {
                let num_propagators = self.loop_lines.iter().map(|x| x.propagators.len()).sum();
                if external_momenta_set {
                    self.socp_problem =
                        SOCPProblem::new(unique_ellipsoids, num_propagators, self.n_loops);
                } else {
                    self.socp_problem = SOCPProblem::new(
                        total_unique_ellipsoids / 2,
                        num_propagators,
                        self.n_loops,
                    );
                }
            }
        }

        if self.settings.general.deformation_strategy == DeformationStrategy::Constant
            && self.constant_deformation.is_none()
        {
            panic!("Constant deformation strategy selected but none was provided.");
        }

        if external_momenta_set {
            if self.settings.general.derive_overlap_structure {
                self.fixed_deformation = self.determine_ellipsoid_overlap_structure(true);
            }
            self.check_fixed_deformation();
        }
    }

    /// Update the ellipsoids when the external momenta have changed.
    pub fn update_ellipsoids(&mut self) {
        assert!(self.e_cm_squared > 0.);

        // now check which ellipsoids exist and update the shifts
        for s in self.surfaces.iter_mut() {
            if s.surface_type == SurfaceType::Hyperboloid {
                continue;
            }

            // determine the shift
            let mut surface_shift: LorentzVector<f64> = self.loop_lines[s.onshell_ll_index]
                .propagators[s.onshell_prop_index]
                .q
                .cast();

            let mut mass_sum: f64 = self.loop_lines[s.onshell_ll_index].propagators
                [s.onshell_prop_index]
                .m_squared
                .sqrt()
                .into();

            let mut cut_index = 0;
            for (cut, ll) in s.cut.iter().zip_eq(self.loop_lines.iter()) {
                if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut {
                    if s.sig_ll_in_cb[cut_index] != 0 {
                        mass_sum +=
                            Into::<f64>::into(ll.propagators[*cut_prop_index].m_squared.sqrt());
                        surface_shift -= ll.propagators[*cut_prop_index]
                            .q
                            .cast()
                            .multiply_sign(s.sig_ll_in_cb[cut_index]);
                    }
                    cut_index += 1;
                }
            }

            // now see if the ellipsoid exists
            // 1. surface_shift != 0
            // 2. surface_shift^2 - (sum_i m_i)^2 >= 0
            // 3. all signs need to be the same (except 0)
            // 4. surface_shift.t needs to have the opposite sign as in step 3.
            let mut exists = true;
            let mut is_pinch = false;
            let size = surface_shift.square() - mass_sum.powi(2);

            if size >= Into::<f64>::into(-1e-13 * self.e_cm_squared)
                && surface_shift.t.multiply_sign(s.delta_sign) < 0.
            {
                if self.settings.general.debug > 0 {
                    if size > 1e-8 * self.e_cm_squared
                        && size < 1e-6 * self.e_cm_squared
                        && !mass_sum.is_zero()
                    {
                        println!("Small ellipsoid {} detected: {}", s.group, size);
                    }
                }

                if size < Into::<f64>::into(1e-8 * self.e_cm_squared) {
                    is_pinch = true;
                }
            } else {
                exists = false;
            }

            if is_pinch {
                s.surface_type = SurfaceType::Pinch;
            } else {
                s.surface_type = SurfaceType::Ellipsoid;
            }
            s.exists = exists;
            s.shift = surface_shift;
        }
    }

    /// Check if the deformation sources satisfy their constraints.
    pub fn check_fixed_deformation(&self) {
        let unique_ellipsoids = self
            .surfaces
            .iter()
            .enumerate()
            .filter(|(i, s)| *i == s.group && s.surface_type == SurfaceType::Ellipsoid && s.exists)
            .count();

        if self.settings.general.deformation_strategy == DeformationStrategy::Fixed
            && unique_ellipsoids > 0
            && self.fixed_deformation.is_empty()
        {
            panic!("Fixed deformation strategy selected but none was provided.");
        }

        for d_lim in &self.fixed_deformation {
            for d in &d_lim.deformation_per_overlap {
                if let Some(overlap) = &d.overlap {
                    if overlap.len() + d.excluded_surface_indices.len() != unique_ellipsoids {
                        println!("Number of ellipsoids between fixed deformation and Rust is different: {} vs {}",
                    overlap.len() + d.excluded_surface_indices.len(), unique_ellipsoids);
                    }
                }

                let loop_momenta = (0..self.n_loops)
                    .map(|i| d.deformation_sources[i].map(|x| Complex::new(x, 0.)))
                    .collect::<Vec<_>>();
                for (surf_index, surf) in self.surfaces.iter().enumerate() {
                    if surf.group == surf_index
                        && surf.exists
                        && surf.surface_type == SurfaceType::Ellipsoid
                        && !d.excluded_surface_indices.contains(&surf_index)
                    {
                        let r = self.evaluate_surface_complex(surf, &loop_momenta);
                        if surf.delta_sign > 0 && r.re >= 0. || surf.delta_sign < 0 && r.re <= 0. {
                            println!(
                                "Deformation source {:?} is not on the inside of surface {}: {}",
                                d.deformation_sources, surf_index, r.re
                            );
                        }
                    }
                }

                // now test if the source is equal to the shift of the propagators
                for (ll_index, prop_index) in &d_lim.excluded_propagators {
                    // determine the sum of sources using the signature
                    let mut source_sum: LorentzVector<f64> = LorentzVector::default();
                    for (sign, source) in self.loop_lines[*ll_index]
                        .signature
                        .iter()
                        .zip_eq(&d.deformation_sources)
                    {
                        source_sum += source.multiply_sign(*sign);
                    }

                    let diff: LorentzVector<f64> =
                        source_sum + self.loop_lines[*ll_index].propagators[*prop_index].q;

                    if diff.spatial_squared() > 1e-10 * self.e_cm_squared {
                        panic!(
                            "Deformation source {:?} is not on a focus of propagator {:?}",
                            d.deformation_sources,
                            (ll_index, prop_index)
                        );
                    }
                }
            }
        }
    }

    #[inline]
    pub fn get_expansion_threshold(&self) -> f64 {
        match self.settings.deformation.scaling.expansion_check_strategy {
            ExpansionCheckStrategy::Ratio => {
                if self.settings.deformation.scaling.expansion_threshold < 0.
                    || self.maximum_ratio_expansion_threshold < 0.
                {
                    self.settings.deformation.scaling.expansion_threshold.abs()
                } else {
                    if self.settings.deformation.scaling.expansion_threshold
                        > 0.5 * self.maximum_ratio_expansion_threshold
                    {
                        0.5 * self.maximum_ratio_expansion_threshold
                    } else {
                        self.settings.deformation.scaling.expansion_threshold
                    }
                }
            }
            _ => self.settings.deformation.scaling.expansion_threshold.abs(),
        }
    }

    #[inline]
    pub fn compute_min_mij(&self) -> f64 {
        // TODO make this quantity static as it does not need to be recomputed statically every
        // time.
        if self.settings.deformation.fixed.m_ij < 0. {
            1.0
        } else {
            let d = self.settings.deformation.fixed.delta;
            let e = self.get_expansion_threshold();
            e * e / ((2.0 - e * e) * (d / (1. - d)).sqrt())
        }
    }

    /// Map a vector in the unit hypercube to the infinite hypercube.
    /// Also compute the Jacobian.
    pub fn parameterize<T: FloatLike>(
        x: &[T],
        e_cm_squared: T,
        loop_index: usize,
        settings: &Settings,
    ) -> ([T; 3], T) {
        let e_cm =
            e_cm_squared.sqrt() * Into::<T>::into(settings.parameterization.shifts[loop_index].0);
        let mut l_space = [T::zero(); 3];
        let mut jac = T::one();

        // rescale the input to the desired range
        let mut x_r = [T::zero(); 3];
        for (xd, xi, &(lo, hi)) in izip!(
            &mut x_r,
            x,
            &settings.parameterization.input_rescaling[loop_index]
        ) {
            let lo = Into::<T>::into(lo);
            let hi = Into::<T>::into(hi);
            *xd = lo + *xi * (hi - lo);
            jac *= Into::<T>::into(hi - lo);
        }

        match settings.parameterization.mode {
            ParameterizationMode::Cartesian => match settings.parameterization.mapping {
                ParameterizationMapping::Log => {
                    for i in 0..3 {
                        let x = x_r[i];
                        l_space[i] = e_cm * (x / (T::one() - x)).ln();
                        jac *= e_cm / (x - x * x);
                    }
                }
                ParameterizationMapping::Linear => {
                    for i in 0..3 {
                        let x = x_r[i];
                        l_space[i] = e_cm * (T::one() / (T::one() - x) - T::one() / x);
                        jac *= e_cm
                            * (T::one() / (x * x) + T::one() / ((T::one() - x) * (T::one() - x)));
                    }
                }
            },
            ParameterizationMode::Spherical => {
                let radius = match settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        // r = e_cm * ln(1 + b*x/(1-x))
                        let x = x_r[0];
                        let b = Into::<T>::into(settings.parameterization.b);
                        let radius = e_cm * (T::one() + b * x / (T::one() - x)).ln();
                        jac *= e_cm * b / (T::one() - x) / (T::one() + x * (b - T::one()));

                        radius
                    }
                    ParameterizationMapping::Linear => {
                        // r = e_cm * b * x/(1-x)
                        let b = Into::<T>::into(settings.parameterization.b);
                        let radius = e_cm * b * x_r[0] / (T::one() - x_r[0]);
                        jac *= <T as num_traits::Float>::powi(e_cm * b + radius, 2) / e_cm / b;
                        radius
                    }
                };
                let phi = Into::<T>::into(2.) * <T as FloatConst>::PI() * x_r[1];
                jac *= Into::<T>::into(2.) * <T as FloatConst>::PI();

                let cos_theta = -T::one() + Into::<T>::into(2.) * x_r[2]; // out of range
                jac *= Into::<T>::into(2.);
                let sin_theta = (T::one() - cos_theta * cos_theta).sqrt();

                l_space[0] = radius * sin_theta * phi.cos();
                l_space[1] = radius * sin_theta * phi.sin();
                l_space[2] = radius * cos_theta;

                jac *= radius * radius; // spherical coord
            }
        }

        // add a shift such that k=l is harder to be picked up by integrators such as cuhre
        l_space[0] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].1);
        l_space[1] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].2);
        l_space[2] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].3);

        (l_space, jac)
    }

    pub fn inv_parametrize<T: FloatLike>(
        mom: &LorentzVector<T>,
        e_cm_squared: T,
        loop_index: usize,
        settings: &Settings,
    ) -> ([T; 3], T) {
        if settings.parameterization.mode != ParameterizationMode::Spherical {
            panic!("Inverse mapping is only implemented for spherical coordinates");
        }

        let mut jac = T::one();
        let e_cm =
            e_cm_squared.sqrt() * Into::<T>::into(settings.parameterization.shifts[loop_index].0);

        let x: T = mom.x - e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].1);
        let y: T = mom.y - e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].2);
        let z: T = mom.z - e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].3);

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

        let x1 = match settings.parameterization.mapping {
            ParameterizationMapping::Log => {
                let b = Into::<T>::into(settings.parameterization.b);
                let x1 = T::one() - b / (-T::one() + b + (k_r / e_cm).exp());
                jac /= e_cm * b / (T::one() - x1) / (T::one() + x1 * (b - T::one()));
                x1
            }
            ParameterizationMapping::Linear => {
                let b = Into::<T>::into(settings.parameterization.b);
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
            .zip_eq(&settings.parameterization.input_rescaling[loop_index])
        {
            *xi = (*xi - Into::<T>::into(lo)) / Into::<T>::into(hi - lo);
            jac /= Into::<T>::into(hi - lo);
        }

        (x, jac)
    }

    #[inline]
    fn compute_lambda_factor<T: FloatLike, const N: usize>(
        x: Hyperdual<T, N>,
        y: Hyperdual<T, N>,
    ) -> Hyperdual<T, N> {
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
    fn determine_lambda<T: FloatLike, const N: usize>(
        &self,
        kappas: &[LorentzVector<Hyperdual<T, N>>],
        lambda_max: f64,
        cache: &mut LTDCache<T>,
    ) -> Hyperdual<T, N>
    where
        LTDCache<T>: CacheSelector<T, N>,
    {
        let mut lambda_sq =
            Hyperdual::<T, N>::from_real(<T as Float>::powi(Into::<T>::into(lambda_max), 2));

        let sigma = Hyperdual::<T, N>::from_real(Into::<T>::into(
            self.settings.deformation.scaling.softmin_sigma,
        ));
        let mut smooth_min_num = lambda_sq * (-lambda_sq / sigma).exp();
        let mut smooth_min_den = (-lambda_sq / sigma).exp();

        // update the cut information with the kappa-dependent information
        let cache = cache.get_cache_mut();

        for ll in &self.loop_lines {
            // skip loop lines without loop dependence
            if ll.signature.iter().all(|&x| x == 0) {
                continue;
            }

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

                if self.settings.deformation.scaling.pole_check_strategy != PoleCheckStrategy::None
                {
                    info.a = (info.kappa_sq * info.real_energy.powi(2)
                        - info.kappa_dot_mom.powi(2))
                        / info.real_energy.powi(3)
                        * Into::<T>::into(-0.5);
                    info.b = info.kappa_dot_mom / info.real_energy;
                    info.c = info.real_energy;
                }

                if info.kappa_sq.is_zero() {
                    continue;
                }

                // we have to make sure that our linear expansion for the deformation vectors is reasonable
                // for that we need:
                // old check: lambda < c * (q_i^2^cut+m_i^2)/|kappa_i^cut * q_i^cut|
                // new check: lambda^2 < (-2*kappa_i.q_i^2 + sqrt(4*kappa_i.q_i^4 + kappa^4 c^2 (q_i^2+m_i^2)^2))/kappa^4
                if self.settings.deformation.scaling.expansion_check_strategy
                    != ExpansionCheckStrategy::None
                {
                    let c = Hyperdual::<T, N>::from_real(Into::<T>::into(
                        self.get_expansion_threshold(),
                    ));

                    let lambda_exp_sq =
                        match self.settings.deformation.scaling.expansion_check_strategy {
                            ExpansionCheckStrategy::FirstLambdaOrder => {
                                let lambda_exp =
                                    c * info.spatial_and_mass_sq / info.kappa_dot_mom.abs(); // note: not holomorphic
                                lambda_exp * lambda_exp
                            }
                            ExpansionCheckStrategy::FullLambdaDependence => {
                                let a = info.kappa_sq * info.kappa_sq;
                                let b = info.kappa_dot_mom * info.kappa_dot_mom;
                                let d = info.spatial_and_mass_sq * info.spatial_and_mass_sq;
                                (-b * Into::<T>::into(2.)
                                    + (b * b * Into::<T>::into(4.) + a * c * c * d).sqrt())
                                    / a
                            }
                            ExpansionCheckStrategy::MagicFudge => {
                                let _a = info.kappa_sq * info.kappa_sq;
                                let b = info.kappa_dot_mom * info.kappa_dot_mom;
                                let d = info.spatial_and_mass_sq * info.spatial_and_mass_sq;
                                c * c * (d / (b + d))
                            }
                            ExpansionCheckStrategy::MagicFudgeWithMin => {
                                let a = info.kappa_sq * info.kappa_sq;
                                let b = info.kappa_dot_mom * info.kappa_dot_mom;
                                let d = info.spatial_and_mass_sq * info.spatial_and_mass_sq;
                                if a < b {
                                    c * c * (d / (b + d))
                                } else {
                                    c * c * (d / (a + d))
                                }
                            }
                            ExpansionCheckStrategy::Ratio => {
                                let a = info.kappa_sq;
                                let b = info.kappa_dot_mom;
                                let d = info.spatial_and_mass_sq;
                                // note that if c < 1, the branch cut check is always satisfied
                                c * c * Into::<T>::into(2.) * d * d / (a * d - b * b)
                            }
                            ExpansionCheckStrategy::None => unreachable!(),
                        };

                    if sigma.is_zero() {
                        // add an abs, since due to numerical instability we may overshoot when doing x-y
                        // in the numerator or denominator.
                        if lambda_exp_sq.abs() < lambda_sq {
                            lambda_sq = lambda_exp_sq.abs();
                        }
                    } else {
                        let e = (-lambda_exp_sq.abs() / sigma).exp();
                        smooth_min_num += lambda_exp_sq.abs() * e;
                        smooth_min_den += e;
                    }
                }

                // prevent a discontinuity in the cut delta by making sure that the real part of the cut propagator is > 0
                // for that we need lambda_sq < (q_i^2^cut+m_i^2)/(kappa_i^cut^2)
                if self.settings.deformation.scaling.branch_cut_check {
                    // use min(mass, UV mass) instead of the actual propagator mass,
                    // such that UV propagators are also protected from crossing the branch cut
                    let mut lambda_disc_sq = info.spatial_and_uv_mass_sq / info.kappa_sq
                        * Hyperdual::<T, N>::from_real(Into::<T>::into(0.95));

                    if self.settings.deformation.scaling.branch_cut_m > 0. {
                        let branch_cut_check_m = Hyperdual::<T, N>::from_real(Into::<T>::into(
                            self.settings.deformation.scaling.branch_cut_m,
                        ));
                        lambda_disc_sq *= Hyperdual::<T, N>::one()
                            + (info.kappa_dot_mom * info.kappa_dot_mom)
                                / (branch_cut_check_m
                                    * Into::<T>::into(self.e_cm_squared)
                                    * Into::<T>::into(self.e_cm_squared));
                    }

                    lambda_disc_sq = lambda_disc_sq.pow(Into::<T>::into(
                        self.settings.deformation.scaling.branch_cut_alpha,
                    ));

                    if sigma.is_zero() {
                        if lambda_disc_sq.abs() < lambda_sq {
                            lambda_sq = lambda_disc_sq.abs();
                        }
                    } else {
                        let e = (-lambda_disc_sq.abs() / sigma).exp();
                        smooth_min_num += lambda_disc_sq.abs() * e;
                        smooth_min_den += e;
                    }
                }
            }
        }

        if self.settings.deformation.scaling.pole_check_strategy != PoleCheckStrategy::None {
            // TODO: rename
            // process all existing and non-existing ellipsoids and pinches
            for (surf_index, s) in self.surfaces.iter().enumerate() {
                if surf_index == s.group && s.surface_type != SurfaceType::Hyperboloid {
                    let mut a = Hyperdual::<T, N>::from_real(T::zero());
                    let mut b = Hyperdual::<T, N>::from_real(T::zero());
                    let mut c = Hyperdual::<T, N>::from_real(T::zero());
                    for ((ll_index, prop_index), delta_sign, os_sign) in &s.id {
                        let info =
                            &cache.cut_info[self.loop_lines[*ll_index].propagators[*prop_index].id];
                        a += info.a;
                        b += info.b;
                        c += info.c + info.shift.t.multiply_sign(os_sign * delta_sign);
                    }

                    let prop_lambda_sq = match self.settings.deformation.scaling.pole_check_strategy
                    {
                        PoleCheckStrategy::Exact => {
                            assert_eq!(self.n_loops, 1);

                            let cut_info = &cache.cut_info
                                [self.loop_lines[(s.id[0].0).0].propagators[(s.id[0].0).1].id];
                            let on_shell_info = &cache.cut_info
                                [self.loop_lines[(s.id[1].0).0].propagators[(s.id[1].0).1].id];

                            // now determine the exact solution for one-loop surfaces
                            let p0 = on_shell_info.shift.t - cut_info.shift.t;
                            let a = cut_info.spatial_and_mass_sq
                                - on_shell_info.spatial_and_mass_sq
                                - p0.powi(2);

                            let b = kappas[0]
                                .spatial_dot(&(cut_info.momentum - on_shell_info.momentum))
                                * Into::<T>::into(2.);

                            let aa = kappas[0].spatial_squared() * p0.powi(2) * Into::<T>::into(4.)
                                - b * b;
                            let bb = (a * b
                                - kappas[0].spatial_dot(&on_shell_info.momentum)
                                    * p0.powi(2)
                                    * Into::<T>::into(4.))
                                * Into::<T>::into(2.);
                            let cc = a * a
                                - on_shell_info.spatial_and_mass_sq
                                    * p0.powi(2)
                                    * Into::<T>::into(4.);

                            let mut lambda_plus = (-Complex::new(Hyperdual::<T, N>::zero(), bb)
                                + (Complex::new(
                                    -bb * bb - aa * cc * Into::<T>::into(4.),
                                    Hyperdual::<T, N>::zero(),
                                ))
                                .sqrt())
                                / (aa * Into::<T>::into(2.));
                            let mut lambda_min = (-Complex::new(Hyperdual::<T, N>::zero(), bb)
                                - (Complex::new(
                                    -bb * bb - aa * cc * Into::<T>::into(4.),
                                    Hyperdual::<T, N>::zero(),
                                ))
                                .sqrt())
                                / (aa * Into::<T>::into(2.));

                            if aa.real().is_zero() {
                                lambda_plus = Complex::new(Hyperdual::<T, N>::zero(), cc / bb);
                                lambda_min = Complex::new(Hyperdual::<T, N>::zero(), cc / bb);
                            }

                            // evaluate the surface with lambda
                            let mut prop_lambda_sq = lambda_sq;
                            for lambda in &[lambda_plus, lambda_min] {
                                let def = kappas[0].real().to_complex(false)
                                    * Complex::new(lambda.re.real(), lambda.im.real());
                                let surf = ((cut_info.momentum.real().to_complex(true) + def)
                                    .spatial_squared()
                                    + cut_info.m_sq.real())
                                .sqrt()
                                    + ((on_shell_info.momentum.real().to_complex(true) + def)
                                        .spatial_squared()
                                        + on_shell_info.m_sq.real())
                                    .sqrt()
                                    + p0.real();

                                if surf.norm() < Into::<T>::into(1e-5) {
                                    let new_lambda = lambda.norm_sqr();
                                    if new_lambda < prop_lambda_sq {
                                        prop_lambda_sq = new_lambda
                                    }
                                }
                            }
                            prop_lambda_sq
                        }
                        PoleCheckStrategy::RealSolution => {
                            let x = (b / a).powi(2) * Into::<T>::into(0.25);
                            let y = -c / a;
                            let prop_lambda_sq = Topology::compute_lambda_factor(x, y);
                            if a.real().is_zero() {
                                continue;
                            }
                            prop_lambda_sq
                        }
                        PoleCheckStrategy::TangentCheck => {
                            let t_out =
                                c / Into::<T>::into(self.settings.deformation.scaling.theta_r_out);
                            let t_in =
                                c * Into::<T>::into(self.settings.deformation.scaling.theta_r_in);

                            let t_c = t_out
                                + b / Into::<T>::into(self.settings.deformation.scaling.theta_c);

                            let rhs = if s.exists {
                                t_out.min(t_in).min(t_c)
                            } else {
                                t_out
                            };

                            (rhs - c) / a
                        }
                        PoleCheckStrategy::None => unreachable!(),
                    };

                    if sigma.is_zero() {
                        if prop_lambda_sq.abs() < lambda_sq {
                            lambda_sq = prop_lambda_sq.abs();
                        }
                    } else {
                        let e = (-prop_lambda_sq.abs() / sigma).exp();
                        smooth_min_num += prop_lambda_sq.abs() * e;
                        smooth_min_den += e;
                    }
                }
            }
        }

        // dampen on soft points
        if self.settings.deformation.scaling.soft_dampening_power > 0. {
            for ll in &self.loop_lines {
                // skip loop lines without loop dependence
                if ll.signature.iter().all(|&x| x == 0) {
                    continue;
                }

                for p in &ll.propagators {
                    let info = &mut cache.cut_info[p.id];

                    if info.m_sq.real() == T::zero() {
                        let t = (info.real_energy.real() * info.real_energy.real()
                            / Into::<T>::into(self.e_cm_squared))
                        .powf(Into::<T>::into(
                            self.settings.deformation.scaling.soft_dampening_power,
                        ));

                        let dampening = Into::<T>::into(lambda_max) * t
                            / (t + Into::<T>::into(self.settings.deformation.fixed.m_ij.powi(2)));

                        if dampening * dampening < lambda_sq.real() {
                            lambda_sq = Hyperdual::from_real(dampening * dampening);
                        }
                    }
                }
            }
        }

        // dampen every pinch and massless E-surface on the line segment between the focal points
        if self.settings.deformation.fixed.dampen_on_pinch
            && self.settings.deformation.fixed.dampen_on_pinch_after_lambda
            && !self.fixed_deformation.is_empty()
        {
            let mij_min_sq = Into::<T>::into(self.compute_min_mij().powi(2));
            let mij_sq =
                Into::<T>::into(self.settings.deformation.fixed.m_ij.abs().powi(2)) * mij_min_sq;
            for (surf_index, s) in self.surfaces.iter().enumerate() {
                if (s.surface_type != SurfaceType::Pinch
                    && s.surface_type != SurfaceType::Ellipsoid)
                    || !s.exists
                    || !s.massless
                {
                    continue;
                }

                let e = cache.ellipsoid_eval[surf_index].unwrap().powi(2)
                    + cache.ellipsoid_normal_norm_eval[surf_index]
                        .unwrap()
                        .powi(2);
                let n = Into::<T>::into(self.surfaces[surf_index].shift.t.powi(2));

                let t = (e
                    / (Into::<T>::into(
                        self.settings.deformation.fixed.pinch_dampening_k_com * self.e_cm_squared,
                    ) + Into::<T>::into(
                        self.settings.deformation.fixed.pinch_dampening_k_shift,
                    ) * n))
                    .pow(Into::<T>::into(
                        self.settings.deformation.fixed.pinch_dampening_alpha,
                    ));

                let sup = t / (t + mij_sq) * Into::<T>::into(lambda_max);

                if sup * sup < lambda_sq {
                    lambda_sq = sup * sup;
                }
            }
        }

        // Setup a veto region at the intersection of pinched and non-pinched E-surfaces.
        if self.settings.deformation.fixed.ir_handling_strategy != IRHandling::None {
            let mut min_ellipse = Hyperdual::<T, N>::from_real(Into::<T>::into(1e99));
            let mut min_pinch = Hyperdual::<T, N>::from_real(Into::<T>::into(1e99));
            for (surf_index, s) in self.surfaces.iter().enumerate() {
                if s.surface_type == SurfaceType::Hyperboloid || !s.exists || surf_index != s.group
                {
                    continue;
                }
                if s.surface_type == SurfaceType::Pinch
                    && self.settings.deformation.fixed.ir_beta_pinch < 0.
                {
                    continue;
                }
                if s.surface_type == SurfaceType::Ellipsoid
                    && self.settings.deformation.fixed.ir_beta_ellipse < 0.
                {
                    continue;
                }

                // TODO: the unwrap_or below is necessary for when a deformation is active but
                // exactly zero because no E-surfaces exist! This case must be better handled!
                let e = cache.ellipsoid_eval[surf_index]
                    .unwrap_or(Hyperdual::<T, N>::from_real(Into::<T>::into(1e99)))
                    .powi(2);
                let n = Into::<T>::into(self.surfaces[surf_index].shift.t.powi(2));

                let t = (e
                    / (Into::<T>::into(
                        self.settings.deformation.fixed.ir_k_com * self.e_cm_squared,
                    ) + Into::<T>::into(self.settings.deformation.fixed.ir_k_shift) * n))
                    .pow(Into::<T>::into(self.settings.deformation.fixed.ir_alpha));

                if s.surface_type == SurfaceType::Ellipsoid && min_ellipse > t {
                    min_ellipse = t;
                }

                if s.surface_type == SurfaceType::Pinch && min_pinch > t {
                    min_pinch = t;
                }
            }

            let mut min_e = Hyperdual::<T, N>::from_real(Into::<T>::into(1e99));
            if self.settings.deformation.fixed.ir_beta_energy >= 0. {
                for ll in self.loop_lines.iter() {
                    if ll.signature.iter().all(|&x| x == 0) {
                        continue;
                    }

                    for p in ll.propagators.iter() {
                        let normalised_e = cache.cut_info[p.id].spatial_and_mass_sq
                            / Into::<T>::into(self.e_cm_squared);
                        if min_e > normalised_e {
                            min_e = normalised_e;
                        }
                    }
                }
            }

            let ir_proximity = if min_pinch.powi(2)
                * Into::<T>::into(self.settings.deformation.fixed.ir_beta_pinch.powi(2))
                > min_e.powi(2)
                    * Into::<T>::into(self.settings.deformation.fixed.ir_beta_energy.powi(2))
            {
                min_e.powi(2)
                    * Into::<T>::into(self.settings.deformation.fixed.ir_beta_energy.powi(2))
            } else {
                min_pinch.powi(2)
                    * Into::<T>::into(self.settings.deformation.fixed.ir_beta_pinch.powi(2))
            } + min_ellipse.powi(2)
                * Into::<T>::into(self.settings.deformation.fixed.ir_beta_ellipse.powi(2));

            if self.settings.deformation.fixed.ir_handling_strategy
                == IRHandling::DismissDeformation
                && self.settings.deformation.fixed.ir_interpolation_length > 0.0
            {
                let sup = (((ir_proximity
                    / Into::<T>::into(self.settings.deformation.fixed.ir_threshold.powi(2)))
                    - Into::<T>::into(1.0))
                    / Into::<T>::into(self.settings.deformation.fixed.ir_interpolation_length))
                .min(Hyperdual::<T, N>::from_real(Into::<T>::into(0.0)))
                .max(Hyperdual::<T, N>::from_real(Into::<T>::into(1.0)));
                if sup * sup < lambda_sq {
                    lambda_sq = sup * sup;
                }
            } else {
                if ir_proximity
                    < Into::<T>::into(self.settings.deformation.fixed.ir_threshold.powi(2))
                {
                    match self.settings.deformation.fixed.ir_handling_strategy {
                        IRHandling::DismissPoint => {
                            return Hyperdual::<T, N>::from_real(T::nan()); // TODO: improve into a proper escalation
                        }
                        IRHandling::DismissDeformation => {
                            return Hyperdual::<T, N>::from_real(Into::<T>::into(0.));
                        }
                        IRHandling::None => {
                            unreachable!();
                        }
                    }
                }
            }
        }

        if sigma.is_zero() {
            lambda_sq.sqrt()
        } else {
            (smooth_min_num / smooth_min_den).sqrt()
        }
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
    fn deform_constant<T: FloatLike, const N: usize>(
        &self,
        loop_momenta: &[LorentzVector<Hyperdual<T, N>>],
    ) -> [LorentzVector<Hyperdual<T, N>>; MAX_LOOP]
    where
        LTDCache<T>: CacheSelector<T, N>,
    {
        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let cd = self.constant_deformation.as_ref().unwrap();

        for i in 0..self.n_loops {
            kappas[i] = loop_momenta[i];

            let rot_matrix = &self.rotation_matrix;

            // rotate the k back to the original frame
            let old_x = kappas[i].x;
            let old_y = kappas[i].y;
            let old_z = kappas[i].z;
            kappas[i].x = old_x * Into::<T>::into(rot_matrix[0][0])
                + old_y * Into::<T>::into(rot_matrix[1][0])
                + old_z * Into::<T>::into(rot_matrix[2][0]);
            kappas[i].y = old_x * Into::<T>::into(rot_matrix[0][1])
                + old_y * Into::<T>::into(rot_matrix[1][1])
                + old_z * Into::<T>::into(rot_matrix[2][1]);
            kappas[i].z = old_x * Into::<T>::into(rot_matrix[0][2])
                + old_y * Into::<T>::into(rot_matrix[1][2])
                + old_z * Into::<T>::into(rot_matrix[2][2]);

            kappas[i] = kappas[i].comp_mul(&cd.alpha[i].cast()) + cd.beta[i].cast();

            // now rotate the constant deformation for rotated topologies
            let old_x = kappas[i].x;
            let old_y = kappas[i].y;
            let old_z = kappas[i].z;
            kappas[i].x = old_x * Into::<T>::into(rot_matrix[0][0])
                + old_y * Into::<T>::into(rot_matrix[0][1])
                + old_z * Into::<T>::into(rot_matrix[0][2]);
            kappas[i].y = old_x * Into::<T>::into(rot_matrix[1][0])
                + old_y * Into::<T>::into(rot_matrix[1][1])
                + old_z * Into::<T>::into(rot_matrix[1][2]);
            kappas[i].z = old_x * Into::<T>::into(rot_matrix[2][0])
                + old_y * Into::<T>::into(rot_matrix[2][1])
                + old_z * Into::<T>::into(rot_matrix[2][2]);
        }

        let mut normalization = Hyperdual::<T, N>::zero();
        for i in 0..self.n_loops {
            normalization += kappas[i].spatial_squared_impr();
        }
        normalization = normalization.sqrt().inv();

        for i in 0..self.n_loops {
            kappas[i] *= normalization;
        }

        kappas
    }

    fn deform_fixed<T: FloatLike, const N: usize>(
        &self,
        loop_momenta: &[LorentzVector<Hyperdual<T, N>>],
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<Hyperdual<T, N>>; MAX_LOOP]
    where
        LTDCache<T>: CacheSelector<T, N>,
    {
        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let mut kappa_source = [LorentzVector::default(); MAX_LOOP];
        let mij_min_sq = Into::<T>::into(self.compute_min_mij().powi(2));
        let mut mij_sq =
            Into::<T>::into(self.settings.deformation.fixed.m_ij.abs().powi(2)) * mij_min_sq;

        if self.fixed_deformation.is_empty() {
            return kappas;
        }

        let cache = cache.get_cache_mut();

        // evaluate all excluded ellipsoids and pinches
        for (i, surf) in self.surfaces.iter().enumerate() {
            if surf.surface_type != SurfaceType::Pinch
                || !self.settings.deformation.fixed.dampen_on_pinch
            {
                if (self.settings.deformation.fixed.ir_handling_strategy == IRHandling::None)
                    && (surf.surface_type != SurfaceType::Ellipsoid
                        || surf.group != i
                        || !surf.exists)
                {
                    continue;
                }
            }

            let mut cut_energy = Hyperdual::from_real(Into::<T>::into(surf.shift.t));
            let mut cut_dirs = [LorentzVector::default(); MAX_LOOP];

            for &((li, pi), _, _) in &surf.id {
                let prop = &self.loop_lines[li].propagators[pi];
                cut_energy += cache.cut_info[prop.id].real_energy;

                for (cut_dir, c) in cut_dirs.iter_mut().zip(&prop.signature) {
                    *cut_dir += cache.cut_info[prop.id].momentum
                        / cache.cut_info[prop.id].real_energy.multiply_sign(*c);
                }
            }

            let mut norm_der = Hyperdual::<T, N>::zero();
            for c in &cut_dirs[..self.n_loops] {
                norm_der += c.spatial_squared_impr();
            }

            cache.ellipsoid_eval[i] = Some(cut_energy);
            cache.ellipsoid_normal_norm_eval[i] = Some(norm_der);
        }

        let mut source_scaling = Hyperdual::<T, N>::zero();
        if self.settings.deformation.fixed.source_dampening_factor > 0. {
            let aij: Hyperdual<T, N> =
                NumCast::from(self.settings.deformation.fixed.source_dampening_factor).unwrap();

            for (surf1_index, surf1) in self.surfaces.iter().enumerate() {
                for (surf2_index, surf2) in self.surfaces.iter().enumerate() {
                    if surf1.group == surf1_index
                        && surf1.surface_type == SurfaceType::Ellipsoid
                        && surf2.group == surf2_index
                        && surf2.surface_type == SurfaceType::Ellipsoid
                        && surf2_index > surf1_index
                        && surf1.exists
                        && surf2.exists
                    {
                        let eval1 = cache.ellipsoid_eval[surf1_index].unwrap().powi(2);
                        let eval2 = cache.ellipsoid_eval[surf2_index].unwrap().powi(2);
                        source_scaling +=
                            (-(eval1 + eval2) / (aij * Into::<T>::into(self.e_cm_squared))).exp();
                    }
                }
            }
        } else {
            source_scaling = Hyperdual::<T, N>::one();
        }

        for (source_index, d_lim) in self.fixed_deformation.iter().enumerate() {
            for k in &mut kappa_source {
                *k = LorentzVector::default();
            }

            for (overlap_index, d) in d_lim.deformation_per_overlap.iter().enumerate() {
                let lambda = Hyperdual::<T, N>::one();
                // TODO: index
                /*if i < self.settings.deformation.lambdas.len() {
                    lambda = NumCast::from(self.settings.deformation.lambdas[i]).unwrap();
                }*/

                if self.settings.general.debug > 2 {
                    println!(
                        "Deformation contribution from source {} and overlap {}:",
                        source_index, overlap_index
                    );
                }

                let mut s = Hyperdual::<T, N>::one();
                let mut softmin_num = Hyperdual::<T, N>::zero();
                let mut softmin_den = Hyperdual::<T, N>::zero();

                // dampen every excluded surface and pinch
                for surf_index in d.excluded_surface_indices.iter().cloned().chain(
                    self.surfaces.iter().enumerate().filter_map(|(i, su)| {
                        if self.settings.deformation.fixed.dampen_on_pinch
                            && !self.settings.deformation.fixed.dampen_on_pinch_after_lambda
                            && su.surface_type == SurfaceType::Pinch
                        {
                            Some(i)
                        } else {
                            None
                        }
                    }),
                ) {
                    let unique_ellipsoid_index = self.surfaces[surf_index]
                        .unique_ellipsoid_id
                        .unwrap_or(self.settings.deformation.fixed.m_ijs.len());
                    if unique_ellipsoid_index < self.settings.deformation.fixed.m_ijs.len() {
                        mij_sq = Into::<T>::into(
                            self.settings.deformation.fixed.m_ijs[unique_ellipsoid_index]
                                .abs()
                                .powi(2),
                        ) * mij_min_sq;
                    }
                    // t is the weighing factor that is 0 if we are on both cut_i and cut_j
                    // at the same time and goes to 1 otherwise

                    // The surface shift should never be zero at this stage.
                    let t = cache.ellipsoid_eval[surf_index].unwrap().powi(2)
                        / Into::<T>::into(self.surfaces[surf_index].shift.t.powi(2));

                    let sup = if self.settings.deformation.fixed.mode == AdditiveMode::SoftMin {
                        if self.settings.deformation.fixed.sigma.is_zero() {
                            s = s.min(t);
                            s
                        } else {
                            let e =
                                (-t / Into::<T>::into(self.settings.deformation.fixed.sigma)).exp();
                            softmin_num += t * e;
                            softmin_den += e;
                            e
                        }
                    } else {
                        let sup = t / (t + mij_sq);
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

                        if !self.settings.deformation.fixed.m_ij.abs().is_zero() {
                            s = s / (s + mij_sq);
                        }
                    }
                }

                // normalize the deformation vector per source
                if self.settings.deformation.fixed.normalize_per_source {
                    let mut normalization = Hyperdual::<T, N>::zero();
                    for ii in 0..self.n_loops {
                        let dir = loop_momenta[ii] - d.deformation_sources[ii].cast();
                        normalization += dir.spatial_squared_impr();
                    }

                    normalization = normalization.sqrt();

                    for ii in 0..self.n_loops {
                        let dir = loop_momenta[ii] - d.deformation_sources[ii].cast();
                        kappa_source[ii] += -dir / normalization * s * lambda * source_scaling;
                    }
                } else {
                    for ii in 0..self.n_loops {
                        let dir = loop_momenta[ii] - d.deformation_sources[ii].cast();
                        // the kappa returned by this function is expected to be dimensionless
                        kappa_source[ii] += -dir * s * lambda * source_scaling
                            / Hyperdual::<T, N>::from_real(Into::<T>::into(
                                self.e_cm_squared.sqrt(),
                            ));
                    }
                }

                if self.settings.general.debug > 2 {
                    for (ll, l) in loop_momenta.iter().enumerate() {
                        println!(
                            "  | k{}={:e}\n  | source={:e}\n  | suppression={:e}\n  | contribution kappa{}={:e}",
                            ll + 1,
                            l.real(),
                            &d.deformation_sources[ll],
                            s.real(),
                            ll + 1,
                            kappa_source[ll].real(),
                        );
                    }
                }
            }

            // make sure the lambda growth coming from multiple sources is under control
            if self
                .settings
                .deformation
                .fixed
                .normalisation_of_subspace_components
            {
                for ii in 0..self.n_loops {
                    kappa_source[ii] *= Hyperdual::<T, N>::from_real(Into::<T>::into(
                        1. / d_lim.deformation_per_overlap.len() as f64,
                    ));
                }
            }

            // now do the branch cut check per non-excluded loop line
            // this allows us to have a final non-zero kappa in l-space if we are on the focus in k-space
            let mut lambda_sq = Hyperdual::<T, N>::from_real(Into::<T>::into(
                self.settings
                    .deformation
                    .scaling
                    .source_branch_cut_threshold,
            ))
            .powi(2);
            for (ll_index, ll) in self.loop_lines.iter().enumerate() {
                let mut kappa_cut = LorentzVector::default();
                for (kappa, &sign) in kappa_source[..self.n_loops]
                    .iter()
                    .zip_eq(ll.signature.iter())
                {
                    kappa_cut += kappa.multiply_sign(sign);
                }

                for (prop_index, p) in ll.propagators.iter().enumerate() {
                    if d_lim.excluded_propagators.contains(&(ll_index, prop_index)) {
                        continue;
                    }
                    let mut lambda_disc_sq = cache.cut_info[p.id].spatial_and_uv_mass_sq
                        / kappa_cut.spatial_squared_impr()
                        * Hyperdual::<T, N>::from_real(Into::<T>::into(0.95));
                    if self.settings.deformation.scaling.source_branch_cut_m > 0. {
                        let branch_cut_check_m = Hyperdual::<T, N>::from_real(Into::<T>::into(
                            self.settings.deformation.scaling.source_branch_cut_m,
                        ));
                        lambda_disc_sq *= Hyperdual::<T, N>::one()
                            + (cache.cut_info[p.id].kappa_dot_mom
                                * cache.cut_info[p.id].kappa_dot_mom)
                                / (branch_cut_check_m
                                    * Into::<T>::into(self.e_cm_squared)
                                    * Into::<T>::into(self.e_cm_squared));
                    }

                    lambda_disc_sq = lambda_disc_sq.pow(Into::<T>::into(
                        self.settings.deformation.scaling.branch_cut_alpha,
                    )) * Hyperdual::<T, N>::from_real(Into::<T>::into(
                        self.settings
                            .deformation
                            .scaling
                            .source_branch_cut_multiplier,
                    ));

                    if lambda_disc_sq < lambda_sq {
                        lambda_sq = lambda_disc_sq;
                    }
                }
            }
            let lambda = if self.settings.deformation.fixed.local {
                //make the deformation localised around the surfaces
                //Set threshold for lambda
                let mut lambda0: Hyperdual<T, N> = Hyperdual::<T, N>::zero();
                let a = Into::<T>::into(self.settings.deformation.fixed.a_ijs[0]);
                let eij = &mut cache.ellipsoid_eval;
                for (surf_index, surf) in self.surfaces.iter().enumerate() {
                    if surf.surface_type != SurfaceType::Ellipsoid
                        || surf.group != surf_index
                        || !surf.exists
                    {
                        continue;
                    } else {
                        // We do not want to normalize by e_cm_squared anymore here
                        // let t = eij[surf_index].unwrap().powi(2) / Into::<T>::into(self.e_cm_squared);
                        let t = eij[surf_index].unwrap().powi(2)
                            / Into::<T>::into(self.surfaces[surf_index].shift.t.powi(2));
                        lambda0 += (-t / a).exp();
                    }
                }
                if Into::<T>::into(1.0) < lambda0.real() {
                    lambda0 = Hyperdual::<T, N>::one();
                }
                lambda_sq.sqrt() * lambda0
            } else {
                lambda_sq.sqrt()
            };

            for ii in 0..self.n_loops {
                let a = kappa_source[ii] * lambda;
                if self.settings.general.debug > 2 {
                    println!(
                        "  | lambda{}={}\n  | kappa{} contrib={:e}",
                        ii + 1,
                        lambda.real(),
                        ii + 1,
                        a.real()
                    );
                }

                kappas[ii] += a;
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{}={:e}", ii + 1, kappas[ii].real());
            }
        }

        kappas
    }

    fn normalize_on_e_surfaces<T: FloatLike, const N: usize>(
        &self,
        kappas: &mut [LorentzVector<Hyperdual<T, N>>],
        selector_m: f64,
        cache: &mut LTDCache<T>,
    ) where
        LTDCache<T>: CacheSelector<T, N>,
    {
        // Start cache
        let cache = cache.get_cache_mut();

        let mut e_surfaces_selection = Hyperdual::<T, N>::one();

        let mij_min = Into::<T>::into(self.compute_min_mij());

        // First evaluate all unique E-surfaces and multiply them into the selector
        for (i, surf) in self.surfaces.iter().enumerate() {
            if surf.surface_type != SurfaceType::Ellipsoid || i != surf.group || !surf.exists {
                continue;
            }

            if cache.ellipsoid_eval[i].is_none() {
                let mut cut_counter = 0;
                let mut cut_energy = Hyperdual::<T, N>::default();
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
                cache.ellipsoid_eval[i] = Some(
                    cut_energy
                        + Into::<T>::into(surf.shift.t)
                        + cache.cut_info[surface_prop.id]
                            .real_energy
                            .multiply_sign(surf.delta_sign),
                );
            }

            // take the square so that the number is always positive
            // We do not want to normalize by e_cm_squared anymore here
            // let t = cache.ellipsoid_eval[i].unwrap().powi(2) / Into::<T>::into(self.e_cm_squared);
            let t = cache.ellipsoid_eval[i].unwrap().powi(2)
                / Into::<T>::into(self.surfaces[i].shift.t.powi(2));

            let m = Into::<T>::into(selector_m.abs());

            e_surfaces_selection *= t / (t + m * m * mij_min * mij_min);
        }

        // Sum of \vec{kappa}^2 for the kappa of each loop
        let mut current_norm = Hyperdual::<T, N>::zero();
        for ii in 0..self.n_loops {
            current_norm += kappas[ii].spatial_squared_impr();
        }
        current_norm = (current_norm
            / Hyperdual::<T, N>::from_real(Into::<T>::into(self.n_loops as f64)))
        .sqrt();

        let normalisation = Hyperdual::<T, N>::one()
            / (e_surfaces_selection * (Hyperdual::<T, N>::one() - current_norm) + current_norm);
        //println!("normalisation_on_E_surfaces={:e}\n", normalisation.real());

        for k in kappas[..self.n_loops].iter_mut() {
            *k *= normalisation;
        }
    }

    fn deform_generic<T: FloatLike, const N: usize>(
        &self,
        loop_momenta: &[LorentzVector<Hyperdual<T, N>>],
        cache: &mut LTDCache<T>,
    ) -> ([LorentzVector<T>; MAX_LOOP], Complex<T>)
    where
        LTDCache<T>: CacheSelector<T, N>,
    {
        // compute all cut energies
        if self.settings.general.deformation_strategy != DeformationStrategy::Constant
            || self.settings.deformation.scaling.lambda > 0.
        {
            let cache = cache.get_cache_mut();

            // reset surface evaluations
            for e in &mut cache.ellipsoid_eval {
                *e = None;
            }

            for ll in &self.loop_lines {
                let mut mom = LorentzVector::default();
                for (l, &c) in izip!(loop_momenta.iter(), ll.signature.iter()) {
                    mom += l.multiply_sign(c);
                }

                // update the cut info
                for p in &ll.propagators {
                    let q: LorentzVector<Hyperdual<T, N>> = p.q.cast();
                    let m_sq = Hyperdual::<T, N>::from_real(Into::<T>::into(p.m_squared));

                    let mom_sq = (mom + q).spatial_squared_impr();
                    let cm = mom_sq + m_sq;
                    let energy = cm.sqrt();
                    cache.cut_energies[p.id] = energy;

                    let info = &mut cache.cut_info[p.id];
                    info.id = p.id;
                    info.momentum = mom + q;
                    info.real_energy = energy;
                    info.spatial_and_mass_sq = cm;
                    info.spatial_and_uv_mass_sq = mom_sq
                        + m_sq.min(Hyperdual::<T, N>::from_real(Into::<T>::into(
                            self.settings.cross_section.m_uv_sq,
                        )));
                    // TODO: these two never change and can be set at the start!
                    info.shift = q;
                    info.m_sq = m_sq;
                }
            }
        }

        // add a warning if we are close to a focus
        /*
                for (i, x) in cache.get_cache().cut_energies.iter().enumerate() {
                    if *x < Into::<T>::into(0.0000000001) {
                        println!(
                            "id={}, x={}, point={:?}",
                            i,
                            x.real(),
                            loop_momenta[..self.n_loops].iter().map(|x| x.real())
                        );
                    }
                }
        */
        let mut kappas = match self.settings.general.deformation_strategy {
            DeformationStrategy::Constant => self.deform_constant(loop_momenta),
            DeformationStrategy::Fixed => self.deform_fixed(loop_momenta, cache),
            DeformationStrategy::None => unreachable!(),
        };
        // Force normalisation of the kappa on the E-surfaces
        if self.settings.deformation.normalize_on_E_surfaces_m > 0. {
            self.normalize_on_e_surfaces(
                &mut kappas[..self.n_loops],
                self.settings.deformation.scaling.lambda,
                cache,
            );
        }

        // make sure the kappa has the right dimension by multiplying in the scale
        let scale = Hyperdual::<T, N>::from_real(
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
                    *kappa *= k_scale * Into::<T>::into(2.)
                        / (Hyperdual::<T, N>::one() + (k_scale / scale).exp());
                }
                OverallDeformationScaling::ExpDampening => {
                    *kappa *= (-k.spatial_squared_impr() / (scale * scale)).exp();
                }
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("kappa{} scaled (prelambda)={:e}", ii + 1, kappas[ii].real());
            }
        }

        let lambda = if self.settings.deformation.scaling.lambda > 0. {
            self.determine_lambda(
                &kappas[..self.n_loops],
                self.settings.deformation.scaling.lambda,
                cache,
            )
        } else {
            NumCast::from(self.settings.deformation.scaling.lambda.abs()).unwrap()
        };

        cache.overall_lambda = lambda.real();

        if self.settings.general.debug > 2 {
            println!("lambda={:e}", lambda.real());
        }

        for k in kappas[..self.n_loops].iter_mut() {
            *k *= lambda;
            k.t = Hyperdual::<T, N>::zero(); // make sure we do not have a left-over deformation
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
        cache: &mut LTDCache<T>,
    ) -> ([LorentzVector<T>; MAX_LOOP], Complex<T>) {
        if DeformationStrategy::None == self.settings.general.deformation_strategy
            || self.n_loops == 0
        {
            let r = [LorentzVector::default(); MAX_LOOP];
            return (r, Complex::new(T::one(), T::zero()));
        }

        match self.n_loops {
            1 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];

                r[0] = loop_momenta[0].map(|x| Hyperdual::<T, 4>::from_real(x));
                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                }

                return self.deform_generic(&r[..self.n_loops], cache);
            }
            2 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Hyperdual::<T, 7>::from_real(x));
                r[1] = loop_momenta[1].map(|x| Hyperdual::<T, 7>::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            _ => {}
        };

        #[cfg(feature = "higher_loops")]
        match self.n_loops {
            3 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Hyperdual::<T, 10>::from_real(x));
                r[1] = loop_momenta[1].map(|x| Hyperdual::<T, 10>::from_real(x));
                r[2] = loop_momenta[2].map(|x| Hyperdual::<T, 10>::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                    r[2][i + 1][i + 7] = T::one();
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            4 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                for j in 0..self.n_loops {
                    r[j] = loop_momenta[j].map(|x| Hyperdual::<T, 13>::from_real(x));

                    for i in 0..3 {
                        r[j][i + 1][i + 1 + j * 3] = T::one();
                    }
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            5 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                for j in 0..self.n_loops {
                    r[j] = loop_momenta[j].map(|x| Hyperdual::<T, 16>::from_real(x));

                    for i in 0..3 {
                        r[j][i + 1][i + 1 + j * 3] = T::one();
                    }
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            6 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                for j in 0..self.n_loops {
                    r[j] = loop_momenta[j].map(|x| Hyperdual::<T, 19>::from_real(x));

                    for i in 0..3 {
                        r[j][i + 1][i + 1 + j * 3] = T::one();
                    }
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            _ => {}
        }

        panic!(
            "Binding for deformation at {} loops is not implemented",
            self.n_loops
        )
    }
}
