use arrayvec::ArrayVec;
use colored::Colorize;
use dual_num::{DimName, DualN};
use float;
use itertools::Itertools;
use num::Complex;
use num_traits::ops::inv::Inv;
use num_traits::{Float, FloatConst, FromPrimitive, NumCast, One, Pow, Signed, Zero};
use partial_fractioning::{PartialFractioning, PartialFractioningMultiLoops};
use rand::rngs::StdRng;
use rand::seq::IteratorRandom;
use rand::{Rng, SeedableRng};
use std::cmp::Ordering;
use topologies::{
    CacheSelector, Cut, CutList, LTDCache, LTDNumerator, LoopLine, SOCPProblem, Surface,
    SurfaceType, Topology,
};
use utils::Signum;
use vector::LorentzVector;
use {
    AdditiveMode, DeformationStrategy, ExpansionCheckStrategy, FloatLike, IRHandling,
    OverallDeformationScaling, ParameterizationMapping, ParameterizationMode, PoleCheckStrategy,
    Settings, MAX_LOOP,
};

use utils;

type Dual4<T> = DualN<T, dual_num::U4>;
type Dual7<T> = DualN<T, dual_num::U7>;
#[cfg(feature = "higher_loops")]
type Dual10<T> = DualN<T, dual_num::U10>;
#[cfg(feature = "higher_loops")]
type Dual13<T> = DualN<T, dual_num::U13>;
#[cfg(feature = "higher_loops")]
type Dual16<T> = DualN<T, dual_num::U16>;
#[cfg(feature = "higher_loops")]
type Dual19<T> = DualN<T, dual_num::U19>;

impl LoopLine {
    /// Return the inverse of the evaluated loop line
    fn evaluate<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<num::Complex<T>>],
        cut: &Cut,
        topo: &Topology,
        cache: &mut LTDCache<T>,
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
                    let mut r = cache.complex_cut_energies[p.id] * Into::<T>::into(2.);
                    cache.propagators_eval[p.id] = r;
                    r = utils::powi(r, p.power);
                    res *= r;

                    if topo.settings.general.debug > 3 {
                        println!("  | prop x {}={}", i, r);
                    }
                }
                _ => {
                    // multiply dual propagator
                    let mut r = utils::powi(e + Into::<T>::into(p.q.t), 2)
                        - cache.complex_prop_spatial[p.id];
                    cache.propagators_eval[p.id] = r;
                    r = utils::powi(r, p.power);

                    if topo.settings.general.debug > 3 {
                        println!("  | prop   {}={}", i, r);
                    }

                    if !r.re.is_finite()
                        || !r.im.is_finite()
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
    /// If `external_momenta_set` is true, the existence conditions
    /// for ellipsoids and e_cm are determined as well.
    pub fn process(&mut self, external_momenta_set: bool) {
        // construct the numerator
        self.numerator = if self.numerator_tensor_coefficients.len() == 0
            && self.numerator_tensor_coefficients_sparse.len() == 0
        {
            LTDNumerator::one(self.n_loops)
        } else {
            if self.numerator_tensor_coefficients_sparse.len() > 0 {
                LTDNumerator::from_sparse(self.n_loops, &self.numerator_tensor_coefficients_sparse)
            } else {
                LTDNumerator::new(
                    self.n_loops,
                    &self
                        .numerator_tensor_coefficients
                        .iter()
                        .map(|x| Complex::new(x.0, x.1))
                        .collect::<Vec<_>>(),
                )
            }
        };
        // Prepare the partial fractioning map at one loop if the threshold is set to a positive number
        if self.n_loops == 1
            && self.settings.general.partial_fractioning_multiloop
            && self.settings.general.partial_fractioning_threshold > 0.0
        {
            let num_propagators_deg_1l: usize = self
                .loop_lines
                .iter()
                .filter(|x| x.signature == &[1])
                .map(|x| x.propagators.iter().map(|p| p.power).sum::<usize>())
                .sum();

            self.partial_fractioning = PartialFractioning::new(num_propagators_deg_1l, 10000);
        }
        if self.n_loops > 0 && self.settings.general.partial_fractioning_threshold > 0.0 {
            self.partial_fractioning_multiloops =
                PartialFractioningMultiLoops::new(&self.loop_lines, 100000);
        }
        // set the identity rotation matrix
        self.rotation_matrix = [
            [float::one(), float::zero(), float::zero()],
            [float::zero(), float::one(), float::zero()],
            [float::zero(), float::zero(), float::one()],
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
                self.e_cm_squared = self.external_kinematics[0].spatial_squared_impr();
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
                let mut cut_shift: Vec<LorentzVector<float>> = vec![]; // qs of cut
                let mut cut_mass = vec![];
                for (cut, ll) in cut_option.iter().zip_eq(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut
                    {
                        cut_shift.push(ll.propagators[*cut_prop_index].q.cast());
                        cut_mass.push(
                            float::from_f64(ll.propagators[*cut_prop_index].m_squared.sqrt())
                                .unwrap(),
                        );
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

                        let mut surface_shift: LorentzVector<float> = onshell_prop.q.cast();
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

                        let mut cut_mass_sum = float::zero();
                        for (ss, mass) in surface_signs.iter().zip_eq(cut_mass.iter()) {
                            cut_mass_sum += mass.multiply_sign(*ss);
                        }

                        let surface_mass = float::from_f64(onshell_prop.m_squared.sqrt()).unwrap();

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
                                        >= Into::<float>::into(-1e-13 * self.e_cm_squared)
                                        && surface_shift.t.multiply_sign(delta_sign) < float::zero()
                                    {
                                        if surface_shift.square()
                                            - (cut_mass_sum.abs() + surface_mass).powi(2)
                                            < Into::<float>::into(1e-10 * self.e_cm_squared)
                                        {
                                            if surface_mass.is_zero() {
                                                is_pinch = true;
                                            } else {
                                                // we do not consider pointlike ellipsoids to exist
                                                exists = false;
                                            }
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
                                            || neg_surface_signs_count == 1 && eval < float::zero()
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
            let mut surface_shift: LorentzVector<float> = self.loop_lines[s.onshell_ll_index]
                .propagators[s.onshell_prop_index]
                .q
                .cast();

            let mut mass_sum: float = self.loop_lines[s.onshell_ll_index].propagators
                [s.onshell_prop_index]
                .m_squared
                .sqrt()
                .into();

            let mut cut_index = 0;
            for (cut, ll) in s.cut.iter().zip_eq(self.loop_lines.iter()) {
                if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut {
                    mass_sum +=
                        Into::<float>::into(ll.propagators[*cut_prop_index].m_squared.sqrt());
                    surface_shift -= ll.propagators[*cut_prop_index]
                        .q
                        .cast()
                        .multiply_sign(s.sig_ll_in_cb[cut_index]);
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

            if size >= Into::<float>::into(-1e-13 * self.e_cm_squared)
                && surface_shift.t.multiply_sign(s.delta_sign) < float::zero()
            {
                if self.settings.general.debug > 0 {
                    if size > 1e-8 * self.e_cm_squared
                        && size < 1e-6 * self.e_cm_squared
                        && !mass_sum.is_zero()
                    {
                        println!("Small ellipsoid {} detected: {}", s.group, size);
                    }
                }

                if size < Into::<float>::into(1e-8 * self.e_cm_squared) {
                    if mass_sum.is_zero() {
                        is_pinch = true;
                    } else {
                        // we do not consider pointlike ellipsoids to exist
                        exists = false;
                    }
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
                            panic!(
                                "Deformation source {:?} is not on the inside of surface {}: {}",
                                d.deformation_sources, surf_index, r.re
                            );
                        }
                    }
                }

                // now test if the source is equal to the shift of the propagators
                for (ll_index, prop_index) in &d_lim.excluded_propagators {
                    // determine the sum of sources using the signature
                    let mut source_sum: LorentzVector<float> = LorentzVector::default();
                    for (sign, source) in self.loop_lines[*ll_index]
                        .signature
                        .iter()
                        .zip_eq(&d.deformation_sources)
                    {
                        source_sum += source.multiply_sign(*sign);
                    }

                    let diff: LorentzVector<float> =
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
        x: &[f64],
        e_cm_squared: f64,
        loop_index: usize,
        settings: &Settings,
    ) -> ([T; 3], T) {
        let e_cm = Into::<T>::into(e_cm_squared).sqrt()
            * Into::<T>::into(settings.parameterization.shifts[loop_index].0);
        let mut l_space = [T::zero(); 3];
        let mut jac = T::one();

        // rescale the input to the desired range
        let mut x_r = [0.; 3];
        for (xd, xi, &(lo, hi)) in izip!(
            &mut x_r,
            x,
            &settings.parameterization.input_rescaling[loop_index]
        ) {
            *xd = lo + xi * (hi - lo);
            jac *= Into::<T>::into(hi - lo);
        }

        match settings.parameterization.mode {
            ParameterizationMode::Cartesian => match settings.parameterization.mapping {
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
                let radius = match settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        // r = e_cm * ln(1 + b*x/(1-x))
                        let x = Into::<T>::into(x_r[0]);
                        let b = Into::<T>::into(settings.parameterization.b);
                        let radius = e_cm * (T::one() + b * x / (T::one() - x)).ln();
                        jac *= e_cm * b / (T::one() - x) / (T::one() + x * (b - T::one()));

                        radius
                    }
                    ParameterizationMapping::Linear => {
                        // r = e_cm * b * x/(1-x)
                        let b = Into::<T>::into(settings.parameterization.b);
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
        l_space[0] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].1);
        l_space[1] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].2);
        l_space[2] += e_cm * Into::<T>::into(settings.parameterization.shifts[loop_index].3);

        (l_space, jac)
    }

    pub fn inv_parametrize<T: FloatLike>(
        mom: &LorentzVector<T>,
        e_cm_squared: f64,
        loop_index: usize,
        settings: &Settings,
    ) -> ([T; 3], T) {
        if settings.parameterization.mode != ParameterizationMode::Spherical {
            panic!("Inverse mapping is only implemented for spherical coordinates");
        }

        let mut jac = T::one();
        let e_cm = Into::<T>::into(e_cm_squared).sqrt()
            * Into::<T>::into(settings.parameterization.shifts[loop_index].0);

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
    ) -> DualN<T, U>
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
                    let c = DualN::from_real(Into::<T>::into(self.get_expansion_threshold()));

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
                        * DualN::from_real(Into::<T>::into(0.95));

                    if self.settings.deformation.scaling.branch_cut_m > 0. {
                        let branch_cut_check_m = DualN::from_real(Into::<T>::into(
                            self.settings.deformation.scaling.branch_cut_m,
                        ));
                        lambda_disc_sq *= DualN::one()
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
                    let mut a = DualN::from_real(T::zero());
                    let mut b = DualN::from_real(T::zero());
                    let mut c = DualN::from_real(T::zero());
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

                            let mut lambda_plus = (-Complex::new(DualN::zero(), bb)
                                + (Complex::new(
                                    -bb * bb - aa * cc * Into::<T>::into(4.),
                                    DualN::zero(),
                                ))
                                .sqrt())
                                / (aa * Into::<T>::into(2.));
                            let mut lambda_min = (-Complex::new(DualN::zero(), bb)
                                - (Complex::new(
                                    -bb * bb - aa * cc * Into::<T>::into(4.),
                                    DualN::zero(),
                                ))
                                .sqrt())
                                / (aa * Into::<T>::into(2.));

                            if aa.real().is_zero() {
                                lambda_plus = Complex::new(DualN::zero(), cc / bb);
                                lambda_min = Complex::new(DualN::zero(), cc / bb);
                            }

                            // evaluate the surface with lambda
                            let mut prop_lambda_sq = lambda_sq;
                            for lambda in &[lambda_plus, lambda_min] {
                                let def = kappas[0].real().to_complex(false)
                                    * Complex::new(lambda.re.real(), lambda.im.real());
                                let surf = ((cut_info.momentum.real().to_complex(true) + def)
                                    .spatial_squared()
                                    + cut_info.mass.real()) // note, this is mass^2
                                .sqrt()
                                    + ((on_shell_info.momentum.real().to_complex(true) + def)
                                        .spatial_squared()
                                        + on_shell_info.mass.real())
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

        // dampen every pinch
        if self.settings.deformation.fixed.dampen_on_pinch
            && self.settings.deformation.fixed.dampen_on_pinch_after_lambda
            && !self.fixed_deformation.is_empty()
        {
            let mij_min_sq = Into::<T>::into(self.compute_min_mij().powi(2));
            let mij_sq =
                Into::<T>::into(self.settings.deformation.fixed.m_ij.abs().powi(2)) * mij_min_sq;
            for (surf_index, s) in self.surfaces.iter().enumerate() {
                if s.surface_type != SurfaceType::Pinch {
                    continue;
                }

                let e = cache.ellipsoid_eval[surf_index].unwrap().powi(2);
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

                let sup = t / (t + mij_sq);

                if sup * sup < lambda_sq {
                    lambda_sq = sup * sup;
                }
            }
        }

        // Setup a veto region at the intersection of pinched and non-pinched E-surfaces.
        if self.settings.deformation.fixed.ir_handling_strategy != IRHandling::None {
            let mut min_ellipse = DualN::from_real(Into::<T>::into(1e99));
            let mut min_pinch = DualN::from_real(Into::<T>::into(1e99));
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
                    .unwrap_or(DualN::from_real(Into::<T>::into(1e99)))
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

            let mut min_e = DualN::from_real(Into::<T>::into(1e99));
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

            if ( self.settings.deformation.fixed.ir_handling_strategy==IRHandling::DismissDeformation && 
                 self.settings.deformation.fixed.ir_interpolation_length > 0.0) {
                return ((
                            (
                                ir_proximity/Into::<T>::into(self.settings.deformation.fixed.ir_threshold.powi(2))
                            )-Into::<T>::into(1.0)
                    )/Into::<T>::into(self.settings.deformation.fixed.ir_interpolation_length))
                    .min(DualN::from_real(Into::<T>::into(0.0))).max(DualN::from_real(Into::<T>::into(1.0)));
            } else {
                if ir_proximity < Into::<T>::into(self.settings.deformation.fixed.ir_threshold.powi(2))
                {
                    match self.settings.deformation.fixed.ir_handling_strategy {
                        IRHandling::DismissPoint => {
                            return DualN::from_real(T::nan()); // TODO: improve into a proper escalation
                        }
                        IRHandling::DismissDeformation => {
                            return DualN::from_real(Into::<T>::into(0.));
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

    /// Compute the normal vector for each ellipsoid and evaluate each ellipsoid.
    fn compute_ellipsoid_deformation_vector<U: DimName, T: FloatLike>(
        &self,
        normalize: bool,
        cache: &mut LTDCache<T>,
    ) where
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
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            // only deform the set of different ellipsoids
            if surf.surface_type != SurfaceType::Ellipsoid
                || surf.group != surf_index
                || !surf.exists
            {
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

            inv_surf_prop[surf_index] = Some(
                cut_energy
                    + Into::<T>::into(surf.shift.t)
                    + cache.cut_info[surface_prop.id]
                        .real_energy
                        .multiply_sign(surf.delta_sign),
            );

            let inv_normalization = if normalize {
                let mut normalization = DualN::zero();

                for d in deform_dirs[..self.n_loops].iter() {
                    normalization += d.spatial_squared_impr();
                }

                normalization.sqrt().inv()
            } else {
                DualN::one()
            };

            for (loop_index, dir) in deform_dirs[..self.n_loops].iter().enumerate() {
                // note the sign
                kappa_surf[surf_index * self.n_loops + loop_index] = -dir * inv_normalization;
            }
        }
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
        self.compute_ellipsoid_deformation_vector(true, cache);

        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let cache = cache.get_cache_mut();
        let kappa_surf = &mut cache.deform_dirs;
        let inv_surf_prop = &mut cache.ellipsoid_eval;

        // now combine the kappas from the surface using the chosen strategy
        let mut aij: DualN<T, U> = NumCast::from(self.settings.deformation.additive.a_ij).unwrap();
        let mut lambda: DualN<T, U> = DualN::one();

        for (i, &inv) in inv_surf_prop.iter().enumerate() {
            // not a unique ellipsoid
            if inv.is_none() {
                continue;
            }
            let inv = inv.unwrap();

            if i < self.settings.deformation.additive.a_ijs.len() {
                aij = NumCast::from(self.settings.deformation.additive.a_ijs[i]).unwrap();
            }

            if i < self.settings.deformation.lambdas.len() {
                lambda = NumCast::from(self.settings.deformation.lambdas[i]).unwrap();
            }

            let dampening = match self.settings.deformation.additive.mode {
                AdditiveMode::Exponential => {
                    (-inv * inv / (aij * Into::<T>::into(self.e_cm_squared))).exp()
                }
                AdditiveMode::Hyperbolic => {
                    let t = inv * inv / Into::<T>::into(self.e_cm_squared);
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
                        i,
                        dir,
                        dampening.real()
                    );
                }

                *kappa += dir * dampening * lambda;
            }
        }

        kappas
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
    fn deform_constant<U: DimName, T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
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

        let mut normalization = DualN::zero();
        for i in 0..self.n_loops {
            normalization += kappas[i].spatial_squared_impr();
        }
        normalization = normalization.sqrt().inv();

        for i in 0..self.n_loops {
            kappas[i] *= normalization;
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
        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let mut kappa_source = [LorentzVector::default(); MAX_LOOP];
        let mij_min_sq = Into::<T>::into(self.compute_min_mij().powi(2));
        let mut mij_sq =
            Into::<T>::into(self.settings.deformation.fixed.m_ij.abs().powi(2)) * mij_min_sq;

        if self.fixed_deformation.is_empty() {
            return kappas;
        }

        if self.settings.deformation.fixed.include_normal_source {
            self.compute_ellipsoid_deformation_vector(
                self.settings.deformation.fixed.normalize_per_source,
                cache,
            );
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
                        || !surf.exists
                        || (!self.all_excluded_surfaces[i]
                            && !self.settings.deformation.fixed.local)
                        || cache.ellipsoid_eval[i].is_some())
                {
                    continue;
                }
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
            let eval = Some(
                cut_energy
                    + Into::<T>::into(surf.shift.t)
                    + cache.cut_info[surface_prop.id]
                        .real_energy
                        .multiply_sign(surf.delta_sign),
            );

            cache.ellipsoid_eval[i] = eval;
        }

        let mut source_scaling = DualN::zero();
        if self.settings.deformation.fixed.source_dampening_factor > 0. {
            let aij: DualN<T, U> =
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
            source_scaling = DualN::one();
        }

        for (source_index, d_lim) in self.fixed_deformation.iter().enumerate() {
            for k in &mut kappa_source {
                *k = LorentzVector::default();
            }

            // add the normal for each surface to the deformation
            // TODO: call the additive function instead
            if self.settings.deformation.fixed.include_normal_source
                && d_lim.excluded_propagators.len() == 0
            {
                let kappa_surf = &mut cache.deform_dirs;
                let inv_surf_prop = &mut cache.ellipsoid_eval;

                // now combine the kappas from the surface using the chosen strategy
                let mut aij: DualN<T, U> =
                    NumCast::from(self.settings.deformation.additive.a_ij).unwrap();
                let mut lambda: DualN<T, U> = DualN::one();

                for (i, &inv) in inv_surf_prop.iter().enumerate() {
                    if inv.is_none() {
                        continue;
                    }
                    let inv = inv.unwrap();

                    let unique_ellipsoid_index = self.surfaces[i].unique_ellipsoid_id.unwrap();
                    if unique_ellipsoid_index < self.settings.deformation.additive.a_ijs.len() {
                        aij = NumCast::from(
                            self.settings.deformation.additive.a_ijs[unique_ellipsoid_index],
                        )
                        .unwrap();
                    }

                    if unique_ellipsoid_index < self.settings.deformation.lambdas.len() {
                        lambda = NumCast::from(
                            self.settings.deformation.lambdas[unique_ellipsoid_index],
                        )
                        .unwrap();
                    }

                    let mut dampening = match self.settings.deformation.additive.mode {
                        AdditiveMode::Exponential => {
                            (-inv * inv / (aij * Into::<T>::into(self.e_cm_squared))).exp()
                        }
                        AdditiveMode::Hyperbolic => {
                            let t = inv * inv / Into::<T>::into(self.e_cm_squared);
                            t / (t + aij)
                        }
                        AdditiveMode::Unity => DualN::one(),
                        AdditiveMode::SoftMin => unimplemented!(),
                    };

                    // now exclude all other surfaces
                    for (surf_index, surf) in self.surfaces.iter().enumerate() {
                        if i == surf_index
                            || (surf.surface_type != SurfaceType::Ellipsoid
                                && (surf.surface_type != SurfaceType::Pinch
                                    || !self.settings.deformation.fixed.dampen_on_pinch))
                            || surf_index != surf.group
                            || !surf.exists
                        {
                            continue;
                        }

                        let t = inv_surf_prop[surf_index].unwrap().powi(2)
                            / Into::<T>::into(self.surfaces[surf_index].shift.t.powi(2));

                        let unique_ellipsoid_index = surf
                            .unique_ellipsoid_id
                            .unwrap_or(self.settings.deformation.fixed.m_ijs.len());
                        if unique_ellipsoid_index < self.settings.deformation.fixed.m_ijs.len() {
                            dampening *= t
                                / (t + Into::<T>::into(
                                    self.settings.deformation.fixed.m_ijs[unique_ellipsoid_index]
                                        .powi(2),
                                ));
                        } else {
                            dampening *= t / (t + mij_sq);
                        }
                    }

                    for (loop_index, kappa) in kappa_source[..self.n_loops].iter_mut().enumerate() {
                        let dir = kappa_surf[i * self.n_loops + loop_index];
                        *kappa += dir * dampening * lambda;
                    }
                }
            }

            for (overlap_index, d) in d_lim.deformation_per_overlap.iter().enumerate() {
                let lambda = DualN::one();
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

                let mut s = DualN::one();
                let mut softmin_num = DualN::zero();
                let mut softmin_den = DualN::zero();

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
                if self.settings.general.debug > 2 {
                    println!(
                        "  | k={:e}\n  | dirs={:e}\n  | suppression={:e}\n  | contribution={:e}",
                        loop_momenta[0].real(),
                        &d.deformation_sources[0],
                        s.real(),
                        &d.deformation_sources[0].cast() * s.real(),
                    );
                }

                // normalize the deformation vector per source
                if self.settings.deformation.fixed.normalize_per_source {
                    let mut normalization = DualN::zero();
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
                            / DualN::from_real(Into::<T>::into(self.e_cm_squared.sqrt()));
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
                    kappa_source[ii] *= DualN::from_real(Into::<T>::into(
                        1. / d_lim.deformation_per_overlap.len() as f64,
                    ));
                }
            }

            // now do the branch cut check per non-excluded loop line
            // this allows us to have a final non-zero kappa in l-space if we are on the focus in k-space
            let mut lambda_sq = DualN::from_real(Into::<T>::into(
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
                        * DualN::from_real(Into::<T>::into(0.95));
                    if self.settings.deformation.scaling.source_branch_cut_m > 0. {
                        let branch_cut_check_m = DualN::from_real(Into::<T>::into(
                            self.settings.deformation.scaling.source_branch_cut_m,
                        ));
                        lambda_disc_sq *= DualN::one()
                            + (cache.cut_info[p.id].kappa_dot_mom
                                * cache.cut_info[p.id].kappa_dot_mom)
                                / (branch_cut_check_m
                                    * Into::<T>::into(self.e_cm_squared)
                                    * Into::<T>::into(self.e_cm_squared));
                    }

                    lambda_disc_sq = lambda_disc_sq.pow(Into::<T>::into(
                        self.settings.deformation.scaling.branch_cut_alpha,
                    )) * DualN::from_real(Into::<T>::into(
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
                let mut lambda0: DualN<T, U> = DualN::zero();
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
                    lambda0 = DualN::one();
                }
                lambda_sq.sqrt() * lambda0
            } else {
                lambda_sq.sqrt()
            };

            for ii in 0..self.n_loops {
                kappas[ii] += kappa_source[ii] * lambda;
            }
        }

        if self
            .settings
            .deformation
            .fixed
            .normalisation_per_number_of_sources
        {
            for ii in 0..self.n_loops {
                kappas[ii] *=
                    DualN::from_real(Into::<T>::into(1. / self.fixed_deformation.len() as f64));
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{}={:e}", ii + 1, kappas[ii].real());
            }
        }

        kappas
    }

    fn normalize_on_e_surfaces<U: DimName, T: FloatLike>(
        &self,
        kappas: &mut [LorentzVector<DualN<T, U>>],
        selector_m: f64,
        cache: &mut LTDCache<T>,
    ) where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        // Start cache
        let cache = cache.get_cache_mut();

        let mut e_surfaces_selection = DualN::one();

        let mij_min = Into::<T>::into(self.compute_min_mij());

        // First evaluate all unique E-surfaces and multiply them into the selector
        for (i, surf) in self.surfaces.iter().enumerate() {
            if surf.surface_type != SurfaceType::Ellipsoid || i != surf.group || !surf.exists {
                continue;
            }

            if cache.ellipsoid_eval[i].is_none() {
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
        let mut current_norm = DualN::zero();
        for ii in 0..self.n_loops {
            current_norm += kappas[ii].spatial_squared_impr();
        }
        current_norm =
            (current_norm / DualN::from_real(Into::<T>::into(self.n_loops as f64))).sqrt();

        let normalisation =
            DualN::one() / (e_surfaces_selection * (DualN::one() - current_norm) + current_norm);
        //println!("normalisation_on_E_surfaces={:e}\n", normalisation.real());

        for k in kappas[..self.n_loops].iter_mut() {
            *k *= normalisation;
        }
    }

    fn deform_generic<U: DimName, T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
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
                    let q: LorentzVector<DualN<T, U>> = p.q.cast();
                    let mass = DualN::from_real(Into::<T>::into(p.m_squared));

                    let mom_sq = (mom + q).spatial_squared_impr();
                    let cm = mom_sq + mass;
                    let energy = cm.sqrt();
                    cache.cut_energies[p.id] = energy;

                    let info = &mut cache.cut_info[p.id];
                    info.id = p.id;
                    info.momentum = mom + q;
                    info.real_energy = energy;
                    info.spatial_and_mass_sq = cm;
                    info.spatial_and_uv_mass_sq = mom_sq
                        + mass.min(DualN::from_real(Into::<T>::into(
                            self.settings.cross_section.m_uv_sq,
                        )));
                    // TODO: these two never change and can be set at the start!
                    info.shift = q;
                    info.mass = mass;
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
            DeformationStrategy::Additive => self.deform_ellipsoids(cache),
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
                OverallDeformationScaling::ExpDampening => {
                    *kappa *= (-k.spatial_squared_impr() / (scale * scale)).exp();
                }
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("kappa{} scaled={:e}", ii + 1, kappas[ii].real());
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

                r[0] = loop_momenta[0].map(|x| Dual4::from_real(x));
                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                }

                return self.deform_generic(&r[..self.n_loops], cache);
            }
            2 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual7::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual7::from_real(x));

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
                r[0] = loop_momenta[0].map(|x| Dual10::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual10::from_real(x));
                r[2] = loop_momenta[2].map(|x| Dual10::from_real(x));

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
                    r[j] = loop_momenta[j].map(|x| Dual13::from_real(x));

                    for i in 0..3 {
                        r[j][i + 1][i + 1 + j * 3] = T::one();
                    }
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            5 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                for j in 0..self.n_loops {
                    r[j] = loop_momenta[j].map(|x| Dual16::from_real(x));

                    for i in 0..3 {
                        r[j][i + 1][i + 1 + j * 3] = T::one();
                    }
                }
                return self.deform_generic(&r[..self.n_loops], cache);
            }
            6 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                for j in 0..self.n_loops {
                    r[j] = loop_momenta[j].map(|x| Dual19::from_real(x));

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

    /// Set the energy component of the loop momenta according to
    /// `cut`. It takes the cut energies from the cache.
    pub fn set_loop_momentum_energies<T: FloatLike>(
        &self,
        k_def: &mut [LorentzVector<Complex<T>>],
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

    /// Evaluate a higher order cut by differentiating all non-cut propagators and the numerator.
    /// The power of the cut propagators should be 0.
    pub fn evaluate_higher_order_cut_2<T: FloatLike>(
        &self,
        cut: &Vec<Cut>,
        cut_propagators: &[usize],
        derivative_map: &mut [usize],   // in the cut basis
        numerator_powers: &mut [usize], // in the cut basis
        mat: &Vec<i8>,
        index: usize,
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        if index == derivative_map.len() {
            let mut result = Complex::one();

            // we are done, so evaluate the numerator and denominators
            for (p_e, p_pow) in cache
                .propagators_eval
                .iter()
                .zip_eq(cache.propagator_powers.iter())
            {
                result *= utils::powi(*p_e, *p_pow);
            }
            result = result.inv();

            // evaluate the numerator
            for (cut_prop, num) in cut_propagators.iter().zip_eq(numerator_powers.iter()) {
                if *num > 0 {
                    result *= match cut[self.propagator_id_to_ll_id[*cut_prop].0] {
                        Cut::PositiveCut(_) => cache.complex_cut_energies[*cut_prop],
                        Cut::NegativeCut(_) => -cache.complex_cut_energies[*cut_prop],
                        Cut::NoCut => panic!("Trying to multiply the energy for uncut propagator!"),
                    }
                    .powi(*num as i32);
                }
            }
            return result;
        }

        if derivative_map[index] == 0 {
            // go to the next cut
            return self.evaluate_higher_order_cut_2(
                cut,
                cut_propagators,
                derivative_map,
                numerator_powers,
                mat,
                index + 1,
                cache,
            );
        }

        let mut result = Complex::zero();

        derivative_map[index] -= 1;
        // Check if numerator is differentiable
        if numerator_powers[index] > 0 {
            numerator_powers[index] -= 1;
            result += self.evaluate_higher_order_cut_2(
                cut,
                cut_propagators,
                derivative_map,
                numerator_powers,
                mat,
                index,
                cache,
            ) * Into::<T>::into(numerator_powers[index] as f64 + 1.);
            numerator_powers[index] += 1;
        }

        // Take derivative of denominator
        for prop_id in 0..cache.propagator_powers.len() {
            // skip denominators that are effectively absent (such as cut propagators)
            if cache.propagator_powers[prop_id] == 0 {
                continue;
            }

            let (ll_index, p_index) = self.propagator_id_to_ll_id[prop_id];

            let signature_in_ll = &self.loop_lines[ll_index].signature;

            // convert the signature to the cut basis
            // TODO: cache
            let mut shift = Complex::new(
                Into::<T>::into(self.loop_lines[ll_index].propagators[p_index].q.t),
                T::zero(),
            );
            let mut signature_in_cb = [0; MAX_LOOP];
            for (sig_ll, row) in signature_in_ll
                .iter()
                .zip_eq(mat.chunks_exact(self.n_loops))
            {
                for (col, sig) in row
                    .iter()
                    .zip_eq(signature_in_cb[..self.n_loops].iter_mut())
                {
                    *sig += col * sig_ll;
                }
            }

            // compute the shifts of the propagator in the cut basis
            for (cp, scp) in cut_propagators.iter().zip(&signature_in_cb[..self.n_loops]) {
                let (cp_ll_index, cp_prop_index) = self.propagator_id_to_ll_id[*cp];
                shift -=
                    Into::<T>::into(self.loop_lines[cp_ll_index].propagators[cp_prop_index].q.t)
                        .multiply_sign(*scp);
            }

            let d_coeff = Into::<T>::into(cache.propagator_powers[prop_id] as f64)
                * Into::<T>::into(-2. * signature_in_cb[index] as f64);
            cache.propagator_powers[prop_id] += 1;

            for (i, &sig) in signature_in_cb.iter().enumerate() {
                if sig != 0 {
                    if derivative_map[i] > 0 {
                        // we need to keep track of the power of this variable for further differentiation
                        numerator_powers[i] += 1;
                        result += self.evaluate_higher_order_cut_2(
                            cut,
                            cut_propagators,
                            derivative_map,
                            numerator_powers,
                            mat,
                            index,
                            cache,
                        ) * d_coeff.multiply_sign(sig);
                        numerator_powers[i] -= 1;
                    } else {
                        // the variable does not need to be differentiated, so accumulate it in the shift
                        shift += match cut[i] {
                            Cut::PositiveCut(_) | Cut::NoCut => {
                                cache.complex_cut_energies[cut_propagators[i]].multiply_sign(sig)
                            }
                            Cut::NegativeCut(_) => {
                                cache.complex_cut_energies[cut_propagators[i]].multiply_sign(-sig)
                            }
                        };
                    }
                }
            }

            // evaluate the shift part
            result += self.evaluate_higher_order_cut_2(
                cut,
                cut_propagators,
                derivative_map,
                numerator_powers,
                mat,
                index,
                cache,
            ) * shift
                * d_coeff;

            cache.propagator_powers[prop_id] -= 1;
        }

        derivative_map[index] += 1;
        result
    }

    /// Evaluate a higher order cut. The first stage is this recursive function that
    /// goes through the combinatorics of deriving the cut propagators.
    /// The second recursive function is `evaluate_higher_order_cut_2`.
    ///     - cut_propagators: which propagators ids are being cut
    ///     - derivative_map: order of derivatives expressed in cut basis
    ///     - numerator_powes: powers of the monomial in cut basis
    ///     - mat : if chuncked is a matrix with the cut ll signatures as rows
    pub fn evaluate_higher_order_cut<T: FloatLike>(
        &self,
        cut: &Vec<Cut>,
        cut_propagators: &[usize],
        derivative_map: &mut [usize],
        numerator_powers: &mut [usize],
        mat: &Vec<i8>,
        index: usize,
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        if index == derivative_map.len() {
            return self.evaluate_higher_order_cut_2(
                cut,
                cut_propagators,
                derivative_map,
                numerator_powers,
                mat,
                0,
                cache,
            );
        }
        // Compute the energy corresponding to the current cut propagator
        //let energy = cache.propagators_eval[cut_propagators[index]];
        let energy = match cut[self.propagator_id_to_ll_id[cut_propagators[index]].0] {
            Cut::PositiveCut(_) => {
                cache.complex_cut_energies[cut_propagators[index]] * Into::<T>::into(2.0)
            }
            Cut::NegativeCut(_) => {
                -cache.complex_cut_energies[cut_propagators[index]] * Into::<T>::into(2.0)
            }
            Cut::NoCut => panic!("Trying to multiply the energy for uncut propagator!"),
        };

        // for every cut, we first take out the derivatives wrt the cut propagators
        // as they transform slightly differently
        let mut result = Complex::default();
        let mut coeff = T::one();
        let mut norm = T::one();
        let n = derivative_map[index];
        for i in (0..=n).rev() {
            derivative_map[index] = i;

            result += self.evaluate_higher_order_cut(
                cut,
                cut_propagators,
                derivative_map,
                numerator_powers,
                mat,
                index + 1,
                cache,
            ) * Float::powi(-T::one(), (n - i) as i32)
                * coeff
                / utils::powi(energy, 1 + 2 * n - i);

            if i != 0 {
                coeff *=
                    T::from_usize(i * (2 * n - i + 1)).unwrap() / T::from_usize(n - i + 1).unwrap();
                norm *= T::from_usize(i).unwrap();
            }
        }
        derivative_map[index] = n;
        result / norm
    }

    pub fn evaluate_cut<T: FloatLike>(
        &self,
        k_def: &mut [LorentzVector<Complex<T>>],
        numerator: &LTDNumerator,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
        overwrite_propagator_powers: bool,
        num_id: usize,
    ) -> Result<Complex<T>, &'static str> {
        let mut cut_propagators = [0; MAX_LOOP];
        let mut derivative_map = [0; MAX_LOOP];
        let mut numerator_powers = [0; MAX_LOOP];

        // check if we need to take derivatives
        let mut index = 0;
        for (c, ll) in cut.iter().zip(&self.loop_lines) {
            if overwrite_propagator_powers {
                for p in &ll.propagators {
                    cache.propagator_powers[p.id] = p.power;
                }
            }
            if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = c {
                // Using the value form the cache to allow setting externally
                // the powers of the propagators
                let p_id = ll.propagators[*i].id;
                if cache.propagator_powers[p_id] == 0 || ll.propagators[*i].power == 0 {
                    return Ok(Complex::zero());
                }
                derivative_map[index] = cache.propagator_powers[p_id] - 1;
                cache.propagator_powers[p_id] = 0; // do not consider these linear propagators in the recursion
                cut_propagators[index] = p_id;

                index += 1;
            }
        }

        //dbg!(&cache.propagator_powers, &derivative_map);

        // Set energy on the cut
        self.set_loop_momentum_energies(k_def, cut, mat, cache);

        // Prepare numerator
        //numerator.evaluate_reduced_in_lb(k_def, cache); //NOTE: Only necessary when k_vec is changed
        let r = if derivative_map.iter().all(|x| *x == 0) {
            // we are in the case where we do not need derivatives
            let mut r = Complex::one();

            if overwrite_propagator_powers {
                for (i, (ll_cut, ll)) in cut.iter().zip_eq(self.loop_lines.iter()).enumerate() {
                    // get the loop line result from the cache if possible
                    r *= match ll_cut {
                        Cut::PositiveCut(j) => cache.complex_loop_line_eval[i][*j][0],
                        Cut::NegativeCut(j) => cache.complex_loop_line_eval[i][*j][1],
                        _ => ll.evaluate(&k_def, ll_cut, &self, cache)?,
                    };
                }
            } else {
                // Cache all propagator evaluations by evaluating the loop line
                for (ll_cut, ll) in cut.iter().zip_eq(self.loop_lines.iter()) {
                    ll.evaluate(&k_def, ll_cut, &self, cache)?;
                }

                // Compute the factor coming from the cut direction
                for ll_cut in cut.iter() {
                    r *= match ll_cut {
                        Cut::NegativeCut(_) => -Complex::one(),
                        _ => Complex::one(),
                    };
                }

                // if overwrite_propagator_power is false it means that the powers used are not the
                // one coming from the topology so it's not possible to use the cached loop line result
                for (p_e, p_pow) in cache
                    .propagators_eval
                    .iter()
                    .zip_eq(cache.propagator_powers.iter())
                {
                    r *= utils::powi(*p_e, *p_pow);
                }
                for cut_prop in cut_propagators[..self.n_loops].iter() {
                    r *= match cut[self.propagator_id_to_ll_id[*cut_prop].0] {
                        Cut::PositiveCut(_) => {
                            cache.complex_cut_energies[*cut_prop] * Into::<T>::into(2.0)
                        }
                        Cut::NegativeCut(_) => {
                            -cache.complex_cut_energies[*cut_prop] * Into::<T>::into(2.0)
                        }
                        Cut::NoCut => panic!("Trying to multiply the energy for uncut propagator!"),
                    };
                }
            }
            //numerator.evaluate(&cut_propagators[..self.n_loops], cache) * r.inv()
            numerator.evaluate_lb(k_def, cache, num_id) * r.inv()
        } else {
            numerator.evaluate_reduced_in_cb(cut, &self.loop_lines, mat, cache, num_id);
            // Cache all propagator evaluations by evaluating the loop line
            for (ll_cut, ll) in cut.iter().zip_eq(self.loop_lines.iter()) {
                ll.evaluate(&k_def, ll_cut, &self, cache)?;
            }

            // Compute the factor coming from the cut direction
            let mut cut_factor = 1;
            for ll_cut in cut.iter() {
                cut_factor *= match ll_cut {
                    Cut::NegativeCut(_) => -1,
                    _ => 1,
                };
            }

            // Compute higher order derivative looping over all numerator monomials
            let mut res = Complex::default();
            for (i, monomial_powers) in numerator
                .reduced_coefficient_index_to_powers
                .iter()
                .enumerate()
            {
                // Take derivative looping over the numerator monomials
                for (num_pow, &mon_pow) in numerator_powers.iter_mut().zip(monomial_powers.iter()) {
                    *num_pow = mon_pow as usize;
                }
                res += cache.reduced_coefficient_cb[i]
                    * self.evaluate_higher_order_cut(
                        cut,
                        &cut_propagators[..self.n_loops],
                        &mut derivative_map[..self.n_loops],
                        &mut numerator_powers[..self.n_loops],
                        mat,
                        0,
                        cache,
                    );
            }
            res.multiply_sign(cut_factor)
            //res
        };

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
        k_def: &[LorentzVector<Complex<T>>],
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

                if self.settings.deformation.scaling.branch_cut_check
                    && cm.re < T::zero()
                    && cm.im < T::zero()
                    && self.settings.deformation.scaling.branch_cut_m < 0.
                    && self.settings.deformation.scaling.source_branch_cut_m < 0.
                {
                    eprintln!(
                        "{}: {} for prop {}, ll sig={:?}, ks={:?}: {}",
                        self.name,
                        "Branch cut detected".red(),
                        p.id,
                        ll.signature,
                        k_def,
                        cm
                    );
                }

                cache.complex_prop_spatial[p.id] = cm;
                cache.complex_cut_energies[p.id] = cm.sqrt();

                // TODO: cache all required powers of the energies
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
        channel: Option<isize>,
    ) -> Result<Complex<T>, &'static str> {
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>;
        let mut shifted_k: ArrayVec<[LorentzVector<T>; MAX_LOOP]> = (0..self.n_loops)
            .map(|_| LorentzVector::default())
            .collect();
        let mut cut_shifts: ArrayVec<[LorentzVector<T>; MAX_LOOP]> = (0..self.n_loops)
            .map(|_| LorentzVector::default())
            .collect();

        let num_cuts = self.ltd_cut_options.iter().map(|c| c.len()).sum();

        let mut sampled_cuts: Option<Vec<usize>> = None;
        if let Some(c) = channel {
            if c < 0 {
                // randomly sample |c| cuts
                // TODO: prevent allocation?
                let seed: [u8; 32] = if let Some(global_seed) = self.global_seed {
                    global_seed
                } else {
                    rand::thread_rng().gen()
                };
                let mut rng: StdRng = SeedableRng::from_seed(seed);
                sampled_cuts = Some(
                    (0..num_cuts)
                        .into_iter()
                        .choose_multiple(&mut rng, -c as usize),
                );
            }
        }

        let mut res: Complex<T> = Complex::default();

        let mut cut_counter = 0;
        for (cuts, mat) in self
            .ltd_cut_options
            .iter()
            .zip_eq(self.cb_to_lmb_mat.iter())
        {
            for cut in cuts {
                if let Some(channel_id) = channel {
                    if channel_id >= 0 && channel_id != cut_counter as isize {
                        cut_counter += 1;
                        continue;
                    }
                }

                if let Some(sc) = &sampled_cuts {
                    if !sc.contains(&cut_counter) {
                        cut_counter += 1;
                        continue;
                    }
                }

                if self.settings.general.debug > 2 {
                    println!("Evaluating channel {}", cut_counter);
                }

                // compute the deformation vector for this propagator
                // shift and rotate the origin to the focus
                for (new_k, sign_map) in shifted_k.iter_mut().zip(mat.chunks_exact(self.n_loops)) {
                    let mut cut_index = 0;
                    for (ll_cut_sig, ll_cut) in cut.iter().zip_eq(&self.loop_lines) {
                        if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = ll_cut_sig {
                            cut_shifts[cut_index] = ll_cut.propagators[*i].q.cast();
                            cut_index += 1;
                        }
                    }

                    *new_k = LorentzVector::default();
                    for ((j, kk), shift) in sign_map[..self.n_loops]
                        .iter()
                        .zip_eq(&k[..self.n_loops])
                        .zip_eq(&cut_shifts[..self.n_loops])
                    {
                        *new_k += (kk - shift).multiply_sign(*j);
                    }
                    new_k.t = T::zero();
                }

                let (kappas, jac) = self.deform(&shifted_k, cache);
                k_def = (0..self.n_loops)
                    .map(|i| {
                        shifted_k[i].map(|x| Complex::new(x, T::zero()))
                            + kappas[i].map(|x| Complex::new(T::zero(), x))
                    })
                    .collect();

                self.compute_complex_cut_energies(&k_def[..self.n_loops], cache)?;

                // compute all multi-channeling factors
                let mut other_cut_counter = 0;
                let mut normalization: Complex<T> = Complex::zero();
                let mut suppresion_factor: Complex<T> = Complex::zero();
                for other_cuts in &self.ltd_cut_options {
                    for other_cut in other_cuts {
                        let mut other_suppresion_factor: Complex<T> = Complex::one();
                        for (other_ll_cut_sig, other_ll_cut) in
                            other_cut.iter().zip_eq(&self.loop_lines)
                        {
                            for (prop_index, prop) in other_ll_cut.propagators.iter().enumerate() {
                                if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = other_ll_cut_sig
                                {
                                    if *i != prop_index {
                                        other_suppresion_factor *=
                                            cache.complex_cut_energies[prop.id];
                                    }
                                } else {
                                    other_suppresion_factor *= cache.complex_cut_energies[prop.id];
                                }
                            }
                        }

                        if other_cut_counter == cut_counter {
                            suppresion_factor = other_suppresion_factor;
                        }

                        normalization += other_suppresion_factor;

                        other_cut_counter += 1;
                    }
                }

                cut_counter += 1;

                let channel_factor: Complex<T> = suppresion_factor / normalization;
                if channel_factor.norm_sqr() < Into::<T>::into(self.e_cm_squared * 1e-16) {
                    continue;
                }

                let r = self.evaluate_all_dual_integrands(&mut k_def, cache)?;

                res += r * jac * channel_factor;
            }
        }

        if let Some(sc) = sampled_cuts {
            res *= Into::<T>::into((num_cuts as f64) / (sc.len() as f64));
        }

        Ok(res)
    }

    pub fn evaluate_all_dual_integrands<T: FloatLike>(
        &self,
        k_def: &mut [LorentzVector<Complex<T>>],
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &'static str> {
        let mut result = Complex::default();
        // Check if partial fraction needs to be used
        let use_partial_fractioning = self.settings.general.partial_fractioning_threshold >= 0.0
            && k_def.len() > 0
            && k_def
                .iter()
                .map(|k| k.real().spatial_squared())
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                .unwrap()
                > Into::<T>::into(self.settings.general.partial_fractioning_threshold)
                    * Into::<T>::into(self.e_cm_squared);
        // Check if need to use multiloop version of partial fractioning
        let use_new_pf = self.settings.general.partial_fractioning_multiloop;

        // Partial fractioning
        if use_partial_fractioning {
            // evaluate the loop lines without a signature, since they will not be part of the fractioning
            let mut inv_result_static = Complex::one();

            // When using the old version it needs to manually compute the porpagators with no loop energy dependence
            // WARNING: only works at 1-loop
            if !use_new_pf {
                for ll in &self.loop_lines {
                    if ll.signature == &[0] {
                        inv_result_static *= ll.evaluate(&k_def, &Cut::NoCut, &self, cache)?;
                    }
                }
            }
            if self.settings.general.use_amplitude {
                // find the loop line that has a loop momentum
                let ll_with_loop = self
                    .loop_lines
                    .iter()
                    .filter(|ll| ll.signature == &[1])
                    .next()
                    .unwrap();
                // Precompute all the ellipsoids that will appear in the partial fractioned expression
                // TODO: limit it to only those that are called
                Topology::evaluate_ellipsoids_matrix_1l(ll_with_loop, cache);
                // Evaluate in the case of an amplitude
                for diag in self.amplitude.diagrams.iter() {
                    if diag.name == "born" || diag.name.contains("UV_") {
                        continue;
                    }

                    //Whole evaluation is taken care by LTDNumerator::evaluate once the correct
                    //powers for the propagators are set and the corresponding numerator is selected
                    for prop_pow in cache.propagator_powers.iter_mut() {
                        *prop_pow = 0;
                    }
                    for (ll_id, &pow) in diag.denominators.iter().zip_eq(diag.pows.iter()) {
                        cache.propagator_powers[self.loop_lines[ll_id.0].propagators[ll_id.1].id] =
                            pow;
                    }

                    // Cache the reduced polynomial for the current diagram and evaluate
                    diag.numerator.evaluate_reduced_in_lb(
                        &k_def,
                        0,
                        cache,
                        0,
                        true,
                        use_partial_fractioning,
                    );
                    result += if diag.ct {
                        -self.partial_fractioning.evaluate(
                            &diag.numerator,
                            std::slice::from_ref(ll_with_loop),
                            cache,
                        )
                    } else {
                        self.partial_fractioning.evaluate(
                            &diag.numerator,
                            std::slice::from_ref(ll_with_loop),
                            cache,
                        )
                    };
                }
            } else {
                // reset the powers, since `evaluate_cut` sets the on-shell propagator power to 0
                for prop_pow in cache.propagator_powers.iter_mut() {
                    *prop_pow = 0;
                }
                for ll in self.loop_lines.iter() {
                    for p in ll.propagators.iter() {
                        cache.propagator_powers[p.id] = p.power;
                    }
                }

                // Evaluate in the case of a topology
                self.numerator.evaluate_reduced_in_lb(
                    &k_def,
                    0,
                    cache,
                    0,
                    true,
                    use_partial_fractioning,
                );

                result = if !use_new_pf {
                    // find the loop line that has a loop momentum
                    let ll_with_loop = self
                        .loop_lines
                        .iter()
                        .filter(|ll| ll.signature == &[1])
                        .next()
                        .unwrap();
                    // Precompute all the ellipsoids that will appear in the partial fractioned expression
                    // TODO: limit it to only those that are called
                    Topology::evaluate_ellipsoids_matrix_1l(ll_with_loop, cache);
                    self.partial_fractioning.evaluate(
                        &self.numerator,
                        std::slice::from_ref(ll_with_loop),
                        cache,
                    )
                } else {
                    self.partial_fractioning_multiloops.evaluate(
                        &self.loop_lines,
                        &self.propagator_id_to_ll_id,
                        cache,
                    )
                };
            }

            result /= inv_result_static;
        } else {
            // Compute all the necessary reduced polynomial using the value of k_vec
            if self.settings.general.use_amplitude {
                for (num_id, diag) in self.amplitude.diagrams.iter().enumerate() {
                    diag.numerator.evaluate_reduced_in_lb(
                        &k_def,
                        0,
                        cache,
                        num_id,
                        true,
                        use_partial_fractioning,
                    );
                }
            } else {
                self.numerator.evaluate_reduced_in_lb(
                    &k_def,
                    0,
                    cache,
                    0,
                    true,
                    use_partial_fractioning,
                );
            }

            // Sum over all the cuts
            let mut cut_counter = 0;
            for (cuts, mat) in self
                .ltd_cut_options
                .iter()
                .zip_eq(self.cb_to_lmb_mat.iter())
            {
                // for each cut coming from the same cut structure
                for cut in cuts.iter() {
                    if self.settings.general.cut_filter.is_empty()
                        || self.settings.general.cut_filter.contains(&cut_counter)
                    {
                        let v = if self.settings.general.use_amplitude {
                            self.evaluate_amplitude_cut(
                                &mut k_def[..self.n_loops],
                                cut,
                                mat,
                                cache,
                            )?
                        } else {
                            let v = self.evaluate_cut(
                                &mut k_def[..self.n_loops],
                                &self.numerator,
                                cut,
                                mat,
                                cache,
                                true,
                                0,
                            )?;
                            // Assuming that there is no need for the residue energy or the cut_id
                            let ct = if self.settings.general.use_ct {
                                self.counterterm(
                                    &k_def[..self.n_loops],
                                    Complex::default(),
                                    0,
                                    cache,
                                )
                            } else {
                                Complex::default()
                            };

                            v * (ct + T::one())
                        };

                        // k_def has the correct energy component at this stage
                        result += v
                    }
                }
                cut_counter += 1;
            }
        }
        Ok(result)
    }

    pub fn evaluate<'a, T: FloatLike>(
        &mut self,
        x: &'a [f64],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        // parameterize
        let mut k = [LorentzVector::default(); MAX_LOOP];
        let mut jac_para = T::one();
        for i in 0..self.n_loops {
            // set the loop index to i + 1 so that we can also shift k
            let (l_space, jac) = Topology::parameterize(
                &x[i * 3..(i + 1) * 3],
                self.e_cm_squared,
                i,
                &self.settings,
            );

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
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>;
        let mut jac_def = Complex::one();

        let r = if self.settings.general.multi_channeling {
            self.evaluate_multi_channel(&k, cache, self.settings.general.multi_channeling_channel)
        } else {
            let (kappas, jac) = self.deform(&k, cache);
            k_def = (0..self.n_loops)
                .map(|i| {
                    k[i].map(|x| Complex::new(x, T::zero()))
                        + kappas[i].map(|x| Complex::new(T::zero(), x))
                })
                .collect();
            jac_def = jac;

            if self
                .compute_complex_cut_energies(&k_def[..self.n_loops], cache)
                .is_err()
            {
                return Complex::default();
            }

            self.evaluate_all_dual_integrands(&mut k_def[..self.n_loops], cache)
        };

        let mut result = match r {
            Ok(v) => v,
            Err(_) => {
                return Complex::default();
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
            result
                + Complex::new(
                    Into::<T>::into(self.amplitude.int_ct[0].re),
                    Into::<T>::into(self.amplitude.int_ct[0].im),
                )
        } else {
            result
        }
    }

    pub fn evaluate_ellipsoids_matrix_1l<T: FloatLike>(
        ll_with_loop: &LoopLine,
        cache: &mut LTDCache<T>,
    ) {
        for p1 in ll_with_loop.propagators.iter() {
            for p2 in ll_with_loop.propagators.iter() {
                cache.complex_ellipsoids[p1.id][p2.id] = cache.complex_cut_energies[p1.id]
                    + cache.complex_cut_energies[p2.id]
                    + Into::<T>::into(-p1.q.t + p2.q.t);
            }
        }
    }
}

impl LTDNumerator {
    pub fn from_sparse(n_loops: usize, coefficients: &[(Vec<usize>, (f64, f64))]) -> LTDNumerator {
        let rank = coefficients.iter().map(|c| c.0.len()).max().unwrap();

        let dense_length = (0..=rank)
            .map(|r| LTDNumerator::binomial(4 * n_loops + r - 1, r))
            .sum();

        let mut dense = vec![Complex::default(); dense_length];

        let sorted_linear: Vec<Vec<usize>> = (0..=rank)
            .map(|rank| (0..4 * n_loops).combinations_with_replacement(rank))
            .flatten()
            .collect();

        for (exp_map, (re, im)) in coefficients {
            let index = if exp_map.len() == 0 {
                0
            } else {
                sorted_linear[..]
                    .binary_search_by(|x| utils::compare_slice(&x[..], &exp_map[..]))
                    .unwrap()
                    + 1
            };

            dense[index] = Complex::new(*re, *im);
        }

        LTDNumerator::new(n_loops, dense.as_slice())
    }

    pub fn new(n_loops: usize, coefficients: &[Complex<f64>]) -> LTDNumerator {
        // Determine max_rank
        let mut max_rank = 0;
        let mut check_size = 1;
        let mut level_size = 1;

        // TODO: check that we don't get an infinite loop
        while check_size < coefficients.len() {
            level_size = (level_size * (n_loops * 4 + max_rank)) / (max_rank + 1);
            max_rank += 1;
            check_size += level_size;
        }
        // Check if the sizes match
        if max_rank != 0 {
            assert_eq!(
                coefficients.len(),
                check_size,
                "Coefficient list must be complete in each rank! For rank {} expected {} elements but found {}", max_rank, check_size, coefficients.len());
        }
        let sorted_linear: Vec<Vec<usize>> = (0..=max_rank)
            .map(|rank| (0..4 * n_loops).combinations_with_replacement(rank))
            .flatten()
            .collect();

        // Create the map that goes from the tensor index to the reduced repesentation
        // that uses only the energy components
        let mut coefficient_index_to_powers = vec![[0; MAX_LOOP]];
        for indices in sorted_linear.iter() {
            let mut key = [0; MAX_LOOP];
            for &n in indices.iter() {
                if n % 4 == 0 {
                    key[n / 4] += 1;
                }
            }
            coefficient_index_to_powers.push(key);
        }

        let coefficient_index_map = sorted_linear
            .iter()
            .map(|l| {
                if l.len() == 1 {
                    (0, l[l.len() - 1])
                } else {
                    (
                        sorted_linear[..]
                            .binary_search_by(|x| utils::compare_slice(&x[..], &l[..l.len() - 1]))
                            .unwrap()
                            + 1,
                        l[l.len() - 1],
                    )
                }
            })
            .collect();

        // Compute reduced block size
        let mut reduced_blocks = Vec::new();
        for r in 0..=max_rank {
            for l in 1..=n_loops {
                let mut last_pos = 0;
                if r != 0 {
                    last_pos = 1;
                    // compute (l+r)!/r!/l!
                    for r_i in 1..=r {
                        // Note the parentesis to ensure integer division!
                        last_pos = (last_pos * (l + r_i)) / r_i;
                    }
                    last_pos -= 1;
                }
                reduced_blocks.push(last_pos);
            }
        }

        // Create power to position maps for the reduced coefficient list
        let mut reduced_size = 1;
        for i in 1..=max_rank {
            reduced_size = (reduced_size * (n_loops + i)) / i
        }
        let mut reduced_coefficient_index_to_powers = vec![[0; MAX_LOOP]; reduced_size];
        for pow_distribution in (0..=n_loops).combinations_with_replacement(max_rank) {
            let mut key = [0; MAX_LOOP];
            for &n in pow_distribution.iter() {
                if n != n_loops {
                    key[n] += 1;
                }
            }
            reduced_coefficient_index_to_powers
                [LTDNumerator::powers_to_position(&key, n_loops, &reduced_blocks)] = key;
        }

        // Copy numerator coefficients into vector
        let mut numerator_coefficients = Vec::new();
        for &coeff in coefficients.iter() {
            numerator_coefficients.push(coeff);
        }

        // Compute a mapping directly from the non-empty numerator entries to the reduced
        // numerator for a varying number of fixed energies
        let mut non_empty_coeff_map_to_reduced_numerator = vec![];
        let mut coeff_map_to_reduced_numerator = vec![];
        for absorb_n_energies in 0..=n_loops {
            let mut numerator_map = vec![];
            let mut full_numerator_map = vec![];
            let mut red_pows = [0; MAX_LOOP];
            for (i, (&c, powers)) in numerator_coefficients
                .iter()
                .zip(coefficient_index_to_powers.iter())
                .enumerate()
            {
                for i in absorb_n_energies..n_loops {
                    red_pows[i] = powers[i];
                }

                let reduced_pos = LTDNumerator::powers_to_position(
                    &red_pows[..n_loops],
                    n_loops,
                    &reduced_blocks,
                );

                if !c.is_zero() {
                    numerator_map.push((i, reduced_pos, c.clone()));
                }
                full_numerator_map.push(reduced_pos);
            }

            non_empty_coeff_map_to_reduced_numerator.push(numerator_map);
            coeff_map_to_reduced_numerator.push(full_numerator_map);
        }

        LTDNumerator {
            coefficients: numerator_coefficients,
            n_loops,
            max_rank,
            reduced_size,
            sorted_linear,
            coefficient_index_map,
            coefficient_index_to_powers,
            reduced_coefficient_index_to_powers,
            reduced_blocks,
            non_empty_coeff_map_to_reduced_numerator,
            coeff_map_to_reduced_numerator,
            coefficients_modified: false,
        }
    }

    /// Generate a numerator that is 1.
    pub fn one(n_loops: usize) -> LTDNumerator {
        LTDNumerator::new(n_loops, &[Complex::one()])
    }

    /// Find the position of a set of powers for a set with maximal rank equal to the sum
    /// of all the powers of the input
    /// reduced blocks contains the precomputed sizes for all such group with smaller or equal
    /// number of variables and/or total rank
    pub fn powers_to_position(basis: &[u8], n_variable: usize, reduced_blocks: &[usize]) -> usize {
        let mut rank: u8 = basis.iter().sum();

        if rank == 0 {
            return 0;
        }
        let mut pos = reduced_blocks[(rank as usize - 1) * n_variable + n_variable - 1];
        for (n, pow) in basis[..n_variable].iter().enumerate() {
            rank -= pow;
            pos += 1;

            match rank {
                0 => return pos,
                1 => continue,
                _ => pos += reduced_blocks[(rank as usize - 1) * n_variable + n_variable - (n + 2)],
            };
        }
        return pos;
    }

    pub fn reduced_powers_to_position(&self, basis: &[u8]) -> usize {
        LTDNumerator::powers_to_position(&basis[..self.n_loops], self.n_loops, &self.reduced_blocks)
    }

    /// Perform the rotation to the vector component of the coefficients
    /// C_i1_i2_..._in -> R_i1_j1... R_in_jn C_j1_j2_..._jn
    /// NOTE: untested for multiloops
    pub fn rotate(&self, rotation_matrix: [[float; 3]; 3]) -> LTDNumerator {
        if self.coefficients.len() == 0 {
            return LTDNumerator::one(self.n_loops);
        }
        let mut coefficients = vec![Complex::default(); self.coefficients.len()];
        coefficients[0] = self.coefficients[0];
        for (indices, coeff) in self
            .sorted_linear
            .iter()
            .zip(self.coefficients.iter().skip(1))
        {
            if coeff.is_zero() {
                continue;
            }

            let rank = indices.len();
            let vec_indices: Vec<&usize> = indices.iter().filter(|&x| x % 4 != 0).collect();
            // When there are only energies there is no need to rotate
            if vec_indices.len() == 0 {
                coefficients[self.sorted_linear[..]
                    .binary_search_by(|x| utils::compare_slice(&x[..], &indices[..]))
                    .unwrap()
                    + 1] += *coeff;
            }

            // Rotate the spatial components
            let mut new_indices = vec![0; rank];
            for map_to in (0..vec_indices.len())
                .map(|_| 0..3)
                .multi_cartesian_product()
            {
                let mut res = *coeff;
                for (&i, j) in vec_indices.iter().zip(map_to.iter()) {
                    res *= rotation_matrix[*j][(*i - 1) % 4]
                }

                let mut index = 0;
                for (new_i, &old_i) in new_indices.iter_mut().zip_eq(indices) {
                    *new_i = match old_i % 4 {
                        0 => old_i,
                        _ => {
                            index += 1;
                            4 * (old_i / 4) + map_to[index - 1] + 1
                        }
                    };
                }
                new_indices.sort_unstable(); // painful but necessary step

                // Find position in coefficient list
                coefficients[self.sorted_linear[..]
                    .binary_search_by(|x| utils::compare_slice(&x[..], &new_indices[..]))
                    .unwrap()
                    + 1] += res;
            }
        }
        LTDNumerator::new(self.n_loops, &coefficients)
    }

    pub fn get_monomial_energy_rec<T: FloatLike>(
        &self,
        i: usize,
        loop_momenta: &[LorentzVector<Complex<T>>],
        absorb_first_energies: usize,
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        if i == 0 {
            return Complex::one();
        }
        if cache.numerator_cache_outdated[i] {
            let (cache_index, vec_index) = self.coefficient_index_map[i - 1];
            cache.numerator_momentum_cache[i] = if vec_index % 4 == 0
                && vec_index / 4 >= absorb_first_energies
            {
                // energy component
                self.get_monomial_energy(cache_index, loop_momenta, absorb_first_energies, cache)
            } else {
                // vector or fixed energy component
                self.get_monomial_energy(cache_index, loop_momenta, absorb_first_energies, cache)
                    * loop_momenta[vec_index / 4][vec_index % 4]
            };
        }

        cache.numerator_cache_outdated[i] = false;
        cache.numerator_momentum_cache[i]
    }

    #[inline(always)]
    pub fn get_monomial_energy<T: FloatLike>(
        &self,
        i: usize,
        loop_momenta: &[LorentzVector<Complex<T>>],
        absorb_first_energies: usize,
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        if i == 0 {
            return Complex::one();
        } else if !cache.numerator_cache_outdated[i] {
            cache.numerator_momentum_cache[i]
        } else {
            self.get_monomial_energy_rec(i, loop_momenta, absorb_first_energies, cache)
        }
    }

    /// Change from the loop momentum basis to the cut basis for a polynomial involving only
    /// the energy components of the loop momenta.
    pub fn change_monomial_basis<T: FloatLike>(
        &self,
        powers: &[u8; MAX_LOOP],
        mat: &Vec<i8>,
        shifts: &[Complex<T>; MAX_LOOP],
        basis: &mut [u8; MAX_LOOP],
        coeff: Complex<T>,
        index: usize,
        cache: &mut LTDCache<T>,
    ) {
        // If the index is zero reset the power vector in the new basis
        if index == 0 {
            for i in 0..self.n_loops {
                basis[i] = 0;
            }
        }
        // Change of basis complete
        if index == self.n_loops {
            //println!("\t basis: {:?} [{}]", basis, coeff);
            cache.reduced_coefficient_cb[self.reduced_powers_to_position(basis)] += coeff;
            return;
        }
        // Skip if no power
        if powers[index] == 0 {
            self.change_monomial_basis(powers, mat, shifts, basis, coeff, index + 1, cache);
            return;
        }

        // FIXME: maximum power hardcoded to 10
        // this entry is hard to cache, since we are in a recursive function and we should keep
        // the state
        let mut pow_distribution = [0; 10];
        let power = powers[index] as usize;

        loop {
            let mut monomial = [0; MAX_LOOP];
            for &n in pow_distribution[..power].iter() {
                monomial[n] += 1;
            }
            //Update coefficient
            //println!("\tmonomial: {:?}, shift: {}", monomial, shifts[index]);
            let mut new_coeff = coeff
                * utils::powi(shifts[index], monomial[self.n_loops])
                * T::from_usize(LTDNumerator::multinomial(&monomial[..=self.n_loops])).unwrap();
            for (&mij, &pow_j) in mat
                .chunks_exact(self.n_loops)
                .nth(index)
                .unwrap()
                .iter()
                .zip_eq(monomial[..self.n_loops].iter())
            {
                new_coeff *= utils::powi(Complex::new(T::from_i8(mij).unwrap(), T::zero()), pow_j);
            }
            //Update basis
            for (b, m) in basis.iter_mut().zip(monomial[0..self.n_loops].iter()) {
                *b += *m as u8;
            }
            // Move to next step
            self.change_monomial_basis(powers, mat, shifts, basis, new_coeff, index + 1, cache);
            //Reset basis
            for (b, m) in basis.iter_mut().zip(monomial[0..self.n_loops].iter()) {
                *b -= *m as u8;
            }

            if !utils::next_combination_with_replacement(
                &mut pow_distribution[..power],
                self.n_loops,
            ) {
                break;
            }
        }
    }

    pub fn evaluate_reduced_in_lb<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<Complex<T>>],
        absorb_n_energies: usize,
        cache: &mut LTDCache<T>,
        num_id: usize,
        clear_momentum_cache: bool,
        use_partial_fractioning: bool,
    ) {
        // Update tensor loop dependent part
        // Resize cache vector in case its size is not sufficient
        if cache.numerator_momentum_cache.len() < self.coefficient_index_map.len() + 1 {
            cache
                .numerator_cache_outdated
                .resize(self.coefficient_index_map.len() + 1, true);
            cache.numerator_momentum_cache.resize(
                self.coefficient_index_map.len() + 1,
                Complex::new(T::zero(), T::zero()),
            );
        }

        if clear_momentum_cache {
            for i in &mut cache.numerator_cache_outdated {
                *i = true;
            }

            cache.numerator_momentum_cache[0] = Complex::one();
            cache.numerator_cache_outdated[0] = false;
        }

        // Initialize the reduced_coefficients_lb with the factors coming from evaluating
        // the vectorial part of the loop momenta
        if cache.reduced_coefficient_lb[num_id].len() < self.reduced_size {
            cache.reduced_coefficient_lb[num_id].resize(self.reduced_size, Complex::default());
        }

        for c in &mut cache.reduced_coefficient_lb[num_id][..self.reduced_size] {
            *c = Complex::default();
        }

        if !self.coefficients_modified {
            for (i, reduced_index, c) in
                self.non_empty_coeff_map_to_reduced_numerator[absorb_n_energies].iter()
            {
                let mm = self.get_monomial_energy(*i, loop_momenta, absorb_n_energies, cache);
                cache.reduced_coefficient_lb[num_id][*reduced_index] +=
                    mm * Complex::new(Into::<T>::into(c.re), Into::<T>::into(c.im));
            }
        } else {
            for (i, (c, reduced_index)) in self
                .coefficients
                .iter()
                .zip(self.coeff_map_to_reduced_numerator[absorb_n_energies].iter())
                .enumerate()
            {
                if c.is_zero() {
                    continue;
                }

                let mm = self.get_monomial_energy(i, loop_momenta, absorb_n_energies, cache);
                cache.reduced_coefficient_lb[num_id][*reduced_index] +=
                    mm * Complex::new(Into::<T>::into(c.re), Into::<T>::into(c.im));
            }
        }

        // Store coefficients in a multi variable polynomial
        if use_partial_fractioning {
            cache.reduced_coefficient_lb_mpoly.clear();
            for (c, pows) in cache.reduced_coefficient_lb[num_id]
                .iter()
                .zip(self.reduced_coefficient_index_to_powers[..self.reduced_size].iter())
            {
                cache
                    .reduced_coefficient_lb_mpoly
                    .add(&pows[..self.n_loops], *c);
            }
        }
    }

    pub fn evaluate_reduced_in_cb<T: FloatLike>(
        &self,
        cut: &Vec<Cut>,
        loop_lines: &Vec<LoopLine>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
        num_id: usize,
    ) {
        let mut shifts = [Complex::default(); MAX_LOOP];
        // compute the shifts of the propagator in the cut basis
        let mut index = 0;
        for (ll_cut, ll) in cut.iter().zip_eq(loop_lines.iter()) {
            if let Cut::PositiveCut(j) | Cut::NegativeCut(j) = ll_cut {
                // Map from lb to cb
                for (m_row, shift) in mat
                    .chunks_exact(self.n_loops)
                    .zip_eq(shifts[..self.n_loops].iter_mut())
                {
                    *shift -= Into::<T>::into(ll.propagators[*j].q.t * m_row[index] as f64);
                }
                index += 1;
            }
        }
        // Initialize the reduced_coefficeints_cb with the factors coming from evaluating
        // the vectorial part of the loop momenta
        if cache.reduced_coefficient_cb.len() < self.reduced_size {
            cache
                .reduced_coefficient_cb
                .resize(self.reduced_size, Complex::default());
        }
        for c in &mut cache.reduced_coefficient_cb[..self.reduced_size] {
            *c = Complex::default();
        }

        for (index, powers_lb) in self.reduced_coefficient_index_to_powers.iter().enumerate() {
            if cache.reduced_coefficient_lb[num_id][index].is_zero() {
                continue;
            }

            let mut powers_cb = [0; MAX_LOOP];
            //println!(
            //    "MONOMIAL: {:?} * {:?}",
            //    powers_lb, cache.reduced_coefficient_lb[num_id][index]
            //);
            self.change_monomial_basis(
                powers_lb,
                mat,
                &shifts,
                &mut powers_cb,
                cache.reduced_coefficient_lb[num_id][index],
                0,
                cache,
            );
        }
    }

    pub fn evaluate_lb<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<Complex<T>>],
        cache: &mut LTDCache<T>,
        num_id: usize,
    ) -> Complex<T> {
        let mut res = Complex::default();
        for (numerator_powers, coeff) in self
            .reduced_coefficient_index_to_powers
            .iter()
            .zip_eq(cache.reduced_coefficient_lb[num_id][..self.reduced_size].iter())
        {
            if coeff.is_zero() {
                continue;
            }

            let mut prod = Complex::one();
            for (l_i, pow) in numerator_powers[..self.n_loops].iter().enumerate() {
                prod *= loop_momenta[l_i].t.powi(*pow as i32);
            }
            res += coeff * prod;
        }
        res
    }
    pub fn evaluate<T: FloatLike>(
        &self,
        cut_propagators: &[usize],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        let mut res = Complex::default();
        for (numerator_powers, coeff) in self
            .reduced_coefficient_index_to_powers
            .iter()
            .zip_eq(cache.reduced_coefficient_cb[..self.reduced_size].iter())
        {
            if coeff.is_zero() {
                continue;
            }

            let mut prod = Complex::one();
            for (cut_prop, num) in cut_propagators
                .iter()
                .zip_eq(numerator_powers[..self.n_loops].iter())
            {
                prod *= cache.complex_cut_energies[*cut_prop].powi(*num as i32);
            }
            res += coeff * prod;
        }
        res
    }

    // maximal value of n: 62
    pub fn binomial(n: usize, k: usize) -> usize {
        if k > n {
            return 0;
        }
        if k > n - k {
            return LTDNumerator::binomial(n, n - k);
        }
        let mut res = 1;
        for i in 1..=k {
            res *= n - i + 1;
            res /= i;
        }
        res
    }

    pub fn multinomial(n: &[usize]) -> usize {
        let mut res = 1;
        let mut sum = 0;
        for &n_i in n.iter() {
            sum += n_i;
            res *= LTDNumerator::binomial(sum, n_i);
        }
        res
    }
}

//impl<T: FloatLike> ReducedLTDNumerator<T> {
//    pub fn evaluate_reduced_in_lb(
//        &mut self,
//        numerator: &LTDNumerator,
//        loop_momenta: &[LorentzVector<Complex<T>>],
//        absorb_n_energies: usize,
//        cache: &mut LTDCache<T>,
//    ) {
//        // Get information about the numerator
//        self.max_rank = numerator.max_rank;
//        self.n_loops = numerator.n_loops;
//        self.size = numerator.reduced_size;
//        self.reduced_coefficient_index_to_powers = numerator.reduced_coefficient_index_to_powers;
//
//        // Update tensor loop dependent part
//        numerator.update_numerator_momentum_some_energies(loop_momenta, absorb_n_energies, cache);
//        // Initialize the reduced_coefficeints_lb with the factors coming from evaluating
//        // the vectorial part of the loop momenta
//        if self.coefficients.len() < self.size {
//            self.coefficients.resize(self.size, Complex::default());
//        }
//        for c in &mut self.coefficients[..self.size] {
//            *c = Complex::default();
//        }
//
//        for (i, (&c, powers)) in numerator
//            .coefficients
//            .iter()
//            .zip(numerator.coefficient_index_to_powers.iter())
//            .enumerate()
//        {
//            self.coefficients[numerator.powers_to_position[powers]] += cache
//                .numerator_momentum_cache[i]
//                * Complex::new(Into::<T>::into(c.re), Into::<T>::into(c.im));
//        }
//    }
//    pub fn evaluate_reduced_in_cb<T: FloatLike>(
//        &self,
//        cut: &Vec<Cut>,
//        loop_lines: &Vec<LoopLine>,
//        mat: &Vec<i8>,
//        cache: &mut LTDCache<T>,
//    ) {
//        let mut shifts = [Complex::default(); MAX_LOOP];
//        // compute the shifts of the propagator in the cut basis
//        let mut index = 0;
//        for (ll_cut, ll) in cut.iter().zip_eq(loop_lines.iter()) {
//            if let Cut::PositiveCut(j) | Cut::NegativeCut(j) = ll_cut {
//                for (&mij, shift) in mat
//                    .chunks_exact(self.n_loops)
//                    .nth(index)
//                    .unwrap()
//                    .iter()
//                    .zip_eq(shifts[..self.n_loops].iter_mut())
//                {
//                    *shift -= Into::<T>::into(ll.propagators[*j].q.t * mij as f64);
//                }
//                index += 1;
//            }
//        }
//    }
//}
//
