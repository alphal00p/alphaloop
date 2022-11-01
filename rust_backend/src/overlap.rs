use crate::squared_topologies::{DeformationSubgraph, Propagator, Threshold};
use crate::MAX_LOOP;
use crate::{utils, FloatLike};
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use num::{Complex, ToPrimitive};
use num_traits::{Float, Signed, Zero};
use scs;
use serde::Deserialize;
use std::fmt::Display;
use std::mem;
use std::ptr;
use std::slice;
use std::time::Instant;
use utils::Signum;

#[derive(Debug, Clone, Deserialize)]
pub struct DeformationOverlap {
    pub deformation_sources: Vec<LorentzVector<f64>>,
    pub excluded_surface_ids: Vec<Vec<((usize, usize), i8, i8)>>,
    pub overlap: Option<Vec<usize>>,
    #[serde(default)]
    pub radius: f64,
    #[serde(skip_deserializing)]
    pub excluded_surface_indices: Vec<usize>,
}

impl Threshold {
    /// Evaluate an E-surface, where eval < 0 is the inside.
    pub fn evaluate_surface<T: FloatLike>(&self, loop_momenta: &[LorentzVector<T>]) -> T {
        let mut eval = Into::<T>::into(self.shift.t);
        for f in &self.focal_points {
            // get the loop momentum
            let mut mom = f.shift.cast();
            for (m, sig) in loop_momenta.iter().zip(&f.sig) {
                mom += m.multiply_sign(*sig);
            }
            eval +=
                (mom.spatial_squared() + Into::<T>::into(f.mass) * Into::<T>::into(f.mass)).sqrt();
        }
        eval
    }

    /// Evaluate a surface with complex momenta.
    pub fn evaluate_complex<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<num::Complex<T>>],
    ) -> num::Complex<T> {
        let mut res = Complex::new(Into::<T>::into(self.shift.t), T::zero());

        for f in &self.focal_points {
            res += ((utils::evaluate_signature(&f.sig, loop_momenta)
                + f.shift.map(|x| Complex::new(Into::<T>::into(x), T::zero())))
            .spatial_squared()
                + Into::<T>::into(f.mass) * Into::<T>::into(f.mass))
            .sqrt();
        }
        res
    }

    /// Generate a point inside an E-surface. A point that is guaranteed to be on the inside is
    /// `c_i = -p s_i m_i/(sum_j m_j)`, in the focus basis where `s` is the surface signature in the focus basis.
    pub fn get_interior_point<T: FloatLike>(&self) -> Vec<LorentzVector<T>> {
        let mass_sum: T = self
            .focal_points
            .iter()
            .map(|f| Into::<T>::into(f.mass))
            .sum();

        let mut cb_on_foci = [LorentzVector::default(); MAX_LOOP];

        // compute the shift of the dependent focus in the focus basis
        let mut shift = self.focal_points.last().unwrap().shift;
        for (m, f) in self
            .sig_in_fb
            .iter()
            .zip(&self.focal_points[..self.focal_points.len() - 1])
        {
            shift += f.shift.multiply_sign(-*m);
        }

        // construct the interior point in the focus basis
        for ((m, f), s) in cb_on_foci
            .iter_mut()
            .zip(&self.focal_points[..self.foci.len() - 1])
            .zip(&self.sig_in_fb)
        {
            if !mass_sum.is_zero() && !f.mass.is_zero() {
                *m = -shift.cast().multiply_sign(*s) * Into::<T>::into(f.mass) / mass_sum;
            }
        }

        // do the basis transformation from the focus basis to the lmb
        let mut ll_on_foci = vec![LorentzVector::default(); self.focal_points[0].sig.len()];
        for (r, lm) in self
            .fb_to_cmb
            .chunks(self.focal_points[0].sig.len())
            .zip(&mut ll_on_foci)
        {
            for ((cm, s), f) in cb_on_foci.iter().zip(r).zip(&self.focal_points) {
                if *s != 0. {
                    *lm += (cm - f.shift.cast()) * Into::<T>::into(*s);
                }
            }
        }

        assert!(
            self.evaluate_surface(&ll_on_foci) < Into::<T>::into(1e-8 * self.shift.t.abs()),
            "Point {:?} not in interior of {}",
            ll_on_foci,
            self
        );
        ll_on_foci
    }
}

impl Display for Threshold {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for focus in &self.focal_points {
            write!(
                f,
                "+|({}{})^2{}|",
                focus
                    .sig
                    .iter()
                    .enumerate()
                    .filter_map(|(i, s)| if *s != 0 {
                        Some(format!("{}k{}", if *s == 1 { "+" } else { "-" }, i + 1))
                    } else {
                        None
                    })
                    .join("+"),
                if focus.shift.x != 0. || focus.shift.y != 0. || focus.shift.z != 0. {
                    format!(
                        "+({:.2e},{:.2e},{:.2e})",
                        focus.shift.x, focus.shift.y, focus.shift.z
                    )
                } else {
                    String::new()
                },
                if focus.mass > 0. {
                    format!("+{}^2", focus.mass)
                } else {
                    String::new()
                }
            )?;
        }
        write!(
            f,
            "{:+.2e}; exists={}, pinched={}",
            self.shift.t, self.exists, self.pinched
        )
    }
}

impl SOCPProblem {
    /// Construct the SOCP problem.
    fn construct(
        &mut self,
        ellipsoid_ids: &[usize],
        ellipsoids: &[&Threshold],
        focal_points: &[Propagator],
        minimize: bool,
    ) {
        // TODO: change the n_loops to the number of specific variables
        let multiplier = if minimize { 6 * self.n_loops } else { 1 };
        self.radius_computation = minimize;

        self.cone.l = (ellipsoid_ids.len() * multiplier) as scs::scs_int;

        let mut b_index = 0;

        for _ in 0..multiplier {
            for &e_id in ellipsoid_ids {
                self.b[b_index] = -ellipsoids[e_id].shift.t;
                b_index += 1;
            }
        }

        let mut lm_tag = [false; MAX_LOOP];

        // construct the focus part of the shift
        // TODO: do not generate the focus more than once if k[m/3] does not appear in the focus!!!
        let mut focus_count = 0;
        for m in 0..multiplier {
            for prop in focal_points {
                // check if the focus occurs in an ellipsoid
                'el_loop: for &e_id in ellipsoid_ids {
                    for (foc_index, f) in ellipsoids[e_id].focal_points.iter().enumerate() {
                        if f.id == prop.id {
                            for (tag, &sig) in lm_tag.iter_mut().zip(&f.sig) {
                                *tag |= sig != 0;
                            }

                            self.b[b_index] = 0.;
                            b_index += 1;
                            if f.mass != 0. {
                                self.b[b_index] = f.mass;
                                b_index += 1;
                                self.q[focus_count] = 5;
                            } else {
                                self.q[focus_count] = 4;
                            };
                            self.b[b_index] = f.shift.x;
                            b_index += 1;
                            self.b[b_index] = f.shift.y;
                            b_index += 1;
                            self.b[b_index] = f.shift.z;
                            b_index += 1;

                            self.focus_list[focus_count] = (prop.id, e_id, foc_index, m);
                            focus_count += 1;
                            break 'el_loop;
                        }
                    }
                }
            }
        }

        let mut var_count = 0;
        for (vm, &x) in self.var_map.iter_mut().zip(lm_tag.iter()) {
            if x {
                *vm = var_count;
                var_count += 1;
            } else {
                *vm = 1000;
            }
        }

        let width = if minimize {
            3 * var_count + focus_count + 1
        } else {
            3 * var_count + focus_count
        };

        self.cone.qsize = focus_count as scs::scs_int;
        self.a.m = b_index as scs::scs_int;
        self.a.n = width as scs::scs_int;
        self.data.m = self.a.m;
        self.data.n = self.a.n;

        for aa in &mut self.a_dense[..width * b_index] {
            *aa = 0.;
        }

        for cc in &mut self.c[..width] {
            *cc = 0.;
        }

        let focus_start = if minimize {
            // minimize the radius. ECOS very rarely gives a proof of infeasibility if we maximize instead
            self.c[3 * var_count] = 1.;
            3 * var_count + 1
        } else {
            3 * var_count
        };

        // write the ellipsoid constraints: f1 + f2 < c, etc
        let mut row_counter = 0; // the start of the focus information
        for m in 0..multiplier {
            for &e_id in ellipsoid_ids {
                for f in ellipsoids[e_id].focal_points.iter() {
                    let foc_index = self
                        .focus_list
                        .iter()
                        .position(|x| x.0 == f.id && x.3 == m)
                        .unwrap();
                    self.a_dense[row_counter * width + focus_start + foc_index] = 1.;
                }
                row_counter += 1;
            }
        }

        // now for every focus f, add the norm constraints |k_x+p_x,m| < f
        for (foc_index, (_, foc_e_surf, foc_e_surf_foc_index, radius_index)) in
            self.focus_list[..focus_count].iter().enumerate()
        {
            let focus = &ellipsoids[*foc_e_surf].focal_points[*foc_e_surf_foc_index];

            // set the linear equation
            self.a_dense[row_counter * width + focus_start + foc_index] = -1.0;
            row_counter += 1;
            if focus.mass != 0. {
                // the mass is an empty line
                row_counter += 1;
            }

            let (rad_lm, rad_dir, rad_sign) =
                (radius_index / 6, radius_index % 3, radius_index % 2);

            for dir_index in 0..3 {
                for (var_index, &v) in focus.sig.iter().enumerate() {
                    if v != 0 {
                        self.a_dense
                            [row_counter * width + 3 * self.var_map[var_index] + dir_index] =
                            Into::<f64>::into(-v);

                        if minimize && rad_lm == var_index && rad_dir == dir_index {
                            self.a_dense[row_counter * width + 3 * var_count] =
                                if rad_sign == 0 { -1. } else { 1.0 };
                        }
                    }
                }
                row_counter += 1;
            }
        }

        // now write a in compressed column format
        let mut non_zero_index = 0;
        self.p[0] = 0;
        for col in 0..width {
            for row in 0..b_index {
                if self.a_dense[row * width + col] != 0. {
                    self.x[non_zero_index] = self.a_dense[row * width + col];
                    self.i[non_zero_index] = row as scs::scs_int;
                    non_zero_index += 1;
                }
            }
            self.p[col + 1] = non_zero_index as scs::scs_int;
        }
    }

    /// Generate a next plausible solution to check for overlaps.
    fn next_set(
        &self,
        num_ellipsoids: usize,
        acc: &mut [isize],
        different: &mut [bool],
        initialize: bool,
    ) -> bool {
        let mut i = if initialize {
            acc[0] = -1;
            0
        } else {
            acc.len() - 1
        };

        if acc.len() == 1 {
            // special case for n = 1: check if the ellipsoid does not appear in any overlap structure
            'next_num_loop1: for new_num in ((acc[0] + 1) as usize)..num_ellipsoids {
                for appears in &self.pair_appears_in_overlap_structurex
                    [new_num as usize * num_ellipsoids..(new_num as usize + 1) * num_ellipsoids]
                {
                    if *appears {
                        continue 'next_num_loop1;
                    }
                }
                acc[0] = new_num as isize;
                return true;
            }
            return false;
        }

        loop {
            if acc[i] != num_ellipsoids as isize {
                'next_num_loop: for new_num in
                    ((acc[i] + 1) as usize)..=(num_ellipsoids - acc.len() + i)
                {
                    different[i] = false;
                    for x in acc[..i].iter() {
                        if !self.pair_overlap_matrix[*x as usize * num_ellipsoids + new_num] {
                            continue 'next_num_loop;
                        }

                        if !self.pair_appears_in_overlap_structurex
                            [*x as usize * num_ellipsoids + new_num]
                        {
                            different[i] = true;
                        }
                    }

                    if i > 0 {
                        different[i] |= different[i - 1];
                    }

                    acc[i] = new_num as isize;

                    // check if there is still a possible completion
                    let spots_left = acc.len() - i - 1;
                    if !different[i] {
                        let mut possible = false;
                        'outer_loop: for j in (acc[i] as usize + 1)..num_ellipsoids {
                            for k in 0..j {
                                if self.pair_overlap_matrix[k * num_ellipsoids + j]
                                    && !self.pair_appears_in_overlap_structurex
                                        [k * num_ellipsoids + j]
                                    && ((k as isize > acc[i] && spots_left >= 2)
                                        || (spots_left >= 1
                                            && acc[..i + 1].contains(&(k as isize))))
                                {
                                    possible = true;
                                    break 'outer_loop;
                                }
                            }
                        }

                        if !possible {
                            continue 'next_num_loop;
                        }
                    }

                    if i + 1 == acc.len() {
                        if different[i] {
                            return true;
                        }
                        // ideally this is never reached
                        continue;
                    }

                    // continue the same logic for the next index i + 1
                    acc[i + 1] = new_num as isize;
                    i += 2; // add 2 because we subtract by one later
                    break;
                }
            }

            if i == 0 {
                break;
            }
            i -= 1;
        }
        false
    }

    /// Assuming the test for overlap of `ellipsoid_list` succeeded, build a `FixedDeformationOverlap`.
    fn build_deformation_overlap_info(
        &self,
        ellipsoid_list: &[usize],
        all_ellipsoids: &[&Threshold],
    ) -> DeformationOverlap {
        // TODO: recycle from pre-allocated overlap structures
        let mut var_count = 0;
        let mut sources = vec![LorentzVector::default(); self.n_loops];

        // map the source back from the loop momenta present in the ellipsoids to the full space
        for (source, &i) in sources.iter_mut().zip_eq(&self.var_map[..self.n_loops]) {
            if i < self.n_loops {
                *source = LorentzVector::from_args(
                    0.,
                    self.sol_x[i * 3],
                    self.sol_x[i * 3 + 1],
                    self.sol_x[i * 3 + 2],
                );
                var_count += 1;
            }
        }

        let mut excluded_surface_indices = Vec::with_capacity(all_ellipsoids.len());
        for (i, e) in all_ellipsoids.iter().enumerate() {
            if !ellipsoid_list.contains(&i) {
                excluded_surface_indices.push(e.id);
            }
        }

        let mut radius = 0.;

        if self.radius_computation {
            radius = -self.sol_x[var_count * 3];
        } else {
            // give a meausure for distance
            for (foc_count, foc_val) in self.c[3 * var_count..self.a.n as usize]
                .iter()
                .zip(&self.sol_x[3 * var_count..self.a.n as usize])
            {
                radius += foc_count * foc_val;
            }
        }

        let overlap = ellipsoid_list
            .iter()
            .map(|e_index| all_ellipsoids[*e_index].id)
            .collect();

        DeformationOverlap {
            deformation_sources: sources,
            excluded_surface_ids: vec![],
            excluded_surface_indices,
            overlap: Some(overlap),
            radius,
        }
    }
}

impl DeformationSubgraph {
    /// Find a point on a surface along a line and return the scale:
    /// `E(dir * t + k) = 0`, returning `t` closest to `t_start`.
    fn get_scale_to_reach_surface<T: FloatLike>(
        &self,
        s: &Threshold,
        k: &[LorentzVector<T>],
        dir: &[LorentzVector<T>],
        t_start: T,
    ) -> Option<T> {
        let mut t = t_start;
        let tolerance = Into::<T>::into(1e-15 * self.e_cm_squared.sqrt());
        for it in 0..30 {
            let mut f = Into::<T>::into(s.shift.t);
            let mut df = T::zero();

            for focus in &s.focal_points {
                let shift = focus.shift.cast();
                let k = utils::evaluate_signature(&focus.sig, k);
                let dir = utils::evaluate_signature(&focus.sig, dir);
                let energy = ((k + dir * t + shift).spatial_squared()
                    + Into::<T>::into(focus.mass) * Into::<T>::into(focus.mass))
                .sqrt();

                if !energy.is_zero() {
                    f += energy;
                    df += (t * dir.spatial_squared() + (k + shift).spatial_dot(&dir)) / energy;
                }
            }

            if self.settings.general.debug > 4 {
                println!("  | surf scaling: f={}, df={}, t={}", f, df, t);
            }

            if Float::abs(f) < tolerance {
                if self.settings.general.debug > 2 {
                    println!("  | t = {}", t);
                }

                return Some(t);
            }

            if it == 29 {
                if self.settings.general.debug > 2 {
                    println!(
                        "  | no convergence after {} iterations: f={}, df={}, t={}",
                        it, f, df, t
                    );
                }
                return None;
            }

            t = t - f / df;
        }
        None
    }

    /// Get the normal at the point on an E-surface. We always assume that the surface is positive.
    fn get_normal<T: FloatLike>(
        &self,
        surface: &Threshold,
        k: &[LorentzVector<T>],
        normalize: bool,
    ) -> [LorentzVector<T>; MAX_LOOP] {
        let mut df = [LorentzVector::default(); MAX_LOOP];

        for f in &surface.focal_points {
            let shift = f.shift.cast();
            let mut mom_part = utils::evaluate_signature(&f.sig, k) + shift;
            mom_part.t = T::zero();
            let energy = (mom_part.spatial_squared()
                + Into::<T>::into(f.mass) * Into::<T>::into(f.mass))
            .sqrt();

            if energy == T::zero() || Float::is_nan(energy) {
                continue;
            }

            for (d, s) in df.iter_mut().zip(&f.sig) {
                if *s != 0 {
                    *d += mom_part.multiply_sign(*s) / energy;
                }
            }
        }

        // normalize
        if normalize {
            let mut norm = T::zero();
            for d in &df[..self.n_loops] {
                norm += d.spatial_squared();
            }
            norm = norm.sqrt().inv();

            for d in &mut df[..self.n_loops] {
                *d *= norm;
            }
        }

        df
    }

    /// Based on "On the distance between two ellipsoids" by Anhua Lin and Shih-Ping Han.
    /// It may not always converge, and thus can only be used as a heuristic.
    fn overlap_heuristic(&self, surf1: &Threshold, surf2: &Threshold) -> bool {
        let mut x1 = [LorentzVector::default(); MAX_LOOP];
        let mut x2 = [LorentzVector::default(); MAX_LOOP];
        x1[..self.n_loops].copy_from_slice(&surf1.interior_point);
        x2[..self.n_loops].copy_from_slice(&surf2.interior_point);

        if surf2.evaluate_surface(&x1[..self.n_loops]) < 0.
            || surf1.evaluate_surface(&x2[..self.n_loops]) < 0.
        {
            return true;
        }

        let mut dir = [LorentzVector::default(); MAX_LOOP];

        let mut previous_distance = -100.;
        let mut step_size;
        for it in 0..50 {
            if it < 10 {
                step_size = 0.1 + 0.39 / 10. * (10 - it) as f64;
            } else {
                step_size = 0.05;
            }

            for i in 0..self.n_loops {
                dir[i] = x2[i] - x1[i];
            }

            debug_assert!(dir[..self.n_loops].iter().any(|d| d.spatial_squared() > 0.));

            // find points on the surfaces
            let t1 = self
                .get_scale_to_reach_surface(surf1, &x1[..self.n_loops], &dir[..self.n_loops], 1.0)
                .unwrap();
            let t2 = self
                .get_scale_to_reach_surface(surf2, &x1[..self.n_loops], &dir[..self.n_loops], 0.0)
                .unwrap();

            if t1 < 0. || t2 < t1 {
                // there is overlap!
                return true;
            }

            // construct the points on the surfaces
            for i in 0..self.n_loops {
                x2[i] = x1[i] + dir[i] * t2;
                x1[i] = x1[i] + dir[i] * t1;
            }
            debug_assert!(
                surf1.evaluate_surface(&x1[..self.n_loops]).abs()
                    < 1e-14 * self.e_cm_squared.to_f64().unwrap()
            );
            debug_assert!(
                surf2.evaluate_surface(&x2[..self.n_loops]).abs()
                    < 1e-14 * self.e_cm_squared.to_f64().unwrap()
            );

            // get the distance
            let mut dist = 0.;
            for i in 0..self.n_loops {
                dist += (x1[i] - x2[i]).spatial_squared();
            }
            dist = dist.sqrt();

            if self.settings.general.debug > 3 {
                println!(
                    "Current distance: {}, step size: {}",
                    dist.sqrt(),
                    step_size
                );
            }

            if (dist - previous_distance).abs() / previous_distance.abs() < 1e-4 {
                if self.settings.general.debug > 2 {
                    println!("No overlap: distance={}", dist);
                }
                return false;
            }
            previous_distance = dist;

            // get the normal
            let x1_norm = self.get_normal(surf1, &x1[..self.n_loops], true);
            let x2_norm = self.get_normal(surf2, &x2[..self.n_loops], true);

            // now solve for the other point
            let t1_distance = self
                .get_scale_to_reach_surface(
                    surf1,
                    &x1[..self.n_loops],
                    &x1_norm[..self.n_loops],
                    -self.e_cm_squared.to_f64().unwrap(), // start with an underestimate
                )
                .unwrap();

            let t2_distance = self
                .get_scale_to_reach_surface(
                    surf2,
                    &x2[..self.n_loops],
                    &x2_norm[..self.n_loops],
                    -self.e_cm_squared,
                )
                .unwrap();

            if t1_distance >= -1e-13 || t2_distance >= -1e-13 {
                // TODO: understand why this can happen. Is the ellipsoid very pinched?
                if self.settings.general.debug > 2 {
                    println!(
                        "Negative distance encountered: {} and {}",
                        t1_distance, t2_distance
                    );
                }
                return false;
            }

            for i in 0..self.n_loops {
                x1[i] = x1[i] + x1_norm[i] * t1_distance * step_size;
                x2[i] = x2[i] + x2_norm[i] * t2_distance * step_size;
            }

            debug_assert!(surf1.evaluate_surface(&x1[..self.n_loops]) < 0.);
            debug_assert!(surf2.evaluate_surface(&x2[..self.n_loops]) < 0.);

            if it == 49 {
                if self.settings.general.debug > 2 {
                    println!(
                        "  | no convergence after {} iterations for surface {} and {}",
                        it, surf1.id, surf2.id
                    );
                }
            }
        }
        false
    }

    // TODO: add support for subspaces again!
    pub fn find_overlap_structure(&mut self, _include_subspaces: bool) -> Vec<DeformationOverlap> {
        let mut overlap_structure = vec![];
        // TODO: just skip the non-existing surfaces in this function
        let ellipsoid_list: Vec<_> = self
            .thresholds
            .iter()
            .filter(|t| t.exists && !t.pinched)
            .collect();

        if ellipsoid_list.len() == 0 {
            return vec![];
        }

        let t = Instant::now();

        let mut problem_count = 0;
        let mut pair_overlap_count = 0;

        if self.settings.deformation.fixed.use_heuristic_centers {
            if ellipsoid_list.len() == 1 {
                // put ourselves on a focus
                let r: f64 = -ellipsoid_list[0]
                    .evaluate_surface(&ellipsoid_list[0].interior_point[..self.n_loops]);
                assert!(r > 0.);

                return vec![DeformationOverlap {
                    deformation_sources: ellipsoid_list[0].interior_point.clone(),
                    excluded_surface_ids: vec![],
                    excluded_surface_indices: vec![],
                    overlap: Some(vec![ellipsoid_list[0].id]),
                    radius: r,
                }];
            }

            // the origin is often inside all surfaces, so do a quick test for that
            let mut origin_inside_radius = self.e_cm_squared;
            let origin = [LorentzVector::default(); MAX_LOOP];
            for e in &ellipsoid_list {
                let r: f64 = -e.evaluate_surface(&origin[..self.n_loops]);
                if r <= 0. || r * r < 1e-10 * self.e_cm_squared {
                    origin_inside_radius = -1.;
                    break;
                }
                if r < origin_inside_radius {
                    origin_inside_radius = r;
                }
            }

            if origin_inside_radius > 0. {
                if self.settings.general.debug > 1 {
                    println!("Origin inside for all E-surface");
                }

                return vec![DeformationOverlap {
                    deformation_sources: origin[..self.n_loops].to_vec(),
                    excluded_surface_ids: vec![],
                    excluded_surface_indices: vec![],
                    overlap: Some(ellipsoid_list.iter().map(|e| e.id).collect()),
                    radius: origin_inside_radius,
                }];
            }
        }

        // collect basic overlap info for all pairs
        for i in 0..ellipsoid_list.len() * ellipsoid_list.len() {
            self.socp_problem.pair_appears_in_overlap_structurex[i] = false;
            self.socp_problem.pair_overlap_matrix[i] = false;
        }

        let mut ecos_time = 0;
        let mut ball_time = 0;

        for i in 0..ellipsoid_list.len() {
            for j in i + 1..ellipsoid_list.len() {
                // skip the SCS check if foci overlap or if the surfaces are on different loop lines
                let mut num_overlap = 0;
                let mut loop_line_overlap = false;
                for foc1 in &ellipsoid_list[i].focal_points {
                    for foc2 in &ellipsoid_list[j].focal_points {
                        if foc1
                            .sig
                            .iter()
                            .zip(foc2.sig.iter())
                            .all(|(s1, s2)| s1.abs() == s2.abs())
                        {
                            loop_line_overlap = true;
                        }
                        if foc1.id == foc2.id {
                            num_overlap += 1;
                        }
                    }
                }

                let mut focus_overlap = false;
                if !loop_line_overlap
                    || (num_overlap + 1 >= ellipsoid_list[i].foci.len()
                        && num_overlap + 1 >= ellipsoid_list[j].foci.len())
                {
                    self.socp_problem.pair_overlap_matrix[i * ellipsoid_list.len() + j] = true;
                    self.socp_problem.pair_overlap_matrix[j * ellipsoid_list.len() + i] = true;

                    focus_overlap = true;
                    if self.settings.general.debug < 2 {
                        continue;
                    }
                    pair_overlap_count += 1;
                }

                let ti = Instant::now();
                let has_overlap_heuristic =
                    self.overlap_heuristic(&ellipsoid_list[i], &ellipsoid_list[j]);
                ball_time += Instant::now().duration_since(ti).as_nanos();

                if self.settings.general.debug < 2 && has_overlap_heuristic {
                    self.socp_problem.pair_overlap_matrix[i * ellipsoid_list.len() + j] = true;
                    self.socp_problem.pair_overlap_matrix[j * ellipsoid_list.len() + i] = true;
                    pair_overlap_count += 1;
                    continue;
                }

                self.socp_problem
                    .construct(&[i, j], &ellipsoid_list, &self.propagators, false);

                // perform the ECOS test
                let ti = Instant::now();
                self.socp_problem.initialize_workspace_ecos();
                let r = self.socp_problem.solve_ecos();
                let has_overlap = r == 0 || r == -2 || r == 10;
                self.socp_problem.finish_ecos();
                ecos_time += Instant::now().duration_since(ti).as_nanos();

                if focus_overlap && !has_overlap {
                    panic!(
                        "Foci overlap but SCS says that there is no overlap for {:#?} and {:#?}",
                        ellipsoid_list[i], ellipsoid_list[j]
                    );
                }

                if has_overlap != has_overlap_heuristic {
                    if self.settings.general.debug > 2 {
                        println!(
                            "Heuristic failed: {} vs {} for {} and {}",
                            has_overlap, has_overlap_heuristic, i, j
                        );
                    }
                }

                // verify with SCS
                if self.settings.general.debug > 1 {
                    self.socp_problem.initialize_workspace_scs();
                    let has_overlap_scs = self.socp_problem.solve_scs() > 0;
                    self.socp_problem.finish_scs();

                    if has_overlap != has_overlap_scs {
                        panic!(
                            "Inconsistency between SCS and ECOS: {} vs {} for {:#?} and {:#?}",
                            has_overlap_scs, has_overlap, ellipsoid_list[i], ellipsoid_list[j],
                        );
                    }
                }

                if has_overlap {
                    self.socp_problem.pair_overlap_matrix[i * ellipsoid_list.len() + j] = true;
                    self.socp_problem.pair_overlap_matrix[j * ellipsoid_list.len() + i] = true;
                    pair_overlap_count += 1;
                }

                problem_count += 1;
            }
        }

        if self.settings.general.debug > 0 {
            println!("Time taken = {:#?}", Instant::now().duration_since(t));
            println!("Heuristic time: {}ms", ball_time as f64 / 1000000.);
            println!("ECOS time:      {}ms", ecos_time as f64 / 1000000.);
            println!("Computed {} pair overlaps", problem_count);
            println!("Pair overlap count: {}", pair_overlap_count);
        }

        let mut option = mem::replace(&mut self.socp_problem.option, vec![]);
        let mut option_translated = mem::replace(&mut self.socp_problem.option_translated, vec![]);
        let mut different = mem::replace(&mut self.socp_problem.different, vec![]);

        for n in (1..=ellipsoid_list.len()).rev() {
            if self.settings.general.debug > 0 {
                println!(
                    "Progress: n={}, current structure={:?}",
                    n, overlap_structure
                );
            }

            let mut initialize = true;
            while self.socp_problem.next_set(
                ellipsoid_list.len(),
                &mut option[..n],
                &mut different[..n],
                initialize,
            ) {
                initialize = false;

                for (ot, o) in option_translated[..n].iter_mut().zip(&option[..n]) {
                    *ot = *o as usize;
                }

                // try if a point inside is inside all other surfaces
                let mut best_point: Option<(usize, f64)> = None;
                'on: for &o1 in &option_translated[..n] {
                    let mut radius = self.e_cm_squared;
                    for &o2 in &option_translated[..n] {
                        if o1 == o2 {
                            continue;
                        }

                        let r = -ellipsoid_list[o2]
                            .evaluate_surface(&ellipsoid_list[o1].interior_point[..self.n_loops]);

                        if r <= 0. || r * r < 1e-10 * self.e_cm_squared {
                            continue 'on;
                        }

                        radius = radius.min(r);
                    }

                    if best_point.is_none() || best_point.unwrap().1 < radius {
                        best_point = Some((o1, radius));
                    }

                    // optimize the internal point if we use heuristical centers
                    if !self.settings.deformation.fixed.use_heuristic_centers {
                        break;
                    }
                }

                let mut has_overlap = best_point.is_some();
                if best_point.is_none() || self.settings.general.debug > 1 {
                    self.socp_problem.construct(
                        &option_translated[..n],
                        &ellipsoid_list,
                        &self.propagators,
                        false,
                    );

                    // perform the ECOS test
                    let ti = Instant::now();
                    self.socp_problem.initialize_workspace_ecos();
                    let r = self.socp_problem.solve_ecos();
                    has_overlap = r == 0 || r == -2 || r == 10;

                    if best_point.is_some() && !has_overlap {
                        panic!(
                            "ECOS claims no overlap, but there is a point inside all for {:?}",
                            &option_translated[..n]
                        );
                    }

                    self.socp_problem.finish_ecos();
                    ecos_time += Instant::now().duration_since(ti).as_nanos();

                    // verify with SCS
                    if self.settings.general.debug > 1 {
                        self.socp_problem.initialize_workspace_scs();
                        let has_overlap_scs = self.socp_problem.solve_scs() > 0;
                        self.socp_problem.finish_scs();

                        if has_overlap != has_overlap_scs {
                            panic!(
                                "Inconsistency between SCS and ECOS: {} vs {} for {:?}",
                                has_overlap_scs,
                                has_overlap,
                                &option_translated[..n]
                            );
                        }
                    }

                    problem_count += 1;
                }

                if has_overlap {
                    let overlap = if self.settings.deformation.fixed.use_heuristic_centers
                        && best_point.is_some()
                    {
                        let (point_e_id, radius) = best_point.unwrap();
                        if self.settings.general.debug > 0 {
                            println!(
                                "Using internal point of E-surface {:?}: {}",
                                point_e_id, radius
                            );
                        }

                        let mut excluded_surface_indices = Vec::with_capacity(ellipsoid_list.len());
                        for (i, e) in ellipsoid_list.iter().enumerate() {
                            if !option_translated[..n].contains(&i) {
                                excluded_surface_indices.push(e.id);
                            }
                        }

                        let overlap = option_translated[..n]
                            .iter()
                            .map(|e_index| ellipsoid_list[*e_index].id)
                            .collect();

                        DeformationOverlap {
                            deformation_sources: ellipsoid_list[point_e_id].interior_point
                                [..self.n_loops]
                                .to_vec(),
                            excluded_surface_ids: vec![],
                            excluded_surface_indices,
                            overlap: Some(overlap),
                            radius,
                        }
                    } else {
                        // find centers with maximal radius
                        if self.settings.deformation.fixed.maximize_radius {
                            let ti = Instant::now();

                            self.socp_problem.construct(
                                &option_translated[..n],
                                &ellipsoid_list,
                                &self.propagators,
                                true,
                            );

                            self.socp_problem.initialize_workspace_ecos();
                            let r = self.socp_problem.solve_ecos();
                            if !(r == 0 || r == -2 || r == 10) {
                                println!(
                                    "ECOS cannot opimize radius: r={}. Falling back to default overlap finding.",
                                    r
                                );

                                self.socp_problem.finish_ecos();
                                self.socp_problem.construct(
                                    &option_translated[..n],
                                    &ellipsoid_list,
                                    &self.propagators,
                                    false,
                                );
                                self.socp_problem.initialize_workspace_ecos();
                                let r = self.socp_problem.solve_ecos();
                                assert!(r == 0 || r == -2 || r == 10, "ECOS cannot find overlap");
                            }
                            self.socp_problem.finish_ecos();
                            ecos_time += Instant::now().duration_since(ti).as_nanos();
                            problem_count += 1;
                        } else if best_point.is_some() {
                            let ti = Instant::now();
                            // no ECOS run has been performed yet
                            self.socp_problem.construct(
                                &option_translated[..n],
                                &ellipsoid_list,
                                &self.propagators,
                                false,
                            );
                            self.socp_problem.initialize_workspace_ecos();
                            let r = self.socp_problem.solve_ecos();
                            self.socp_problem.finish_ecos();
                            ecos_time += Instant::now().duration_since(ti).as_nanos();
                            problem_count += 1;
                            assert!(r == 0 || r == -2 || r == 10, "ECOS cannot find overlap");
                        }

                        self.socp_problem.build_deformation_overlap_info(
                            &option_translated[..n],
                            &ellipsoid_list,
                        )
                    };

                    for (i, &ei) in option[..n].iter().enumerate() {
                        for &ej in option[i + 1..n].iter() {
                            self.socp_problem.pair_appears_in_overlap_structurex
                                [ei as usize * ellipsoid_list.len() + ej as usize] = true;
                            self.socp_problem.pair_appears_in_overlap_structurex
                                [ej as usize * ellipsoid_list.len() + ei as usize] = true;
                        }
                    }

                    overlap_structure.push(overlap);
                }
            }
        }

        mem::swap(&mut self.socp_problem.option, &mut option);
        mem::swap(
            &mut self.socp_problem.option_translated,
            &mut option_translated,
        );
        mem::swap(&mut self.socp_problem.different, &mut different);

        if self.settings.general.debug > 0 {
            println!("Problem count  {}", problem_count);
            println!("Solved in      {:#?}", Instant::now().duration_since(t));
            println!("ECOS time      {:#?}ns", ecos_time);
        }

        overlap_structure
    }
}

/// An SOCP problem instance with cached data for SCS and ECOS. The user should never
/// overwrite the vectors with other vectors, but should update them instead.
/// Resizing vectors is also not allowed.
#[derive(Debug)]
pub struct SOCPProblem {
    pub max_ellipsoids: usize,
    pub max_foci: usize,
    pub n_loops: usize,
    pub radius_computation: bool,
    pub pair_overlap_matrix: Vec<bool>,
    pub pair_appears_in_overlap_structurex: Vec<bool>,
    pub focus_list: Vec<(usize, usize, usize, usize)>,
    pub option: Vec<isize>,
    pub option_translated: Vec<usize>,
    pub different: Vec<bool>,
    pub a_dense: Vec<f64>,
    pub var_map: Vec<usize>,
    pub q: Vec<i32>,
    pub x: Vec<f64>,
    pub i: Vec<scs::scs_int>,
    pub p: Vec<scs::scs_int>,
    pub i_ecos: Vec<i64>,
    pub p_ecos: Vec<i64>,
    pub q_ecos: Vec<i64>,
    pub b: Vec<f64>,
    pub c: Vec<f64>,
    pub sol_x: Vec<f64>,
    pub sol_y: Vec<f64>,
    pub sol_s: Vec<f64>,
    a: Box<scs::ScsMatrix>, // we need to box so that the memory location does not change
    data: Box<scs::ScsData>,
    cone: Box<scs::ScsCone>,
    info: Box<scs::ScsInfo>,
    settings: Box<scs::ScsSettings>,
    sol: Box<scs::ScsSolution>,
    workspace: *mut scs::ScsWork,
    workspace_ecos: *mut ecos::pwork,
}

// This is needed to work with the Python bindings. Generally, this is is not safe.
unsafe impl std::marker::Send for SOCPProblem {}

impl Default for SOCPProblem {
    fn default() -> SOCPProblem {
        SOCPProblem::new(1, 1, 1)
    }
}

impl Clone for SOCPProblem {
    fn clone(&self) -> SOCPProblem {
        let mut p = SOCPProblem::new(self.max_ellipsoids, self.max_foci, self.n_loops);

        // now copy the settings from the old problem
        // all other data is not copied
        p.settings = self.settings.clone();
        p
    }
}

impl SOCPProblem {
    pub fn new(mut max_ellipsoids: usize, mut max_foci: usize, n_loops: usize) -> SOCPProblem {
        let mut settings = Box::new(scs::ScsSettings {
            normalize: 1,
            scale: 1.0,
            rho_x: 1e-3,
            max_iters: 50000,
            eps: 1e-8,
            alpha: 1.5,
            cg_rate: 2.0,
            verbose: 0,
            warm_start: 0, // FIXME: this could cause really bad issues: some tests may fail
            acceleration_lookback: 10,
            write_data_filename: ptr::null(),
        });

        // for the center finding, every ellipsoid and focus is shifted
        max_ellipsoids *= 6 * n_loops;
        max_foci *= 6 * n_loops;

        let pair_overlap_matrix = vec![false; max_ellipsoids * max_ellipsoids];
        let pair_appears_in_overlap_structurex = vec![false; max_ellipsoids * max_ellipsoids];

        let option = vec![0isize; max_ellipsoids];
        let option_translated = vec![0; max_ellipsoids];
        let different = vec![false; max_ellipsoids];

        let num_constraints = max_ellipsoids + max_foci * 5;
        let var_length = 3 * n_loops + max_foci + 1;
        let num_non_empty = max_ellipsoids * n_loops + max_foci + 8 * n_loops * max_foci;

        let focus_list = vec![(0, 0, 0, 0); max_foci];
        let a_dense = vec![0.; num_constraints * var_length];
        let mut q = vec![0; max_foci];
        let q_ecos = vec![0; max_foci];
        let mut x = vec![0.; num_non_empty];
        let mut i = vec![0; num_non_empty];
        let i_ecos = vec![0; num_non_empty];
        let mut p = vec![0; var_length + 1];
        let p_ecos = vec![0; var_length + 1];
        let mut b = vec![0.; num_constraints];
        let mut c = vec![0.; var_length];
        let mut sol_x = vec![0.; var_length];
        let mut sol_y = vec![0.; num_constraints];
        let mut sol_s = vec![0.; num_constraints];

        let mut a = Box::new(scs::ScsMatrix {
            x: &mut x[0] as *mut f64,
            i: &mut i[0] as *mut scs::scs_int,
            p: &mut p[0] as *mut scs::scs_int,
            m: 0,
            n: 0,
        });

        let data = Box::new(scs::ScsData {
            m: 0, // rows
            n: 0, // cols
            A: a.as_mut() as *mut scs::ScsMatrix,
            b: &mut b[0] as *mut f64,
            c: &mut c[0] as *mut f64,
            stgs: settings.as_mut() as *mut scs::ScsSettings,
        });

        let info = Box::new(scs::ScsInfo::default());

        let cone = Box::new(scs::ScsCone {
            f: 0,
            l: 0,
            q: &mut q[0] as *mut scs::scs_int,
            qsize: 0,
            s: ptr::null_mut(),
            ssize: 0,
            ep: 0,
            ed: 0,
            p: ptr::null_mut(),
            psize: 0,
        });

        let sol = Box::new(scs::ScsSolution {
            x: &mut sol_x[0] as *mut f64,
            y: &mut sol_y[0] as *mut f64,
            s: &mut sol_s[0] as *mut f64,
        });

        SOCPProblem {
            max_ellipsoids: max_ellipsoids / 6 / n_loops,
            max_foci: max_foci / 6 / n_loops,
            n_loops,
            radius_computation: false,
            pair_overlap_matrix,
            pair_appears_in_overlap_structurex,
            focus_list,
            option,
            option_translated,
            different,
            a_dense,
            var_map: vec![1000; MAX_LOOP],
            q,
            q_ecos,
            x,
            i,
            i_ecos,
            p,
            p_ecos,
            b,
            c,
            sol_x,
            sol_y,
            sol_s,
            a,
            data,
            cone,
            info,
            settings,
            sol,
            workspace: ptr::null_mut(),
            workspace_ecos: ptr::null_mut(),
        }
    }

    /// With all the variables set
    #[inline]
    pub fn initialize_workspace_ecos(&mut self) {
        for (ie, i) in self.i_ecos.iter_mut().zip(&self.i) {
            *ie = *i as i64;
        }

        for (pe, p) in self.p_ecos.iter_mut().zip(&self.p) {
            *pe = *p as i64;
        }

        for (qe, q) in self.q_ecos.iter_mut().zip(&self.q) {
            *qe = *q as i64;
        }

        unsafe {
            self.workspace_ecos = ecos::ECOS_setup(
                self.data.n as i64,
                self.data.m as i64,
                0,
                self.cone.l as i64,
                self.cone.qsize as i64,
                &mut self.q_ecos[0] as *mut i64,
                0,
                self.a.x,
                &mut self.p_ecos[0] as *mut i64,
                &mut self.i_ecos[0] as *mut i64,
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
                self.data.c,
                self.data.b,
                ptr::null_mut(),
            );

            assert_ne!(self.workspace_ecos, ptr::null_mut());

            (*(*self.workspace_ecos).stgs).verbose = 0;
            //(*(*self.workspace_ecos).stgs).maxit = 100;
            //(*(*self.workspace_ecos).stgs).feastol=1e-13;
        }
    }

    #[inline]
    pub fn solve_ecos(&mut self) -> scs::scs_int {
        unsafe {
            let res = ecos::ECOS_solve(self.workspace_ecos) as i32;

            // now copy the values
            for (x, r) in self.sol_x[..self.data.n as usize]
                .iter_mut()
                .zip(slice::from_raw_parts(
                    (*self.workspace_ecos).x,
                    self.data.n as usize,
                ))
            {
                *x = *r;
            }
            res
        }
    }

    #[inline]
    pub fn finish_ecos(&mut self) {
        unsafe {
            ecos::ECOS_cleanup(self.workspace_ecos, 0);
        }
        self.workspace = ptr::null_mut();
    }

    /// With all the variables set
    #[inline]
    pub fn initialize_workspace_scs(&mut self) {
        unsafe {
            self.workspace = scs::scs_init(
                self.data.as_ref() as *const scs::ScsData,
                self.cone.as_ref() as *const scs::ScsCone,
                self.info.as_mut() as *mut scs::ScsInfo,
            );
        }
    }

    #[inline]
    pub fn solve_scs(&mut self) -> scs::scs_int {
        unsafe {
            scs::scs_solve(
                self.workspace,
                self.data.as_ref() as *const scs::ScsData,
                self.cone.as_ref() as *const scs::ScsCone,
                self.sol.as_mut() as *mut scs::ScsSolution,
                self.info.as_mut() as *mut scs::ScsInfo,
            )
        }
    }

    #[inline]
    pub fn solve_full_scs(&mut self) -> scs::scs_int {
        unsafe {
            scs::scs(
                self.data.as_ref() as *const scs::ScsData,
                self.cone.as_ref() as *const scs::ScsCone,
                self.sol.as_mut() as *mut scs::ScsSolution,
                self.info.as_mut() as *mut scs::ScsInfo,
            )
        }
    }

    #[inline]
    pub fn finish_scs(&mut self) {
        unsafe {
            scs::scs_finish(self.workspace);
        }
        self.workspace = ptr::null_mut();
    }
}
