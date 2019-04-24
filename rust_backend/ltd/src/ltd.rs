use arrayvec::ArrayVec;
use disjoint_sets::UnionFind;
use dual_num::{DualN, Scalar};
use float;
use itertools::Itertools;
use num::Complex;
use num_traits::{Float, FloatConst, FromPrimitive, Inv, Num, NumCast, One, Signed, Zero};
use topologies::{CacheSelector, Cut, CutInfo, CutList, LTDCache, LoopLine, Surface, Topology};
use vector::{Field, LorentzVector, RealNumberLike};
use {
    AdditiveMode, DeformationStrategy, OverallDeformationScaling, ParameterizationMapping,
    ParameterizationMode, MAX_LOOP,
};

use utils;

type Dual4<T> = DualN<T, dual_num::U4>;
type Dual7<T> = DualN<T, dual_num::U7>;
type Dual10<T> = DualN<T, dual_num::U10>;
type Dual13<T> = DualN<T, dual_num::U13>;

impl LoopLine {
    /// Get the momenta of the cut with the cut evaluated
    /// `(+/-sqrt(|q_cut^2| + m^2), q_cut_vec)`
    /// where `q_cut` is the cut momentum (including affine term).
    /// The energies are taken from the cache.
    fn get_cut_momentum<T: FromPrimitive + Float + Field>(
        &self,
        loop_momenta: &[LorentzVector<T>],
        cut: &Cut,
        cut_energies: &[T],
    ) -> LorentzVector<T> {
        // construct the loop momentum that flows through this loop line
        let mut mom: LorentzVector<T> = LorentzVector::default();
        for (l, &c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom += l * T::from_i8(c).unwrap();
        }

        match *cut {
            Cut::PositiveCut(j) => {
                let q: LorentzVector<T> = self.propagators[j].q.cast();
                mom[0] = cut_energies[self.propagators[j].id];
                mom.x += q.x;
                mom.y += q.y;
                mom.z += q.z;
                mom
            }
            Cut::NegativeCut(j) => {
                let q: LorentzVector<T> = self.propagators[j].q.cast();
                mom[0] = -cut_energies[self.propagators[j].id];
                mom.x += q.x;
                mom.y += q.y;
                mom.z += q.z;
                mom
            }
            Cut::NoCut => unreachable!("Cannot compute the energy part of a non-cut loop line"),
        }
    }

    /// Return the inverse of the evaluated loop line
    fn evaluate<
        T: From<float>
            + 'static
            + Num
            + FromPrimitive
            + Float
            + Field
            + Signed
            + RealNumberLike
            + FloatConst
            + std::fmt::LowerExp
            + num_traits::float::FloatCore,
    >(
        &self,
        loop_momenta: &[LorentzVector<num::Complex<T>>],
        cut: &Cut,
        topo: &Topology,
        cache: &LTDCache<T>,
    ) -> Result<num::Complex<T>, &'static str> {
        let kinematics_scale = T::from_f64(topo.e_cm_squared).unwrap();
        let threshold = T::from_f64(topo.settings.general.numerical_threshold).unwrap();

        // construct the loop energy
        // if the line is cut, we can get it from the cache
        let e = match cut {
            Cut::PositiveCut(i) => {
                cache.complex_cut_energies[self.propagators[*i].id]
                    - T::from_f64(self.propagators[*i].q.t).unwrap()
            }
            Cut::NegativeCut(i) => {
                -cache.complex_cut_energies[self.propagators[*i].id]
                    - T::from_f64(self.propagators[*i].q.t).unwrap()
            }
            Cut::NoCut => {
                let mut e: num::Complex<T> = Complex::default();
                for (l, &c) in loop_momenta.iter().zip(self.signature.iter()) {
                    e += l.t * num::Complex::new(T::from_i8(c).unwrap(), T::zero());
                }
                e
            }
        };

        let mut res = num::Complex::new(T::one(), T::zero());
        if topo.settings.general.debug > 3 {
            println!(
                "Loop line evaluation for cut {}\n  | signature = {:?}",
                cut, self.signature
            );
        }
        for (i, p) in self.propagators.iter().enumerate() {
            match cut {
                Cut::PositiveCut(j) if i == *j => {
                    let r = (e + T::from_f64(p.q.t).unwrap()) * T::from_f64(2.).unwrap();
                    res *= r;

                    if topo.settings.general.debug > 3 {
                        println!("  | prop +{}={}", i, r);
                    }
                }
                Cut::NegativeCut(j) if i == *j => {
                    let r = (e + T::from_f64(p.q.t).unwrap()) * T::from_f64(-2.).unwrap();
                    res *= r;

                    if topo.settings.general.debug > 3 {
                        println!("  | prop -{}={}", i, r);
                    }
                }
                _ => {
                    // multiply dual propagator
                    let r = utils::powi(e + T::from_f64(p.q.t).unwrap(), 2)
                        - cache.complex_prop_spatial[p.id];

                    if topo.settings.general.debug > 3 {
                        println!("  | prop  {}={}", i, r);
                    }

                    if !r.is_finite() || r.re * r.re < kinematics_scale * threshold {
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
                    .zip(self.loop_lines.iter())
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
            .zip(self.ltd_cut_options.iter())
            .enumerate()
        {
            let mut cut_signatures_matrix = vec![];

            // note that negative cuts will be treated later
            for (c, ll) in residue_sign.iter().zip(self.loop_lines.iter()) {
                if *c != 0 {
                    for mom in &ll.signature {
                        cut_signatures_matrix.push(*mom as f64);
                    }
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
                let mut cut_mass_sum = float::zero();
                for (cut, ll) in cut_option.iter().zip(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut
                    {
                        cut_shift.push(ll.propagators[*cut_prop_index].q.cast());
                        cut_mass_sum += float::from_f64(ll.propagators[*cut_prop_index].m_squared)
                            .unwrap()
                            .sqrt();
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
                        for (&sign, q) in sig_ll_in_cb.iter().zip(cut_shift.iter()) {
                            surface_shift -= q * float::from_i8(sign).unwrap();
                        }

                        // now see if we have an ellipsoid
                        // 1. surface_shift != 0
                        // 2. surface_shift^2 - (sum_i m_i)^2 >= 0
                        // 3. all signs need to be the same (except 0)
                        // 4. surface_shift.t needs to have the opposite sign as in step 3.
                        // TODO: check if the ansatz for the massive case holds at more than one loop
                        if surface_shift.t != float::zero() {
                            // multiply the residue sign
                            let mut surface_signs = sig_ll_in_cb.clone();
                            for (ss, sign) in surface_signs.iter_mut().zip(residue_sign) {
                                if *sign < 0 {
                                    *ss *= -1;
                                }
                            }

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
                                        - (cut_mass_sum
                                            + float::from_f64(onshell_prop.m_squared.sqrt())
                                                .unwrap())
                                        .powi(2)
                                        >= float::zero()
                                        && float::from_i8(delta_sign).unwrap() * surface_shift.t
                                            < float::zero()
                                    {
                                        /*
                                        println!(
                                            "Found ellipsoid for {}: cut={:?}, mom_map={:?}, signs={:?}, marker={}, shift={}",
                                            self.name, cut_option, sig_ll_in_cb.iter().cloned().collect::<Vec<_>>(),
                                            surface_signs.iter().cloned().collect::<Vec<_>>(), delta_sign, surface_shift
                                        );
                                        */
                                        self.surfaces.push(Surface {
                                            group,
                                            ellipsoid: true,
                                            cut_structure_index: cut_index,
                                            cut_option_index,
                                            cut: cut_option.clone(),
                                            onshell_ll_index: ll_index,
                                            onshell_prop_index,
                                            delta_sign,
                                            sig_ll_in_cb: sig_ll_in_cb.iter().cloned().collect(),
                                            signs: surface_signs.iter().cloned().collect(),
                                            shift: surface_shift.clone(),
                                        });
                                    }
                                } else {
                                    // TODO: existence condition for hyperboloids only known up to 2 loops
                                    if surface_shift.square()
                                        - (cut_mass_sum
                                            - float::from_f64(onshell_prop.m_squared.sqrt())
                                                .unwrap())
                                        .powi(2)
                                        <= float::zero()
                                    {
                                        self.surfaces.push(Surface {
                                            group,
                                            ellipsoid: false,
                                            cut_structure_index: cut_index,
                                            cut_option_index,
                                            cut: cut_option.clone(),
                                            onshell_ll_index: ll_index,
                                            onshell_prop_index,
                                            delta_sign,
                                            sig_ll_in_cb: sig_ll_in_cb.iter().cloned().collect(),
                                            signs: surface_signs.iter().cloned().collect(),
                                            shift: surface_shift.clone(),
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
        let mut group_representative: Vec<(Vec<((usize, usize), i8)>, usize)> = vec![];
        for (surf_index, s) in &mut self.surfaces.iter_mut().enumerate() {
            // create a tuple that identifies a surface:
            // (p_i, s_i) where p_i is a cut propagator or an on-shell one
            // and s_i is the surface sign
            // create a list of propagators that are cut or are on-shell
            // we sort it with their respective
            let mut s_cut_sorted = s
                .cut
                .iter()
                .enumerate()
                .filter(|(_, x)| **x != Cut::NoCut)
                .zip(s.signs.iter())
                .filter_map(|((lli, c), s)| match c {
                    Cut::PositiveCut(i) | Cut::NegativeCut(i) if *s != 0 => Some(((lli, *i), *s)),
                    _ => None,
                })
                .collect::<Vec<_>>();

            // now add the surface sign
            s_cut_sorted.push(((s.onshell_ll_index, s.onshell_prop_index), s.delta_sign));
            s_cut_sorted.sort();

            // also try with opposite signs
            let s_cut_sorted_inv = s_cut_sorted
                .iter()
                .map(|(c, s)| (*c, -s))
                .collect::<Vec<_>>();

            let mut is_new = false;
            match group_representative.iter().find(|r| r.0 == s_cut_sorted) {
                Some(i) => s.group = i.1,
                None => {
                    match group_representative
                        .iter()
                        .find(|r| r.0 == s_cut_sorted_inv)
                    {
                        Some(j) => s.group = j.1,
                        None => {
                            is_new = true;
                        }
                    }
                }
            }

            if is_new {
                s.group = surf_index;
                group_representative.push((s_cut_sorted, surf_index));
            }

            if self.settings.general.debug > 1 {
                println!(
                    "  | id={}, group={}, ellipsoid={}, prop={:?} cut={}, mom_map={:?}, signs={:?}, marker={}, shift={}",
                    surf_index, s.group, s.ellipsoid, (s.onshell_ll_index, s.onshell_prop_index), CutList(&s.cut), s.sig_ll_in_cb,
                    s.signs, s.delta_sign, s.shift
                );
            }
        }

        if self.settings.general.debug > 1 {
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
                    if !s.ellipsoid || s.group != surf_index {
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
            if s.ellipsoid {
                continue;
            }

            let mut cs = s
                .cut
                .iter()
                .zip(&self.loop_lines)
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
    }

    /// Map a vector in the unit hypercube to the infinite hypercube.
    /// Also compute the Jacobian.
    pub fn parameterize<
        T: From<float>
            + Num
            + FromPrimitive
            + Float
            + Field
            + FloatConst
            + RealNumberLike
            + num_traits::float::FloatCore,
    >(
        &self,
        x: &[f64],
        loop_index: usize,
    ) -> ([T; 3], T) {
        let e_cm = T::from_f64(self.e_cm_squared).unwrap().sqrt()
            * T::from_f64(self.settings.parameterization.shifts[loop_index].0).unwrap();
        let mut l_space = [T::zero(); 3];
        let mut jac = T::one();

        match self.settings.parameterization.mode {
            ParameterizationMode::Cartesian => match self.settings.parameterization.mapping {
                ParameterizationMapping::Log => {
                    for i in 0..3 {
                        let x = T::from_f64(x[i]).unwrap();
                        l_space[i] = e_cm * (x / (T::one() - x)).ln();
                        jac *= e_cm / (x - x * x);
                    }
                }
                ParameterizationMapping::Linear => {
                    for i in 0..3 {
                        let x = T::from_f64(x[i]).unwrap();
                        l_space[i] = e_cm * (T::one() / (T::one() - x) - T::one() / x);
                        jac *= e_cm
                            * (T::one() / (x * x) + T::one() / ((T::one() - x) * (T::one() - x)));
                    }
                }
            },
            ParameterizationMode::Spherical => {
                let e_cm = T::from_f64(self.e_cm_squared).unwrap().sqrt()
                    * T::from_f64(self.settings.parameterization.shifts[loop_index].0).unwrap();
                let radius = match self.settings.parameterization.mapping {
                    ParameterizationMapping::Log => {
                        // ln(1/(1-x))
                        let radius =
                            e_cm * (T::one() / (T::one() - T::from_f64(x[0]).unwrap())).ln();
                        jac *= e_cm / (T::one() - T::from_f64(x[0]).unwrap());
                        radius
                    }
                    ParameterizationMapping::Linear => {
                        // x/(1-x)
                        let radius = e_cm * T::from_f64(x[0]).unwrap()
                            / (T::one() - T::from_f64(x[0]).unwrap());
                        jac *= <T as num_traits::Float>::powi(e_cm + radius, 2) / e_cm;
                        radius
                    }
                };
                let phi =
                    T::from_f64(2.).unwrap() * <T as FloatConst>::PI() * T::from_f64(x[1]).unwrap();
                jac *= T::from_f64(2.).unwrap() * <T as FloatConst>::PI();

                let cos_theta = -T::one() + T::from_f64(2.).unwrap() * T::from_f64(x[2]).unwrap(); // out of range
                jac *= T::from_f64(2.).unwrap();
                let sin_theta = (T::one() - cos_theta * cos_theta).sqrt();

                l_space[0] = radius * sin_theta * phi.cos();
                l_space[1] = radius * sin_theta * phi.sin();
                l_space[2] = radius * cos_theta;

                jac *= radius * radius; // spherical coord
            }
        }

        // add a shift such that k=l is harder to be picked up by integrators such as cuhre
        l_space[0] +=
            e_cm * T::from_f64(self.settings.parameterization.shifts[loop_index].1).unwrap();
        l_space[1] +=
            e_cm * T::from_f64(self.settings.parameterization.shifts[loop_index].2).unwrap();
        l_space[2] +=
            e_cm * T::from_f64(self.settings.parameterization.shifts[loop_index].3).unwrap();

        (l_space, jac)
    }

    pub fn inv_parametrize<
        T: From<float>
            + Num
            + FromPrimitive
            + Float
            + Field
            + FloatConst
            + RealNumberLike
            + num_traits::float::FloatCore,
    >(
        &self,
        mom: &LorentzVector<f64>,
        loop_index: usize,
    ) -> ([T; 3], T) {
        if self.settings.parameterization.mode != ParameterizationMode::Spherical {
            panic!("Inverse mapping is only implemented for spherical coordinates");
        }

        let mut jac = T::one();
        let e_cm = T::from_f64(self.e_cm_squared).unwrap().sqrt()
            * T::from_f64(self.settings.parameterization.shifts[loop_index].0).unwrap();

        let x: T = T::from_f64(mom.x).unwrap()
            - e_cm * T::from_f64(self.settings.parameterization.shifts[loop_index].1).unwrap();
        let y: T = T::from_f64(mom.y).unwrap()
            - e_cm * T::from_f64(self.settings.parameterization.shifts[loop_index].2).unwrap();
        let z: T = T::from_f64(mom.z).unwrap()
            - e_cm * T::from_f64(self.settings.parameterization.shifts[loop_index].3).unwrap();

        let k_r_sq = x * x + y * y + z * z;
        let k_r = k_r_sq.sqrt();

        jac /= <T as num_traits::Float>::powi(e_cm + k_r, 2) / e_cm;

        let x2 = if y < T::zero() {
            T::one() + T::from_f64(0.5).unwrap() * T::FRAC_1_PI() * T::atan2(y, x)
        } else {
            T::from_f64(0.5).unwrap() * T::FRAC_1_PI() * T::atan2(y, x)
        };

        // cover the degenerate case
        if k_r_sq.is_zero() {
            return ([T::zero(), x2, T::zero()], T::zero());
        }

        let x1 = match self.settings.parameterization.mapping {
            ParameterizationMapping::Log => {
                let x1 = T::one() - (-k_r / e_cm).exp();
                jac /= e_cm / (T::one() - x1);
                x1
            }
            ParameterizationMapping::Linear => {
                jac /= <T as num_traits::Float>::powi(e_cm + k_r, 2) / e_cm;
                k_r / (e_cm + k_r)
            }
        };

        let x3 = T::from_f64(0.5).unwrap() * (T::one() + z / k_r);

        jac /= T::from_f64(2.).unwrap() * <T as FloatConst>::PI();
        jac /= T::from_f64(2.).unwrap();
        jac /= k_r * k_r;

        ([x1, x2, x3], jac)
    }

    #[inline]
    fn compute_lambda_factor<T: Num + Float + Field + PartialOrd>(x: T, y: T) -> T {
        // FIXME: not smooth
        if x * NumCast::from(2.).unwrap() < y {
            y * NumCast::from(0.25).unwrap()
        } else if y < NumCast::from(0.).unwrap() {
            x - y * NumCast::from(0.5).unwrap()
        } else {
            x - y * NumCast::from(0.25).unwrap()
        }
    }

    /// Determine the global lambda scaling by going through each propagator
    /// for each cut and taking the minimum of their respective lambdas
    /// Additionally, check for the expansion condition and make sure that
    /// the real part of the cut propagator is positive
    fn determine_lambda<
        U: dual_num::Dim + dual_num::DimName,
        T: From<float>
            + Scalar
            + FromPrimitive
            + Signed
            + Field
            + RealNumberLike
            + Float
            + FloatConst
            + std::fmt::LowerExp,
    >(
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
        // update the cut information with the kappa-dependent information
        let cache = cache.get_cache_mut();
        for ll in &self.loop_lines {
            let mut kappa_cut = LorentzVector::default();
            for (kappa, &sign) in kappas.iter().zip(ll.signature.iter()) {
                kappa_cut += kappa * DualN::from_real(T::from_i8(sign).unwrap());
            }

            for p in &ll.propagators {
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
            }
        }

        let mut lambda_sq = DualN::from_real(T::from_f64(lambda_max).unwrap().powi(2));

        let sigma =
            DualN::from_real(T::from_f64(self.settings.deformation.scaling.softmin_sigma).unwrap());
        let mut smooth_min_num = lambda_sq * (-lambda_sq / sigma).exp();
        let mut smooth_min_den = (-lambda_sq / sigma).exp();

        let mut cut_infos = [&cache.cut_info[0]; MAX_LOOP];

        for (cuts, cb_to_lmb_mat) in self.ltd_cut_options.iter().zip(self.cb_to_lmb_mat.iter()) {
            for cut in cuts {
                let mut index = 0;
                for (c, ll) in cut.iter().zip(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(i) | Cut::PositiveCut(i) = c {
                        cut_infos[index] = &cache.cut_info[ll.propagators[*i].id];
                        index += 1;
                    }
                }

                // we have to make sure that our linear expansion for the deformation vectors is reasonable
                // for that we need lambda < c * (q_i^2^cut+m_i^2)/|kappa_i^cut * q_i^cut|
                for cut_info in &cut_infos[..self.n_loops] {
                    if cut_info.kappa_sq.is_zero() {
                        continue;
                    }

                    if self.settings.deformation.scaling.expansion_check {
                        let lambda_exp = DualN::from_real(
                            T::from_f64(self.settings.deformation.scaling.expansion_threshold)
                                .unwrap(),
                        ) * cut_info.spatial_and_mass_sq
                            / cut_info.kappa_dot_mom.abs(); // note: not holomorphic

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
                        let lambda_disc_sq = cut_info.spatial_and_mass_sq / cut_info.kappa_sq
                            * DualN::from_real(T::from_f64(0.95).unwrap());

                        if sigma.is_zero() {
                            if lambda_disc_sq < lambda_sq {
                                lambda_sq = lambda_disc_sq;
                            }
                        } else {
                            let e = (-lambda_disc_sq / sigma).exp();
                            smooth_min_num += lambda_disc_sq * e;
                            smooth_min_den += e;
                        }
                    }
                }

                if !self.settings.deformation.scaling.cut_propagator_check
                    && !self.settings.deformation.scaling.non_cut_propagator_check
                {
                    continue;
                }

                for (ll_cut, onshell_ll) in cut.iter().zip(self.loop_lines.iter()) {
                    // TODO: this should be precomputed
                    // determine the map from loop line momentum to cut momenta
                    let mut sig_ll_in_cb = [0; MAX_LOOP];
                    for ii in 0..self.n_loops {
                        for jj in 0..self.n_loops {
                            // note we are taking the transpose of cb_to_lmb_mat
                            sig_ll_in_cb[ii] +=
                                cb_to_lmb_mat[jj * self.n_loops + ii] * onshell_ll.signature[jj];
                        }
                    }

                    // multiply the signature by the sign of the cut
                    let mut onshell_signs = [0; MAX_LOOP];
                    let mut index = 0;
                    for (os, &sig) in onshell_signs.iter_mut().zip(sig_ll_in_cb.iter()) {
                        *os = sig;
                    }
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

                    if ellipse_flag == 3 && self.settings.deformation.scaling.skip_hyperboloids {
                        continue;
                    }

                    let k0sq_inv = DualN::from_real(T::one())
                        / cache.cut_info[onshell_ll.propagators[0].id].kappa_sq;

                    // if the kappa is 0, there is no need for rescaling
                    if !k0sq_inv.is_finite() {
                        continue;
                    }

                    // treat the cut part of the surface equation
                    let mut a = DualN::from_real(T::zero());
                    let mut b = DualN::from_real(T::zero());
                    let mut c = DualN::from_real(T::zero());
                    for (cut_info, &sign) in
                        cut_infos[..self.n_loops].iter().zip(onshell_signs.iter())
                    {
                        a += cut_info.a * T::from_i8(sign).unwrap();
                        b += cut_info.b * T::from_i8(sign).unwrap();
                        c += cut_info.c * T::from_i8(sign).unwrap();
                    }

                    // TODO: we are checking surfaces more than once. The ones that depend on only
                    // one loop momentum for example
                    for (prop_index, onshell_prop) in onshell_ll.propagators.iter().enumerate() {
                        let on_shell_info = &cache.cut_info[onshell_prop.id];

                        // if this is the cut propagator, we only need to use the spatial part
                        // the functional form remains the same

                        if *ll_cut == Cut::NegativeCut(prop_index)
                            || *ll_cut == Cut::PositiveCut(prop_index)
                        {
                            // TODO: move to block where we do cut propagator checks
                            if !self.settings.deformation.scaling.cut_propagator_check {
                                continue;
                            }

                            let x = (on_shell_info.kappa_dot_mom * k0sq_inv).powi(2);
                            let y = on_shell_info.spatial_and_mass_sq * k0sq_inv;

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
                        } else {
                            if !self.settings.deformation.scaling.non_cut_propagator_check {
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

                                let a_tot =
                                    (a + on_shell_info.a * sign) * T::from_f64(-0.5).unwrap();
                                let a_tot_inv = a_tot.inv();
                                let b_tot = b + on_shell_info.b * sign;
                                let c_tot = c + on_shell_info.c * sign + on_shell_info.shift.t;
                                let x = (b_tot * a_tot_inv).powi(2) * T::from_f64(0.25).unwrap();
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

        if sigma.is_zero() {
            lambda_sq.sqrt()
        } else {
            (smooth_min_num / smooth_min_den).sqrt()
        }
    }

    /// Construct a deformation vector by going through all the ellipsoids
    fn deform_ellipsoids<
        U: dual_num::Dim + dual_num::DimName,
        T: From<float>
            + Scalar
            + 'static
            + FromPrimitive
            + Float
            + Signed
            + Field
            + RealNumberLike
            + FloatConst,
    >(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut cut_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut deform_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut kappas = [LorentzVector::default(); MAX_LOOP];

        let cache_spec = cache.get_cache_mut();
        let kappa_surf = &mut cache_spec.deform_dirs;
        let inv_surf_prop = &mut cache_spec.ellipsoid_eval;

        // surface equation:
        // m*|sum_i^L a_i q_i^cut + p_vec| + sum_i^L a_i*r_i * |q_i^cut| - p^0 = 0
        // where m=delta_sign a = sig_ll_in_cb, r = residue_sign
        let mut group_counter = 0;
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            // only deform the set of different ellipsoids
            if !surf.ellipsoid || surf.group != surf_index {
                continue;
            }

            // compute v_i = q_i^cut / sqrt(|q_i^cut|^2 + m^2)
            // construct the normalized 3-momenta that flow through the cut propagators
            // the 0th component is to be ignored
            let mut cut_counter = 0;
            for (c, ll) in surf.cut.iter().zip(self.loop_lines.iter()) {
                if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = c {
                    let mut mom: LorentzVector<DualN<T, U>> = LorentzVector::default();
                    for (l, &c) in loop_momenta.iter().zip(ll.signature.iter()) {
                        mom += l * DualN::from_real(T::from_i8(c).unwrap());
                    }
                    let q: LorentzVector<DualN<T, U>> = ll.propagators[*i].q.cast();
                    mom += q;
                    let length: DualN<T, U> = mom.spatial_squared_impr()
                        + T::from_f64(ll.propagators[*i].m_squared).unwrap();
                    mom = mom / length.sqrt();
                    cut_dirs[cut_counter] = mom;
                    cut_counter += 1;
                }
            }

            // compute w=(sum_i^L a_i q_i^cut + p_vec) / sqrt(|sum_i^L a_i q_i^cut + p_vec|^2 + m^2)
            // determine sum_i^L a_i q_i^cut using the loop-momentum basis
            let mut surface_dir = LorentzVector::default();
            for (l, &sig) in loop_momenta
                .iter()
                .zip(self.loop_lines[surf.onshell_ll_index].signature.iter())
            {
                surface_dir += l * DualN::from_real(T::from_i8(sig).unwrap());
            }

            let surface_prop =
                &self.loop_lines[surf.onshell_ll_index].propagators[surf.onshell_prop_index];
            let q: LorentzVector<DualN<T, U>> = surface_prop.q.cast();
            surface_dir += q;

            let length: DualN<T, U> =
                surface_dir.spatial_squared_impr() + T::from_f64(surface_prop.m_squared).unwrap();
            let surface_dir_norm = surface_dir / length.sqrt();

            // compute n_i = (v_i + a_i * w) * abs(a_i)
            // n_i is the normal vector of the ith entry of the cut momentum basis
            // this n_i is 0 when the surface does not depend on the ith cut loop line
            // we update cut_dirs from v_i to n_i
            for (cut_dir, &sign) in cut_dirs[..self.n_loops]
                .iter_mut()
                .zip(surf.sig_ll_in_cb.iter())
            {
                if sign != 0 {
                    *cut_dir += surface_dir_norm * DualN::from_real(T::from_i8(sign).unwrap());
                }
            }

            // convert from cut momentum basis to loop momentum basis
            for (i, deform_dir) in deform_dirs[..self.n_loops].iter_mut().enumerate() {
                *deform_dir = LorentzVector::default();
                for (&sign, cut_dir) in self.cb_to_lmb_mat[surf.cut_structure_index]
                    [i * self.n_loops..(i + 1) * self.n_loops]
                    .iter()
                    .zip(cut_dirs[..self.n_loops].iter())
                {
                    *deform_dir += cut_dir * DualN::from_real(T::from_i8(sign).unwrap());
                }
            }

            // multiply the direction deform_dirs by the exponential function
            // for the exponential function we need to evaluate the propagator of the ellipsoid
            // and for that we need to determine the cut energies

            // compute the cut momenta
            let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
            let mut index = 0;
            for (ll_cut, ll) in surf.cut.iter().zip(self.loop_lines.iter()) {
                if *ll_cut != Cut::NoCut {
                    cut_momenta[index] =
                        ll.get_cut_momentum(loop_momenta, ll_cut, &cache_spec.cut_energies);
                    index += 1;
                }
            }

            // evaluate the inverse propagator of the surface
            // the momentum map from the cut momenta is in signs and the shift is known as well
            // compute the energies for the loop momenta
            let mut mom: LorentzVector<DualN<T, U>> = LorentzVector::new();
            for (&s, c) in surf.sig_ll_in_cb.iter().zip(cut_momenta.iter()) {
                mom += c * DualN::from_real(T::from_i8(s).unwrap());
            }

            let qs: LorentzVector<DualN<T, U>> = surf.shift.cast();
            let momq: LorentzVector<DualN<T, U>> = mom + qs;
            let inv = momq.square_impr() - T::from_f64(surface_prop.m_squared).unwrap();
            inv_surf_prop[group_counter] = inv;

            for (loop_index, dir) in deform_dirs.iter().enumerate() {
                // note the sign
                kappa_surf[group_counter * MAX_LOOP + loop_index] =
                    -dir / dir.spatial_squared_impr().sqrt();
            }

            group_counter += 1;
        }

        // now combine the kappas from the surface using the chosen strategy
        match self.settings.general.deformation_strategy {
            DeformationStrategy::Additive => {
                let aij: DualN<T, U> =
                    NumCast::from(self.settings.deformation.additive.a_ij).unwrap();

                for (i, &inv) in inv_surf_prop[..group_counter].iter().enumerate() {
                    let dampening = match self.settings.deformation.additive.mode {
                        AdditiveMode::Exponential => (-inv * inv
                            / (aij * T::from_f64(self.e_cm_squared).unwrap().powi(2)))
                        .exp(),
                        AdditiveMode::Hyperbolic => {
                            let t = inv * inv / T::from_f64(self.e_cm_squared).unwrap().powi(2);
                            t / (t + aij)
                        }
                        AdditiveMode::Unity => DualN::from_real(T::one()),
                        AdditiveMode::SoftMin => unimplemented!(),
                    };

                    for (loop_index, kappa) in kappas[..self.n_loops].iter_mut().enumerate() {
                        let dir = kappa_surf[i * MAX_LOOP + loop_index];
                        if self.settings.general.debug > 2 {
                            println!(
                                "Deformation contribution for surface {}:\n  | dir={}\n  | exp={}",
                                group_counter,
                                dir,
                                dampening.real()
                            );
                        }

                        *kappa += dir * dampening;
                    }
                }
            }
            DeformationStrategy::Multiplicative => {
                let m_ij: DualN<T, U> =
                    NumCast::from(self.settings.deformation.multiplicative.m_ij).unwrap();

                for (i, _) in inv_surf_prop[..group_counter].iter().enumerate() {
                    let mut dampening = DualN::from_real(T::one());
                    for (j, &inv) in inv_surf_prop[..group_counter].iter().enumerate() {
                        if i == j {
                            continue;
                        }

                        let t = inv * inv / T::from_f64(self.e_cm_squared).unwrap().powi(2);
                        dampening *= t / (t + m_ij);
                    }

                    for (loop_index, kappa) in kappas[..self.n_loops].iter_mut().enumerate() {
                        let dir = kappa_surf[i * MAX_LOOP + loop_index];
                        if self.settings.general.debug > 2 {
                            println!(
                                "Deformation contribution for surface {}:\n  | dir={}\n  | exp={}",
                                group_counter,
                                dir,
                                dampening.real()
                            );
                        }

                        *kappa += dir * dampening;
                    }
                }
            }
            _ => unreachable!(),
        }
        kappas
    }

    fn get_deformation_for_cut<
        U: dual_num::Dim + dual_num::DimName,
        T: From<float>
            + 'static
            + FromPrimitive
            + Float
            + Field
            + RealNumberLike
            + num_traits::Signed
            + FloatConst,
    >(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
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

        // compute each cut momentum, subtract the shift and normalize
        let mut index = 0;
        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut cut_masses = [DualN::default(); MAX_LOOP];
        for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
            if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = *ll_cut {
                cut_momenta[index] =
                    ll.get_cut_momentum(loop_momenta, ll_cut, &cache.get_cache().cut_energies);

                // subtract the shift from q0
                cut_momenta[index].t -=
                    DualN::from_real(T::from_f64(ll.propagators[i].q.t).unwrap());
                cut_masses[index] =
                    DualN::from_real(T::from_f64(ll.propagators[i].m_squared).unwrap());

                index += 1;
            }
        }

        // convert from cut momentum basis to loop momentum basis
        for i in 0..self.n_loops {
            let mut deform_dir = LorentzVector::default();
            for ((&sign, cut_dir), &cut_mass) in self.cb_to_lmb_mat[cut_structure_index]
                [i * self.n_loops..(i + 1) * self.n_loops]
                .iter()
                .zip(cut_momenta.iter())
                .zip(cut_masses.iter())
            {
                // normalize the spatial components
                let mut mom_normalized = *cut_dir;
                let length = (mom_normalized.spatial_squared_impr() + cut_mass).sqrt();
                mom_normalized *= length.inv();

                deform_dir += mom_normalized * DualN::from_real(T::from_i8(sign).unwrap());
            }

            deform_dirs[i] = deform_dir;
        }

        deform_dirs
    }

    /// Construct a deformation vector by going through all cut options
    /// and make sure it is 0 on all other cut options
    fn deform_cutgroups<
        U: dual_num::Dim + dual_num::DimName,
        T: From<float>
            + Scalar
            + 'static
            + FromPrimitive
            + Float
            + Field
            + RealNumberLike
            + std::fmt::LowerExp
            + num_traits::Signed
            + num_traits::float::FloatCore
            + FloatConst,
    >(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];

        // go through all cuts that contain ellipsoids
        // by going through all ellipsoids and skipping ones that are in the same cut
        // TODO: precompute non-empty cuts
        let mut last_cut = None;
        let mut non_empty_cut_count = 0;
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            if !surf.ellipsoid {
                continue;
            }

            // evaluate ellipsoid representatives
            if surf.group == surf_index {
                // FIXME: we compute this again later
                let mut index = 0;

                // calculate the cut momenta
                let mut res = DualN::from_real(T::zero());
                for (&ll_cut, ll) in surf.cut.iter().zip(self.loop_lines.iter()) {
                    if ll_cut != Cut::NoCut {
                        cut_momenta[index] = ll.get_cut_momentum(
                            loop_momenta,
                            &ll_cut,
                            &cache.get_cache().cut_energies,
                        );
                        index += 1;
                    }
                }

                let mut normalization = DualN::from_real(T::zero());
                let mut q3 = LorentzVector::default();
                for (&s, mom) in surf.sig_ll_in_cb.iter().zip(&cut_momenta) {
                    if s != 0 {
                        q3 += mom * DualN::from_real(T::from_i8(s).unwrap());
                        res += mom.t.abs();
                        normalization += mom.square();
                    }
                }
                let shift: LorentzVector<DualN<T, U>> = surf.shift.cast();
                q3 += shift;

                res += q3.spatial_squared_impr().sqrt();
                res -= shift.t.abs();
                normalization += shift.square();
                cache.get_cache_mut().ellipsoid_eval[surf_index] = res * res / normalization;
            }

            if Some(&surf.cut) != last_cut {
                let mut deform_dirs_cut = self.get_deformation_for_cut(
                    loop_momenta,
                    &surf.cut,
                    surf.cut_structure_index,
                    cache,
                );

                for i in 0..self.n_loops {
                    cache.get_cache_mut().deform_dirs[non_empty_cut_count * self.n_loops + i] =
                        deform_dirs_cut[i];
                }

                cache.get_cache_mut().non_empty_cuts[non_empty_cut_count] =
                    (surf.cut_structure_index, surf.cut_option_index);

                last_cut = Some(&surf.cut);
                non_empty_cut_count += 1;
            }
        }

        let cache = cache.get_cache();
        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        for (i, &(cut_structure_index, cut_option_index)) in cache.non_empty_cuts
            [..non_empty_cut_count]
            .iter()
            .enumerate()
        {
            if self.settings.general.debug > 2 {
                println!(
                    "Surface contributions for cut {}:",
                    CutList(&self.ltd_cut_options[cut_structure_index][cut_option_index])
                );
            }

            let mut s = DualN::from_real(T::one());
            let mut softmin_num = DualN::from_real(T::zero());
            let mut softmin_den = DualN::from_real(T::zero());

            if self.settings.deformation.cutgroups.mode != AdditiveMode::Unity {
                for &surf_index in
                    &self.ellipsoids_not_in_cuts[cut_structure_index][cut_option_index]
                {
                    // t is the weighing factor that is 0 if we are on both cut_i and cut_j
                    // at the same time and goes to 1 otherwise
                    let t = cache.ellipsoid_eval[surf_index];

                    let sup = if self.settings.deformation.cutgroups.mode == AdditiveMode::SoftMin {
                        if self.settings.deformation.cutgroups.sigma.is_zero() {
                            s = s.min(t);
                            s
                        } else {
                            let e = (-t
                                / T::from_f64(self.settings.deformation.cutgroups.sigma).unwrap())
                            .exp();
                            softmin_num += t * e;
                            softmin_den += e;
                            e
                        }
                    } else {
                        let sup = t
                            / (t + DualN::from_real(
                                T::from_f64(self.settings.deformation.cutgroups.m_ij).unwrap(),
                            ));
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
                }

                if self.settings.deformation.cutgroups.mode == AdditiveMode::SoftMin {
                    if !self.settings.deformation.cutgroups.sigma.is_zero() {
                        s = softmin_num / softmin_den;
                    }

                    if !self.settings.deformation.cutgroups.m_ij.is_zero() {
                        s = s
                            / (s + DualN::from_real(
                                T::from_f64(self.settings.deformation.cutgroups.m_ij).unwrap(),
                            ));
                    }
                }
            }

            if self.settings.general.debug > 2 {
                println!(
                    "  | k={:e}\n  | dirs={:e}\n  | suppression={:e}\n  | contribution={:e}",
                    loop_momenta[0].real(),
                    &cache.deform_dirs[i * self.n_loops].real(),
                    s.real(),
                    &cache.deform_dirs[i * self.n_loops].real() * s.real(),
                );
            }

            for ii in 0..self.n_loops {
                kappas[ii] -= cache.deform_dirs[i * self.n_loops + ii] * s;
                //kappas[ii].t = DualN::from_real(T::zero());
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{}={:e}\n", ii + 1, kappas[ii].real());
            }
        }

        kappas
    }

    /// Construct a deformation vector by going through all the constant
    fn deform_constant<
        U: dual_num::Dim + dual_num::DimName,
        T: From<float>
            + Scalar
            + 'static
            + FromPrimitive
            + Float
            + Signed
            + Field
            + RealNumberLike
            + FloatConst,
    >(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        _cache: &mut LTDCache<T>,
    ) -> [LorentzVector<DualN<T, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy,
        LTDCache<T>: CacheSelector<T, U>,
    {
        let mut kappas = [LorentzVector::default(); MAX_LOOP];

        for i in 0..self.n_loops {
            kappas[i] = loop_momenta[i] / loop_momenta[i].spatial_squared_impr().sqrt();
        }

        kappas
    }

    fn deform_generic<
        U: dual_num::Dim + dual_num::DimName,
        T: From<float>
            + 'static
            + FromPrimitive
            + Float
            + Field
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + FloatConst
            + std::fmt::LowerExp,
    >(
        &self,
        loop_momenta: &[LorentzVector<DualN<T, U>>],
        cut: Option<(usize, usize)>,
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
                for (l, &c) in loop_momenta.iter().zip(ll.signature.iter()) {
                    mom += l * DualN::from_real(T::from_i8(c).unwrap());
                }

                // update the cut info
                for p in &ll.propagators {
                    let q: LorentzVector<DualN<T, U>> = p.q.cast();
                    let cm = (mom + q).spatial_squared() + T::from_f64(p.m_squared).unwrap();
                    let energy = cm.sqrt();
                    cache.cut_energies[p.id] = energy;

                    let info = &mut cache.cut_info[p.id];
                    info.id = p.id;
                    info.momentum = mom + q;
                    info.real_energy = energy;
                    info.spatial_and_mass_sq = cm;
                    // TODO: these two never change and can be set at the start!
                    info.shift = q;
                    info.mass = DualN::from_real(T::from_f64(p.m_squared).unwrap());
                }
            }
        }

        let mut kappas = match self.settings.general.deformation_strategy {
            DeformationStrategy::CutGroups => self.deform_cutgroups(loop_momenta, cache),
            DeformationStrategy::Duals => {
                let co = cut.unwrap();
                let cut = &self.ltd_cut_options[co.0][co.1];
                self.get_deformation_for_cut(loop_momenta, cut, co.0, cache)
            }
            DeformationStrategy::Additive | DeformationStrategy::Multiplicative => {
                self.deform_ellipsoids(loop_momenta, cache)
            }
            DeformationStrategy::Constant => self.deform_constant(loop_momenta, cache),
            DeformationStrategy::None => unreachable!(),
        };

        // make sure the kappa has the right dimension by multiplying in the scale
        let scale = DualN::from_real(
            T::from_f64(
                self.e_cm_squared.sqrt() * self.settings.deformation.overall_scaling_constant,
            )
            .unwrap(),
        );
        for (kappa, k) in &mut kappas[..self.n_loops].iter_mut().zip(loop_momenta.iter()) {
            match self.settings.deformation.overall_scaling {
                OverallDeformationScaling::Constant => *kappa *= scale,
                OverallDeformationScaling::Linear => {
                    *kappa *= k.spatial_squared_impr().sqrt()
                        * T::from_f64(self.settings.deformation.overall_scaling_constant).unwrap()
                }
                OverallDeformationScaling::Sigmoid => {
                    let k_scale = k.spatial_squared_impr().sqrt();
                    *kappa *= k_scale * T::from_f64(2.).unwrap()
                        / (DualN::one() + (k_scale / scale).exp());
                }
            }
        }

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{} scaled={:e}\n", ii + 1, kappas[ii].real());
            }
        }

        let lambda = if self.settings.deformation.scaling.lambda > 0. {
            self.determine_lambda(&kappas, self.settings.deformation.scaling.lambda, cache)
        } else {
            NumCast::from(self.settings.deformation.scaling.lambda.abs()).unwrap()
        };

        if self.settings.general.debug > 2 {
            println!("  | lambda={:e}\n", lambda.real());
        }

        for k in kappas[..self.n_loops].iter_mut() {
            *k *= lambda;
            k.t = DualN::from_real(T::zero()); // make sure we do not have a left-over deformation
        }

        let jac_mat = &mut cache.get_cache_mut().deformation_jacobian;
        for i in 0..3 * self.n_loops {
            for j in 0..3 * self.n_loops {
                // first index: loop momentum, second: xyz, third: dual
                jac_mat[i * 3 * self.n_loops + j] =
                    Complex::new(T::zero(), kappas[i / 3][i % 3 + 1][j + 1]);
            }
            jac_mat[i * 3 * self.n_loops + i] += Complex::new(T::one(), T::zero());
        }

        let jac = utils::determinant(jac_mat, 3 * self.n_loops);
        // take the real part
        let mut r = [LorentzVector::default(); MAX_LOOP];
        for (rr, k) in r.iter_mut().zip(&kappas[..self.n_loops]) {
            *rr = k.real();
        }
        (r, jac)
    }

    pub fn deform<
        T: From<float>
            + FromPrimitive
            + 'static
            + Float
            + Field
            + RealNumberLike
            + num_traits::Signed
            + FloatConst
            + std::fmt::LowerExp
            + num_traits::float::FloatCore,
    >(
        &self,
        loop_momenta: &[LorentzVector<T>],
        cut: Option<(usize, usize)>,
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

                return self.deform_generic(&r, cut, cache);
            }
            2 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual7::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual7::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = T::one();
                    r[1][i + 1][i + 4] = T::one();
                }
                self.deform_generic(&r, cut, cache)
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
                self.deform_generic(&r, cut, cache)
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
                self.deform_generic(&r, cut, cache)
            }
            n => panic!("Binding for deformation at {} loops is not implemented", n),
        }
    }

    /// Set the energy component of the loop momenta according to
    /// `cut`. It takes the cut energies from the cache.
    #[inline]
    pub fn set_loop_momentum_energies<
        T: From<float>
            + Signed
            + 'static
            + Num
            + FromPrimitive
            + Float
            + Field
            + RealNumberLike
            + FloatConst,
    >(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &LTDCache<T>,
    ) {
        // compute the cut energy for each loop line
        let mut cut_energy = [Complex::default(); MAX_LOOP];
        let mut index = 0;
        for (&ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
            if let Cut::PositiveCut(j) | Cut::NegativeCut(j) = ll_cut {
                let e = cache.complex_cut_energies[ll.propagators[j].id];
                if let Cut::PositiveCut(_) = ll_cut {
                    cut_energy[index] =
                        e - Complex::new(T::from_f64(ll.propagators[j].q.t).unwrap(), T::zero());
                } else {
                    cut_energy[index] =
                        -e - Complex::new(T::from_f64(ll.propagators[j].q.t).unwrap(), T::zero());
                }
                index += 1;
            }
        }

        // compute the energies for the loop momenta
        for (i, l) in k_def.iter_mut().enumerate() {
            l.t = Complex::default();
            for (c, e) in mat[i * self.n_loops..(i + 1) * self.n_loops]
                .iter()
                .zip(&cut_energy[..self.n_loops])
            {
                l.t += e * T::from_i8(*c).unwrap();
            }
        }
    }

    #[inline]
    pub fn evaluate_cut<
        T: From<float>
            + Num
            + FromPrimitive
            + Float
            + Field
            + RealNumberLike
            + num_traits::Signed
            + FloatConst
            + std::fmt::LowerExp
            + num_traits::float::FloatCore
            + 'static,
    >(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &str> {
        self.set_loop_momentum_energies(k_def, cut, mat, cache);

        let mut r = Complex::new(T::one(), T::zero());
        for (i, (ll_cut, ll)) in cut.iter().zip(self.loop_lines.iter()).enumerate() {
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
                        k_def[0].t
                    );
                }
                2 => {
                    println!(
                        "Cut {}:\n  | result={:e}\n  | k0={:e}\n  | l0={:e}",
                        CutList(cut),
                        r,
                        k_def[0].t,
                        k_def[1].t
                    );
                }
                3 => {
                    println!(
                        "Cut {}:\n  | result={:e}\n  | k0={:e}\n  | l0={:e}\n  | m0={:e}",
                        CutList(cut),
                        r,
                        k_def[0].t,
                        k_def[1].t,
                        k_def[2].t
                    );
                }
                _ => {}
            }
        }

        Ok(r)
    }

    pub fn compute_complex_cut_energies<
        T: From<float>
            + Num
            + FromPrimitive
            + Float
            + Field
            + FloatConst
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + FloatConst
            + std::fmt::LowerExp
            + 'static,
    >(
        &self,
        k_def: &ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        cache: &mut LTDCache<T>,
    ) {
        // compute all complex cut energies
        for (ll_index, ll) in self.loop_lines.iter().enumerate() {
            let mut mom = LorentzVector::default();
            for (l, &c) in k_def.iter().zip(ll.signature.iter()) {
                mom += l * num::Complex::new(T::from_i8(c).unwrap(), T::zero());
            }

            for p in &ll.propagators {
                let q: LorentzVector<Complex<T>> = p.q.cast();
                let cm = (mom + q).spatial_squared() + T::from_f64(p.m_squared).unwrap();

                if self.settings.deformation.scaling.positive_cut_check
                    && cm.re < T::zero()
                    && cm.im < T::zero()
                {
                    eprintln!(
                        "Branch cut detected for prop {}, ll sig={:?}, ks={:?}: {}",
                        p.id, ll.signature, k_def, cm
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
                    cache.complex_loop_line_eval[ll_index][i][0] = ll
                        .evaluate(k_def, &Cut::PositiveCut(i), &self, cache)
                        .unwrap();
                }
                if has_negative_cut {
                    cache.complex_loop_line_eval[ll_index][i][1] = ll
                        .evaluate(k_def, &Cut::NegativeCut(i), &self, cache)
                        .unwrap();
                }
            }
        }
    }

    #[inline]
    pub fn evaluate<
        'a,
        T: From<float>
            + Num
            + FromPrimitive
            + Float
            + Field
            + FloatConst
            + RealNumberLike
            + num_traits::Signed
            + num_traits::float::FloatCore
            + FloatConst
            + std::fmt::LowerExp
            + 'static,
    >(
        &self,
        x: &'a [f64],
        cache: &mut LTDCache<T>,
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
            let (mut l_space, jac) = self.parameterize(&x[i * 3..(i + 1) * 3], i);

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
        let mut ct = Complex::zero();

        if self.settings.general.deformation_strategy != DeformationStrategy::Duals {
            let (kappas, jac) = self.deform(&k, None, cache);
            k_def = (0..self.n_loops)
                .map(|i| {
                    k[i].map(|x| Complex::new(x, T::zero()))
                        + kappas[i].map(|x| Complex::new(T::zero(), x))
                })
                .collect();
            jac_def = jac;
            ct = self.counterterm(&k_def);

            self.compute_complex_cut_energies(&k_def, cache);
        }

        // evaluate all dual integrands
        let mut result = Complex::default();
        let mut cut_counter = 0;
        for (cut_structure_index, (cuts, mat)) in self
            .ltd_cut_options
            .iter()
            .zip(self.cb_to_lmb_mat.iter())
            .enumerate()
        {
            // for each cut coming from the same cut structure
            for (cut_option_index, cut) in cuts.iter().enumerate() {
                if self.settings.general.cut_filter.is_empty()
                    || self.settings.general.cut_filter.contains(&cut_counter)
                {
                    let mut dual_jac_def = Complex::one();
                    if self.settings.general.deformation_strategy == DeformationStrategy::Duals {
                        let (kappas, jac) =
                            self.deform(&k, Some((cut_structure_index, cut_option_index)), cache);
                        k_def = (0..self.n_loops)
                            .map(|i| {
                                k[i].map(|x| Complex::new(x, T::zero()))
                                    + kappas[i].map(|x| Complex::new(T::zero(), x))
                            })
                            .collect();
                        dual_jac_def = jac;
                        ct = self.counterterm(&k_def);
                        self.compute_complex_cut_energies(&k_def, cache);
                    }

                    match self.evaluate_cut(&mut k_def, cut, mat, cache) {
                        Ok(v) => {
                            // Regulate each sub diagram individually
                            result += v * (ct + T::one()) * dual_jac_def
                        }
                        Err(_) => return (x, k_def, jac_para, jac_def, Complex::default()),
                    }
                }

                cut_counter += 1;
            }
        }

        result *= utils::powi(
            num::Complex::new(
                T::zero(),
                T::from_f64(-2.).unwrap() * <T as FloatConst>::PI(),
            ),
            self.n_loops,
        ); // factor of delta cut

        result *= utils::powi(
            num::Complex::new(
                T::from_f64(1.).unwrap()
                    / <T as Float>::powi(T::from_f64(2.).unwrap() * <T as FloatConst>::PI(), 4),
                T::zero(),
            ),
            self.n_loops,
        ); // loop momentum factor
        result *= jac_def * jac_para;

        (x, k_def, jac_para, jac_def, result)
    }
}
