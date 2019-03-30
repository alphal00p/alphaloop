use arrayvec::ArrayVec;
use dual_num::DualN;
use float;
use itertools::Itertools;
use num_traits::FloatConst;
use num_traits::FromPrimitive;
use num_traits::One;
use num_traits::Zero;
use num_traits::{Float, Inv, Num, NumCast};
use topologies::{Cut, CutList, LoopLine, Surface, Topology};
use vector::{Field, LorentzVector};
use {AdditiveMode, DeformationStrategy};

use utils;

const MAX_ELLIPSE_GROUPS: usize = 100;
const MAX_DIM: usize = 3;
const MAX_LOOP: usize = 3;

type Dual4 = DualN<float, dual_num::U4>;
type Dual7 = DualN<float, dual_num::U7>;
type Dual10 = DualN<float, dual_num::U10>;
type Complex = num::Complex<float>;

impl LoopLine {
    /// Get the momenta of the cut with the cut evaluated
    /// `(+/-sqrt(|q_cut^2| + m^2), q_cut_vec)`
    /// where `q_cut` is the cut momentum (including affine term).
    fn get_cut_momentum<T: From<float> + Num + FromPrimitive + Float + Field>(
        &self,
        loop_momenta: &[LorentzVector<T>],
        cut: &Cut,
    ) -> LorentzVector<T> {
        // construct the loop momentum that flows through this loop line
        let mut mom: LorentzVector<T> = LorentzVector::default();
        for (l, &c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom += l * T::from_i8(c).unwrap();
        }

        match *cut {
            Cut::PositiveCut(j) => {
                let q: LorentzVector<T> = self.propagators[j].q.cast();
                let b: T = (mom + q).spatial_squared()
                    + T::from_f64(self.propagators[j].m_squared).unwrap();
                mom[0] = b.sqrt();
                mom.x += q.x;
                mom.y += q.y;
                mom.z += q.z;
                mom
            }
            Cut::NegativeCut(j) => {
                let q: LorentzVector<T> = self.propagators[j].q.cast();
                let b: T = (mom + q).spatial_squared()
                    + T::from_f64(self.propagators[j].m_squared).unwrap();
                mom[0] = -b.sqrt();
                mom.x += q.x;
                mom.y += q.y;
                mom.z += q.z;
                mom
            }
            Cut::NoCut => unreachable!("Cannot compute the energy part of a non-cut loop line"),
        }
    }

    /// Get the energy of a cut with the shift subtracted
    // TODO: remove this function as soon as num::complex implements Float
    fn get_loop_energy(&self, loop_momenta: &[LorentzVector<Complex>], cut: &Cut) -> Complex {
        // construct the loop momentum that flows through this loop line
        let mut mom = LorentzVector::default();
        for (l, &c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom += l * Complex::new(float::from_i8(c).unwrap(), float::zero());
        }

        match *cut {
            Cut::PositiveCut(j) => {
                let q: LorentzVector<Complex> = self.propagators[j].q.cast();
                -q.t + ((mom + q).spatial_squared()
                    + float::from_f64(self.propagators[j].m_squared).unwrap())
                .sqrt()
            }
            Cut::NegativeCut(j) => {
                let q: LorentzVector<Complex> = self.propagators[j].q.cast();
                -q.t - ((mom + q).spatial_squared()
                    + float::from_f64(self.propagators[j].m_squared).unwrap())
                .sqrt()
            }
            Cut::NoCut => unreachable!("Cannot compute the energy part of a non-cut loop line"),
        }
    }

    /// Return the inverse of the evaluated loop line
    fn evaluate(
        &self,
        loop_momenta: &[LorentzVector<Complex>],
        cut: &Cut,
        threshold: float,
        kinematics_scale: float,
    ) -> Result<Complex, &'static str> {
        // construct the loop momentum that flows through this loop line
        let mut mom: LorentzVector<Complex> = LorentzVector::default();
        for (l, &c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom += l * Complex::new(float::from_i8(c).unwrap(), float::zero());
        }

        let mut res = Complex::new(float::one(), float::zero());
        for (i, p) in self.propagators.iter().enumerate() {
            match cut {
                Cut::PositiveCut(j) if i == *j => {
                    res *= (mom.t + float::from_f64(p.q.t).unwrap()) * float::from_f64(2.).unwrap();
                }
                Cut::NegativeCut(j) if i == *j => {
                    res *=
                        (mom.t + float::from_f64(p.q.t).unwrap()) * float::from_f64(-2.).unwrap();
                }
                _ => {
                    // multiply dual propagator
                    let pq: LorentzVector<Complex> = p.q.cast();
                    let m1: LorentzVector<Complex> = mom + pq;
                    let mut r = m1.square() - float::from_f64(p.m_squared).unwrap();

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
        for l in &mut self.loop_lines {
            for p in &mut l.propagators {
                p.signature = l.signature.clone();
            }
        }

        // determine the e_cm_squared and on-shell flag
        if self.external_kinematics.len() > 2 {
            self.e_cm_squared = (self.external_kinematics[0] + self.external_kinematics[1])
                .square_impr()
                .abs();
        } else {
            self.e_cm_squared = self.external_kinematics[0].square_impr().abs();
        }

        /*let mut v = LorentzVector::default();
        for e in &self.external_kinematics {
            if e.t < 0. {
                v = v + e;
            }
        }
        assert!(v.square() != 0.);
        self.e_cm_squared = v.square();
        println!("E_CM_SQ:{}", self.e_cm_squared);
        */

        // determine the on-shell flag
        for (i, e) in self.external_kinematics.iter().enumerate() {
            if e.square_impr() < 1e-10 * self.e_cm_squared {
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
                for (cut, ll) in cut_option.iter().zip(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut
                    {
                        cut_shift.push(ll.propagators[*cut_prop_index].q.cast());
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
                        // 2. surface_shift^2 > 0
                        // 3. all signs need to be the same (except 0)
                        // 4. surface_shift.t needs to have the opposite sign as in step 3.
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
                                    if surface_shift.square() > float::zero()
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
                                    // TODO: add existence condition for hyperboloids
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
                    "Found surface for {}: group={}, ellipsoid={}, prop={:?} cut={}, mom_map={:?}, signs={:?}, marker={}, shift={}",
                    self.name, s.group, s.ellipsoid, (s.onshell_ll_index, s.onshell_prop_index), CutList(&s.cut), s.sig_ll_in_cb,
                    s.signs, s.delta_sign, s.shift
                );
            }
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

                accum.push(ellipsoids_not_in_cut);
            }
            ellipsoids_not_in_cuts.push(accum);
        }
        self.ellipsoids_not_in_cuts = ellipsoids_not_in_cuts;
    }

    /// Map from unit hypercube to infinite hypercube in N-d
    pub fn parameterize(
        &self,
        x: &[f64],
        loop_index: usize,
    ) -> (ArrayVec<[float; MAX_DIM]>, float) {
        let mut jac = float::one();
        let e_cm = float::from_f64(self.e_cm_squared).unwrap().sqrt()
            * (float::from_f64(self.settings.parameterization.rescaling).unwrap())
                .powi(loop_index as i32);
        let radius =
            e_cm * float::from_f64(x[0]).unwrap() / (float::one() - float::from_f64(x[0]).unwrap()); // in [0,inf)
        jac *= (e_cm + radius).powi(2) / e_cm;
        let phi = float::from_f64(2.).unwrap()
            * <float as FloatConst>::PI()
            * float::from_f64(x[1]).unwrap();
        jac *= float::from_f64(2.).unwrap() * <float as FloatConst>::PI();

        match x.len() {
            2 => {
                let mut l_space = ArrayVec::new();
                l_space.push(radius * phi.cos());
                l_space.push(radius * phi.sin());
                jac *= radius;
                (l_space, jac)
            }
            3 => {
                let cos_theta =
                    -float::one() + float::from_f64(2.).unwrap() * float::from_f64(x[2]).unwrap(); // out of range
                jac *= float::from_f64(2.).unwrap();
                let sin_theta = (float::one() - cos_theta * cos_theta).sqrt();
                let mut l_space = ArrayVec::new();
                l_space.push(radius * sin_theta * phi.cos());
                l_space.push(radius * sin_theta * phi.sin());
                l_space.push(radius * cos_theta);

                // add a shift such that k=l is harder to be picked up by integrators such as cuhre
                l_space[0] += e_cm
                    * float::from_f64(self.settings.parameterization.shift.0 * loop_index as f64)
                        .unwrap();
                l_space[1] += e_cm
                    * float::from_f64(self.settings.parameterization.shift.1 * loop_index as f64)
                        .unwrap();
                l_space[2] += e_cm
                    * float::from_f64(self.settings.parameterization.shift.2 * loop_index as f64)
                        .unwrap();

                jac *= radius * radius; // spherical coord
                (l_space, jac)
            }
            x => unimplemented!("Unknown dimension {}", x),
        }
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

    fn determine_lambda<U: dual_num::Dim + dual_num::DimName>(
        &self,
        loop_momenta: &[LorentzVector<DualN<float, U>>],
        kappas: &[LorentzVector<DualN<float, U>>],
        lambda_max: f64,
    ) -> DualN<float, U>
    where
        dual_num::DefaultAllocator: dual_num::Allocator<float, U>,
        dual_num::Owned<float, U>: Copy,
    {
        // now determine the global lambda scaling by going through each propagator
        // for each cut and taking the minimum of their respective lambdas
        let mut lambda_sq = DualN::from_real(float::from_f64(lambda_max).unwrap().powi(2));

        let sigma =
            DualN::from_real(float::from_f64(self.settings.deformation.softmin_sigma).unwrap());
        let mut smooth_min_num = lambda_sq * (-lambda_sq / sigma).exp();
        let mut smooth_min_den = (-lambda_sq / sigma).exp();
        for (cuts, cb_to_lmb_mat) in self.ltd_cut_options.iter().zip(self.cb_to_lmb_mat.iter()) {
            for cut in cuts {
                // compute the real and imaginary part and the mass of the cut momentum
                let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
                let mut cut_shifts = [LorentzVector::default(); MAX_LOOP];
                let mut kappa_cuts = [LorentzVector::default(); MAX_LOOP];
                let mut mass_cuts = [0.; MAX_LOOP];
                let mut index = 0;
                for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
                    if *ll_cut != Cut::NoCut {
                        cut_momenta[index] = ll.get_cut_momentum(&loop_momenta, ll_cut);

                        // construct the complex part of the cut loop line momentum
                        // kappa is expressed in the loop momentum basis
                        let mut kappa_cut = LorentzVector::default();
                        for (kappa, &sign) in kappas.iter().zip(ll.signature.iter()) {
                            kappa_cut += kappa * DualN::from_real(float::from_i8(sign).unwrap());
                        }
                        kappa_cuts[index] = kappa_cut;

                        if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = ll_cut {
                            mass_cuts[index] = ll.propagators[*i].m_squared;
                            cut_shifts[index] = ll.propagators[*i].q;
                        }
                        index += 1;
                    }
                }

                // we have to make sure that our linear expansion for the deformation vectors is reasonable
                // for that we need lambda < c * (q_i^2^cut+m_i^2)/(kappa_i^cut * q_i^cut)
                for (mom_cut, shift_cut, kappa_cut, mass_cut) in izip!(
                    &cut_momenta[..self.n_loops],
                    &cut_shifts,
                    &kappa_cuts,
                    &mass_cuts
                ) {
                    let scf: LorentzVector<DualN<float, U>> = shift_cut.cast();
                    let m = mom_cut + scf;
                    let lambda_exp = DualN::from_real(
                        float::from_f64(self.settings.deformation.expansion_threshold).unwrap(),
                    ) * (m.spatial_squared_impr()
                        + DualN::from_real(float::from_f64(*mass_cut).unwrap()))
                        / kappa_cut.spatial_dot_impr(&m).abs();

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

                for (ll_cut, onshell_ll) in cut.iter().zip(self.loop_lines.iter()) {
                    // construct the complex part of the loop line momentum
                    let mut kappa_onshell = LorentzVector::default();
                    for (kappa, &c) in kappas.iter().zip(onshell_ll.signature.iter()) {
                        kappa_onshell += kappa * DualN::from_real(float::from_i8(c).unwrap());
                    }

                    // determine the map from cut momenta to loop line momentum
                    // FIXME: prevent allocating
                    let sig_ll_in_cb =
                        na::DMatrix::from_row_slice(self.n_loops, self.n_loops, cb_to_lmb_mat)
                            .transpose()
                            * na::DMatrix::from_row_slice(self.n_loops, 1, &onshell_ll.signature);

                    // multiply the signature by the sign of the cut
                    let mut onshell_signs = sig_ll_in_cb.clone();
                    for (ss, sign) in onshell_signs.iter_mut().zip(cut) {
                        if let Cut::NegativeCut(_) = sign {
                            *ss *= -1;
                        }
                    }

                    // solve for the energy of kappa cut momentum
                    // in a linear approximation, suitable for
                    // Weinzierl's lambda scaling
                    kappa_onshell.t = DualN::from_real(float::zero());
                    for (((s, kc), cm), &mass_sq) in onshell_signs
                        .iter()
                        .zip(kappa_cuts[..self.n_loops].iter())
                        .zip(cut_momenta.iter())
                        .zip(mass_cuts.iter())
                    {
                        kappa_onshell.t += kc.spatial_dot_impr(cm)
                            / (cm.spatial_squared_impr() + float::from_f64(mass_sq).unwrap())
                                .sqrt()
                            * float::from_i8(*s).unwrap();
                    }

                    let k0sq_inv = DualN::from_real(float::one()) / kappa_onshell.square_impr();

                    // if the kappa is 0, there is no need for rescaling
                    if !k0sq_inv.is_finite() {
                        continue;
                    }

                    // construct the real part of the loop line momentum
                    let mut mom: LorentzVector<DualN<float, U>> = LorentzVector::default();;
                    for (l, &c, cut_shift) in izip!(&cut_momenta, &sig_ll_in_cb, &cut_shifts) {
                        let shift: LorentzVector<DualN<float, U>> = cut_shift.cast();
                        mom += (l - shift) * DualN::from_real(float::from_i8(c).unwrap());
                    }

                    for (prop_index, onshell_prop) in onshell_ll.propagators.iter().enumerate() {
                        let pq: LorentzVector<DualN<float, U>> = onshell_prop.q.cast();
                        let onshell_prop_mom = mom + pq;

                        // if this is the cut propagator, we only need to use the spatial part
                        // the functional form remains the same
                        let (x, y) = if *ll_cut == Cut::NegativeCut(prop_index)
                            || *ll_cut == Cut::PositiveCut(prop_index)
                        {
                            let k_sp_inv = DualN::from_real(float::one())
                                / kappa_onshell.spatial_squared_impr();
                            let x = (kappa_onshell.spatial_dot_impr(&onshell_prop_mom) * k_sp_inv)
                                .powi(2);
                            let y = (onshell_prop_mom.spatial_squared_impr()
                                - float::from_f64(onshell_prop.m_squared).unwrap())
                                * k_sp_inv;
                            (x, y)
                        } else {
                            if self.settings.deformation.skip_hyperboloids {
                                // skip to check propagators that will dual cancel anyway if they are 0
                                continue;
                            }

                            let x = (kappa_onshell.dot_impr(&onshell_prop_mom) * k0sq_inv).powi(2);
                            let y = (onshell_prop_mom.square_impr()
                                - float::from_f64(onshell_prop.m_squared).unwrap())
                                * k0sq_inv;
                            (x, y)
                        };
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
        }

        if sigma.is_zero() {
            lambda_sq.sqrt()
        } else {
            (smooth_min_num / smooth_min_den).sqrt()
        }
    }

    /// Construct a deformation vector by going through all the ellipsoids
    fn deform_ellipsoids<U: dual_num::Dim + dual_num::DimName>(
        &self,
        loop_momenta: &[LorentzVector<DualN<float, U>>],
    ) -> [LorentzVector<DualN<float, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<float, U>,
        dual_num::Owned<float, U>: Copy,
    {
        let mut cut_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut deform_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut kappas = [LorentzVector::default(); MAX_LOOP];

        // TODO: store globally. How though, as it's type is DualN<float, U>...
        // plus it has self.n_loops components...
        let mut kappa_surf = [LorentzVector::default(); MAX_LOOP * MAX_ELLIPSE_GROUPS];
        let mut inv_surf_prop = [DualN::default(); MAX_ELLIPSE_GROUPS];

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
                    let mut mom: LorentzVector<DualN<float, U>> = LorentzVector::default();
                    for (l, &c) in loop_momenta.iter().zip(ll.signature.iter()) {
                        mom += l * DualN::from_real(float::from_i8(c).unwrap());
                    }
                    let q: LorentzVector<DualN<float, U>> = ll.propagators[*i].q.cast();
                    mom += q;
                    let length: DualN<float, U> = mom.spatial_squared_impr()
                        + float::from_f64(ll.propagators[*i].m_squared).unwrap();
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
                surface_dir += l * DualN::from_real(float::from_i8(sig).unwrap());
            }

            let surface_prop =
                &self.loop_lines[surf.onshell_ll_index].propagators[surf.onshell_prop_index];
            let q: LorentzVector<DualN<float, U>> = surface_prop.q.cast();
            surface_dir += q;

            let length: DualN<float, U> = surface_dir.spatial_squared_impr()
                + float::from_f64(surface_prop.m_squared).unwrap();
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
                    *cut_dir += surface_dir_norm * DualN::from_real(float::from_i8(sign).unwrap());
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
                    *deform_dir += cut_dir * DualN::from_real(float::from_i8(sign).unwrap());
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
                    cut_momenta[index] = ll.get_cut_momentum(loop_momenta, ll_cut);
                    index += 1;
                }
            }

            // evaluate the inverse propagator of the surface
            // the momentum map from the cut momenta is in signs and the shift is known as well
            // compute the energies for the loop momenta
            let mut mom: LorentzVector<DualN<float, U>> = LorentzVector::new();
            for (&s, c) in surf.sig_ll_in_cb.iter().zip(cut_momenta.iter()) {
                mom += c * DualN::from_real(float::from_i8(s).unwrap());
            }

            let qs: LorentzVector<DualN<float, U>> = surf.shift.cast();
            let momq: LorentzVector<DualN<float, U>> = mom + qs;
            let inv = momq.square_impr() - float::from_f64(surface_prop.m_squared).unwrap();
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
                let aij: DualN<float, U> =
                    NumCast::from(self.settings.deformation.additive.a_ij).unwrap();

                for (i, &inv) in inv_surf_prop[..group_counter].iter().enumerate() {
                    let dampening = match self.settings.deformation.additive.mode {
                        AdditiveMode::Exponential => (-inv * inv
                            / (aij * float::from_f64(self.e_cm_squared).unwrap().powi(2)))
                        .exp(),
                        AdditiveMode::Hyperbolic => {
                            let t = inv * inv / float::from_f64(self.e_cm_squared).unwrap().powi(2);
                            t / (t + aij)
                        }
                        AdditiveMode::Unity => DualN::from_real(float::one()),
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
                let m_ij: DualN<float, U> =
                    NumCast::from(self.settings.deformation.multiplicative.m_ij).unwrap();

                for (i, _) in inv_surf_prop[..group_counter].iter().enumerate() {
                    let mut dampening = DualN::from_real(float::one());
                    for (j, &inv) in inv_surf_prop[..group_counter].iter().enumerate() {
                        if i == j {
                            continue;
                        }

                        let t = inv * inv / float::from_f64(self.e_cm_squared).unwrap().powi(2);
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

    /// Construct a deformation vector by going through all cut options
    /// and make sure it is 0 on all other cut options
    fn deform_cutgroups<U: dual_num::Dim + dual_num::DimName>(
        &self,
        loop_momenta: &[LorentzVector<DualN<float, U>>],
    ) -> [LorentzVector<DualN<float, U>>; MAX_LOOP]
    where
        dual_num::DefaultAllocator: dual_num::Allocator<float, U>,
        dual_num::Owned<float, U>: Copy,
    {
        let mut deform_dirs = [LorentzVector::default(); MAX_ELLIPSE_GROUPS * MAX_LOOP];
        // TODO: precompute
        let mut non_empty_cuts = [(0, 0); MAX_ELLIPSE_GROUPS];

        let mut ellipsoid_eval = [DualN::default(); MAX_ELLIPSE_GROUPS * MAX_LOOP];

        // go through all cuts that contain ellipsoids
        // by going through all ellipsoids and skipping ones that are in the same cut
        let mut last_cut = None;
        let mut non_empty_cut_count = 0;
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            if !surf.ellipsoid {
                continue;
            }

            // evaluate ellipsoid representatives
            if surf.group == surf_index {
                // FIXME: we compute this again later
                let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
                let mut index = 0;

                // calculate the cut momenta
                let mut res = DualN::from_real(float::zero());
                for (&ll_cut, ll) in surf.cut.iter().zip(self.loop_lines.iter()) {
                    if ll_cut != Cut::NoCut {
                        cut_momenta[index] = ll.get_cut_momentum(loop_momenta, &ll_cut);
                        index += 1;
                    }
                }

                let mut q3 = LorentzVector::default();
                for (&s, mom) in surf.sig_ll_in_cb.iter().zip(&cut_momenta) {
                    if s != 0 {
                        q3 += mom * DualN::from_real(float::from_i8(s).unwrap());
                        res += mom.t.abs();
                    }
                }
                let shift: LorentzVector<DualN<float, U>> = surf.shift.cast();
                q3 += shift;

                res += q3.spatial_squared_impr().sqrt();
                res -= shift.t.abs();
                ellipsoid_eval[surf_index] = res;
            }

            if Some(&surf.cut) != last_cut {
                // compute each cut momentum, subtract the shift and normalize
                let mut index = 0;
                let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
                let mut cut_masses = [DualN::default(); MAX_LOOP];
                for (ll_cut, ll) in surf.cut.iter().zip(self.loop_lines.iter()) {
                    if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = *ll_cut {
                        cut_momenta[index] = ll.get_cut_momentum(loop_momenta, ll_cut);

                        // subtract the shift from q0
                        cut_momenta[index].t -=
                            DualN::from_real(float::from_f64(ll.propagators[i].q.t).unwrap());
                        cut_masses[index] =
                            DualN::from_real(float::from_f64(ll.propagators[i].m_squared).unwrap());

                        index += 1;
                    }
                }

                // convert from cut momentum basis to loop momentum basis
                for i in 0..self.n_loops {
                    let mut deform_dir = LorentzVector::default();
                    for ((&sign, cut_dir), &cut_mass) in self.cb_to_lmb_mat
                        [surf.cut_structure_index][i * self.n_loops..(i + 1) * self.n_loops]
                        .iter()
                        .zip(cut_momenta.iter())
                        .zip(cut_masses.iter())
                    {
                        // normalize the spatial components
                        let mut mom_normalized = *cut_dir;
                        let length = (mom_normalized.spatial_squared_impr() + cut_mass).sqrt();
                        mom_normalized *= length.inv();

                        deform_dir +=
                            mom_normalized * DualN::from_real(float::from_i8(sign).unwrap());
                    }

                    deform_dirs[non_empty_cut_count * MAX_LOOP + i] = deform_dir;
                }

                non_empty_cuts[non_empty_cut_count] =
                    (surf.cut_structure_index, surf.cut_option_index);

                last_cut = Some(&surf.cut);
                non_empty_cut_count += 1;
            }
        }

        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        for (i, &(cut_structure_index, cut_option_index)) in
            non_empty_cuts[..non_empty_cut_count].iter().enumerate()
        {
            let mut s = DualN::from_real(float::one());
            for &surf_index in &self.ellipsoids_not_in_cuts[cut_structure_index][cut_option_index] {
                // t is the weighing factor that is 0 if we are on both cut_i and cut_j
                // at the same time and goes to 1 otherwise
                let mut t = ellipsoid_eval[surf_index].powi(2);
                s *= t
                    / (t + DualN::from_real(
                        float::from_f64(
                            self.settings.deformation.cutgroups.m_ij * self.e_cm_squared,
                        )
                        .unwrap(),
                    ));
            }

            for ii in 0..self.n_loops {
                kappas[ii] -= deform_dirs[i * MAX_LOOP + ii] * s;
            }
        }
        kappas
    }

    fn deform_generic<U: dual_num::Dim + dual_num::DimName>(
        &self,
        loop_momenta: &[LorentzVector<DualN<float, U>>],
    ) -> ([LorentzVector<float>; MAX_LOOP], Complex)
    where
        dual_num::DefaultAllocator: dual_num::Allocator<float, U>,
        dual_num::Owned<float, U>: Copy,
    {
        let mut kappas = match self.settings.general.deformation_strategy {
            DeformationStrategy::CutGroups => self.deform_cutgroups(loop_momenta),
            _ => self.deform_ellipsoids(loop_momenta),
        };

        // make sure the kappa has the right dimension by multiplying in the scale
        let scale = DualN::from_real(float::from_f64(self.e_cm_squared.sqrt()).unwrap());
        for kappa in &mut kappas[..self.n_loops] {
            *kappa *= scale;
        }

        let lambda = if self.settings.deformation.lambda > 0. {
            self.determine_lambda(loop_momenta, &kappas, self.settings.deformation.lambda)
        } else {
            NumCast::from(self.settings.deformation.lambda.abs()).unwrap()
        };

        for k in kappas[..self.n_loops].iter_mut() {
            *k *= lambda;
            k.t = DualN::from_real(float::zero()); // make sure we do not have a left-over deformation
        }

        let mut jac_mat = ArrayVec::new();
        jac_mat.extend((0..9 * self.n_loops * self.n_loops).map(|_| Complex::default()));
        for i in 0..3 * self.n_loops {
            jac_mat[i * 3 * self.n_loops + i] += Complex::new(float::one(), float::zero());
            for j in 0..3 * self.n_loops {
                // first index: loop momentum, second: xyz, third: dual
                jac_mat[i * 3 * self.n_loops + j] +=
                    Complex::new(float::zero(), kappas[i / 3][i % 3 + 1][j + 1]);
            }
        }

        let jac = utils::determinant(&jac_mat);
        // take the real part
        let mut r = [LorentzVector::default(); MAX_LOOP];
        for (rr, k) in r.iter_mut().zip(&kappas[..self.n_loops]) {
            *rr = k.real();
        }
        (r, jac)
    }

    pub fn deform(
        &self,
        loop_momenta: &[LorentzVector<float>],
    ) -> ([LorentzVector<float>; MAX_LOOP], Complex) {
        if DeformationStrategy::None == self.settings.general.deformation_strategy {
            let r = [LorentzVector::default(); MAX_LOOP];
            return (r, Complex::new(float::one(), float::zero()));
        }

        match self.n_loops {
            1 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];

                r[0] = loop_momenta[0].map(|x| Dual4::from_real(x));
                for i in 0..3 {
                    r[0][i + 1][i + 1] = float::one();
                }

                return self.deform_generic(&r);
            }
            2 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual7::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual7::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = float::one();
                    r[1][i + 1][i + 4] = float::one();
                }
                self.deform_generic(&r)
            }
            3 => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                r[0] = loop_momenta[0].map(|x| Dual10::from_real(x));
                r[1] = loop_momenta[1].map(|x| Dual10::from_real(x));
                r[2] = loop_momenta[2].map(|x| Dual10::from_real(x));

                for i in 0..3 {
                    r[0][i + 1][i + 1] = float::one();
                    r[1][i + 1][i + 4] = float::one();
                    r[2][i + 1][i + 7] = float::one();
                }
                self.deform_generic(&r)
            }
            n => panic!("Binding for deformation at {} loops is not implemented", n),
        }
    }

    /// Set the energy component of the loop momenta according to
    /// `cut`.
    #[inline]
    pub fn set_loop_momentum_energies(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
    ) {
        // compute the cut energy for each loop line
        let mut cut_energy = [Complex::default(); MAX_LOOP];
        let mut index = 0;
        for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
            if *ll_cut != Cut::NoCut {
                cut_energy[index] = ll.get_loop_energy(&k_def, ll_cut);
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
                l.t += e * float::from_i8(*c).unwrap();
            }
        }
    }

    #[inline]
    pub fn evaluate_cut(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
    ) -> Result<Complex, &str> {
        self.set_loop_momentum_energies(k_def, cut, mat);

        let mut r = Complex::new(float::one(), float::zero());
        for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
            r *= ll.evaluate(
                &k_def,
                ll_cut,
                float::from_f64(self.e_cm_squared).unwrap(),
                float::from_f64(self.settings.general.numerical_threshold).unwrap(),
            )?;
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

    #[inline]
    pub fn evaluate<'a>(
        &mut self,
        x: &'a [f64],
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex>; MAX_LOOP]>,
        float,
        Complex,
        Complex,
    ) {
        // parameterize
        let mut k = [LorentzVector::default(); MAX_LOOP];
        let mut jac_para = float::one();
        for i in 0..self.n_loops {
            let (mut l_space, jac) = self.parameterize(&x[i * 3..(i + 1) * 3], i);

            let rot = self.rotation_matrix;
            k[i] = LorentzVector::from_args(
                float::zero(),
                rot[0][0] * l_space[0] + rot[0][1] * l_space[1] + rot[0][2] * l_space[2],
                rot[1][0] * l_space[0] + rot[1][1] * l_space[1] + rot[1][2] * l_space[2],
                rot[2][0] * l_space[0] + rot[2][1] * l_space[1] + rot[2][2] * l_space[2],
            );

            jac_para *= jac;
        }

        // deform
        let (kappas, jac_def) = self.deform(&k);
        let mut k_def: ArrayVec<[LorentzVector<Complex>; MAX_LOOP]> = (0..self.n_loops)
            .map(|i| {
                k[i].map(|x| Complex::new(x, float::zero()))
                    + kappas[i].map(|x| Complex::new(float::zero(), x))
            })
            .collect();

        let mut result = Complex::default();
        for (cuts, mat) in self.ltd_cut_options.iter().zip(self.cb_to_lmb_mat.iter()) {
            // for each cut coming from the same cut structure
            for cut in cuts {
                match self.evaluate_cut(&mut k_def, cut, mat) {
                    Ok(v) => result += v,
                    Err(_) => return (x, k_def, jac_para, jac_def, Complex::default()),
                }
            }
        }

        // add counterterm
        result += self.counterterm();

        result *= utils::powi(Complex::new(
            float::zero(),
            float::from_f64(-2.).unwrap() * <float as FloatConst>::PI(),
        ), self.n_loops); // factor of delta cut

        result *= utils::powi(Complex::new(
            float::from_f64(1.).unwrap()
                / (float::from_f64(2.).unwrap() * <float as FloatConst>::PI()).powi(4),
            float::zero(),
        )
        , self.n_loops); // loop momentum factor
        result *= jac_def * jac_para;

        (x, k_def, jac_para, jac_def, result)
    }
}
