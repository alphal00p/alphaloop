use arrayvec::ArrayVec;
use dual_num::{DualN, U4, U7};
use itertools::Itertools;
use na::Vector3;
use num_traits::Float;
use std::f64::consts::PI;
use topologies::{Cut, LoopLine, Topology};
use vector::{Field, LorentzVector, RealField};

use utils;
use utils::finv;

const MAX_PROP: usize = 8;
const MAX_DIM: usize = 3;
const MAX_LOOP: usize = 3;

type Dual4 = DualN<f64, U4>;
type Dual7 = DualN<f64, U7>;
type Complex = num::Complex<f64>;

impl LoopLine {
    /// Get the momenta of the cut with the cut evaluated
    /// `(+/-sqrt(|q_cut^2| + m^2), q_cut_vec)`
    /// where `q_cut` is the cut momentum (including affine term).
    fn get_cut_momentum<T: Float + RealField>(
        &self,
        loop_momenta: &[LorentzVector<T>],
        cut: &Cut,
    ) -> LorentzVector<T> {
        // construct the loop momentum that flows through this loop line
        let mut mom = LorentzVector::default();
        for (l, c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom = mom + l * From::from(*c as f64);
        }

        match *cut {
            Cut::PositiveCut(j) => {
                let q = self.propagators[j].q.map(|x| From::from(x));
                let b: T = (&mom + &q).spatial_squared() + self.propagators[j].m_squared;
                mom[0] = b.sqrt();
                mom.x += q.x;
                mom.y += q.y;
                mom.z += q.z;
                mom
            }
            Cut::NegativeCut(j) => {
                let q = self.propagators[j].q.map(|x| From::from(x));
                let b: T = (&mom + &q).spatial_squared() + self.propagators[j].m_squared;
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
        for (l, c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom = mom + l * From::from(*c as f64);
        }

        match *cut {
            Cut::PositiveCut(j) => {
                -self.propagators[j].q.t
                    + ((&mom + &self.propagators[j].q).spatial_squared()
                        + self.propagators[j].m_squared)
                        .sqrt()
            }
            Cut::NegativeCut(j) => {
                -self.propagators[j].q.t
                    - ((&mom + &self.propagators[j].q).spatial_squared()
                        + self.propagators[j].m_squared)
                        .sqrt()
            }
            Cut::NoCut => unreachable!("Cannot compute the energy part of a non-cut loop line"),
        }
    }

    fn evaluate(&self, loop_momenta: &[LorentzVector<Complex>], cut: &Cut) -> Complex {
        // construct the loop momentum that flows through this loop line
        let mut mom = LorentzVector::default();
        for (l, c) in loop_momenta.iter().zip(self.signature.iter()) {
            mom = mom + l * Complex::new(*c as f64, 0.);
        }

        let mut res = Complex::new(1., 0.);
        for (i, p) in self.propagators.iter().enumerate() {
            match cut {
                Cut::PositiveCut(j) if i == *j => {
                    res *= 2. * (mom.t + p.q.t);
                }
                Cut::NegativeCut(j) if i == *j => {
                    res *= -2. * (mom.t + p.q.t);
                }
                _ => {
                    // multiply dual propagator
                    res *= (&mom + &p.q).square() - p.m_squared;
                }
            }
        }
        finv(res)
    }
}

impl Topology {
    /// Compute LTD related quantities for this topology.
    pub fn process(&mut self) {
        if self.n_loops > MAX_LOOP {
            panic!("Please increase MAX_LOOP to {}", self.n_loops);
        }

        // copy the signature to the propagators
        for l in &mut self.loop_lines {
            for p in &mut l.propagators {
                p.signature = l.signature.clone();
            }
        }

        // determine the e_cm_squared and on-shell flag
        if self.external_kinematics.len() > 2 {
            self.e_cm_squared = (self.external_kinematics[0] + self.external_kinematics[1])
                .square()
                .abs();
        } else {
            self.e_cm_squared = self.external_kinematics[0].square().abs();
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
            if e.square() < 1e-10 * self.e_cm_squared {
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
            for cut_option in cut_options.iter() {
                let mut cut_shift = vec![]; // qs of cut
                for (cut, ll) in cut_option.iter().zip(self.loop_lines.iter()) {
                    if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut
                    {
                        cut_shift.push(ll.propagators[*cut_prop_index].q);
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

                        let mut surface_shift = onshell_prop.q;
                        for (&sign, q) in sig_ll_in_cb.iter().zip(cut_shift.iter()) {
                            surface_shift = surface_shift - q * sign as f64;
                        }

                        // now see if we have an ellipsoid
                        // 1. surface_shift != 0
                        // 2. surface_shift^2 > 0
                        // 3. all signs need to be the same (except 0)
                        // 4. surface_shift.t needs to have the opposite sign as in step 3.
                        if surface_shift.t != 0. {
                            // multiply the residue sign
                            let mut surface_signs = sig_ll_in_cb.clone();
                            for (ss, sign) in surface_signs.iter_mut().zip(residue_sign) {
                                if *sign < 0 {
                                    *ss *= -1;
                                }
                            }

                            let surface_sign_sum = surface_signs.iter().sum::<i8>();
                            let surface_signs_abs_sum =
                                surface_signs.iter().map(|x| x.abs()).sum::<i8>();
                            // check the two branches Delta+ and Delta-
                            for &delta_sign in &[-1, 1] {
                                // make sure all non-zero coefficients are the same
                                if (surface_sign_sum + delta_sign).abs()
                                    == surface_signs_abs_sum + delta_sign.abs()
                                {
                                    if surface_shift.square() > 0.
                                        && delta_sign as f64 * surface_shift.t < 0.
                                    {
                                        /*
                                        println!(
                                            "Found ellipsoid for {}: cut={:?}, mom_map={:?}, signs={:?}, marker={}, shift={}",
                                            self.name, cut_option, sig_ll_in_cb.iter().cloned().collect::<Vec<_>>(),
                                            surface_signs.iter().cloned().collect::<Vec<_>>(), delta_sign, surface_shift
                                        );
                                        */
                                        self.ellipsoids.push((
                                            cut_index,
                                            cut_option.clone(),
                                            ll_index,
                                            onshell_prop_index,
                                            delta_sign,
                                            sig_ll_in_cb.iter().cloned().collect(),
                                            surface_shift.clone(),
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Map from unit hypercube to infinite hypercube in N-d
    pub fn parameterize(&self, x: &[f64]) -> (ArrayVec<[f64; MAX_DIM]>, f64) {
        let mut jac = 1.;
        let e_cm = self.e_cm_squared.sqrt();
        let radius = e_cm * x[0] / (1. - x[0]); // in [0,inf)
        jac *= (e_cm + radius).powi(2) / e_cm;
        assert!(x.len() > 1);
        let phi = 2. * PI * x[1];
        jac *= 2. * PI;

        match x.len() {
            2 => {
                let mut l_space = ArrayVec::new();
                l_space.push(radius * phi.cos());
                l_space.push(radius * phi.sin());
                jac *= radius;
                (l_space, jac)
            }
            3 => {
                let cos_theta = -1. + 2. * x[2];
                jac *= 2.;
                let sin_theta = (1. - cos_theta * cos_theta).sqrt();
                let mut l_space = ArrayVec::new();
                l_space.push(radius * sin_theta * phi.cos());
                l_space.push(radius * sin_theta * phi.sin());
                l_space.push(radius * cos_theta);
                jac *= radius * radius; // spherical coord
                (l_space, jac)
            }
            x => unimplemented!("Unknown dimension {}", x),
        }
    }

    #[inline]
    fn compute_lambda_factor<T: Field + RealField>(x: T, y: T) -> T {
        if x * 2. < y {
            y * 0.25
        } else if y < 0. {
            x - y * 0.5
        } else {
            x - y * 0.25
        }
    }

    /// Compute the deformation vector and the deformation Jacobian wrt the loop momentum
    /// in the loop line
    pub fn deform_one_loop(
        &self,
        loop_momenta: &[LorentzVector<f64>],
    ) -> ([LorentzVector<f64>; MAX_LOOP], Complex) {
        match self.deformation.as_ref() {
            "none" => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                return (r, Complex::new(1., 0.));
            }
            "generic" => {
                let mut r = [LorentzVector::default(); MAX_LOOP];

                r[0] = loop_momenta[0].map(|x| Dual4::from_real(x));
                for i in 0..3 {
                    r[0][i + 1][i + 1] = 1.0;
                }

                return self.deform_generic(&r);
            }
            _ => {}
        }

        let ll = &self.loop_lines[0];

        let mut mom = Vector3::zeros();
        for (l, c) in loop_momenta.iter().zip(ll.signature.iter()) {
            mom += Vector3::new(l.x, l.y, l.z) * (*c as f64);
        }

        let mut dual_x = mom.map(|y| Dual4::from_real(y));
        for (i, x) in dual_x.iter_mut().enumerate() {
            x[i + 1] = 1.0;
        }

        let ki_spatial: ArrayVec<[_; MAX_PROP]> = ll
            .propagators
            .iter()
            .map(|p| p.q.map(|y| Dual4::from_real(y)))
            .map(|x| Vector3::new(x.x, x.y, x.z) + dual_x)
            .collect();
        let ki_plus: ArrayVec<[_; MAX_PROP]> = ki_spatial
            .iter()
            .zip(&ll.propagators)
            .map(|(x, p)| (x.dot(x) + p.m_squared).sqrt())
            .collect();

        let dirs: ArrayVec<[_; MAX_PROP]> =
            ki_spatial.iter().map(|x| x / x.dot(x).sqrt()).collect();
        let aij = 0.01;

        // TODO: add lambda per kappa_i?
        // TODO: split integration into groups
        let mut kappa = Vector3::zeros();
        for (i, (pi, dplus)) in ll.propagators.iter().zip(ki_plus.iter()).enumerate() {
            for (j, (pj, pplus)) in ll.propagators.iter().zip(ki_plus.iter()).enumerate() {
                if i != j {
                    // only deform if this combination is an ellipsoid singularity
                    //if self.singularity_matrix[i * self.q_and_mass.len() + j] {
                    let inv = *dplus + *pplus + Dual4::from_real(pj.q[0] - pi.q[0]);
                    let e = (-inv * inv / (aij * self.e_cm_squared.powi(2))).exp();
                    kappa += -(dirs[i] + dirs[j]) * e;
                    //}
                }
            }
        }

        // now determine the global lambda scaling
        let mut lambda_sq = Dual4::from_real(1.);
        for (i, ((p_cut, ki_s_cut), &ki_plus_cut)) in ll
            .propagators
            .iter()
            .zip(ki_spatial.iter())
            .zip(ki_plus.iter())
            .enumerate()
        {
            for (j, ((p_dual, ki_s_dual), &ki_plus_dual)) in ll
                .propagators
                .iter()
                .zip(ki_spatial.iter())
                .zip(ki_plus.iter())
                .enumerate()
            {
                let (a, b, c) = if i == j {
                    // "propagator" coming from delta function cut
                    let a = -kappa.dot(&kappa);
                    let b = ki_s_cut.dot(&kappa) * 2.;
                    let c = ki_plus_cut * ki_plus_cut;
                    (a, b, c)
                } else {
                    let e = (p_dual.q[0] - p_cut.q[0]).powi(2);
                    let c1 = ki_plus_cut * ki_plus_cut + e - ki_plus_dual * ki_plus_dual;
                    let b1 = (ki_s_cut - ki_s_dual).dot(&kappa) * 2.0;
                    let a = -b1 * b1 + kappa.dot(&kappa) * e * 4.;
                    let b = b1 * c1 * 2. - ki_s_cut.dot(&kappa) * e * 8.;
                    let c = c1 * c1 - ki_plus_cut * ki_plus_cut * e * 4.;
                    (a, b, c)
                };

                let x = b * b / (a * a * 4.);
                let y = -c / a;
                let ls = Topology::compute_lambda_factor(x, y);

                if ls < lambda_sq {
                    lambda_sq = ls;
                }
            }
        }

        kappa *= lambda_sq.sqrt();

        let mut jac_mat = [[Complex::new(0., 0.); 3]; 3];
        for (i, k) in kappa.iter().enumerate() {
            jac_mat[i][i] = Complex::new(1., 0.);
            for (j, kd) in k.iter().skip(1).enumerate() {
                jac_mat[i][j] += Complex::new(0., *kd);
            }
        }

        let jac = utils::determinant3x3(&jac_mat);

        let mut r = [LorentzVector::default(); MAX_LOOP];
        let r1 = kappa.map(|x| x.real());
        r[0] = LorentzVector::from_args(0., r1.x, r1.y, r1.z);
        (r, jac)
    }

    fn deform_two_loops(
        &self,
        loop_momenta: &[LorentzVector<f64>],
    ) -> ([LorentzVector<f64>; MAX_LOOP], Complex) {
        let mut r = [LorentzVector::default(); MAX_LOOP];

        match self.deformation.as_ref() {
            "none" => {
                // Below corresponds to no deformation
                (r, Complex::new(1., 0.))
            }
            "generic" => {
                let mut rr = [LorentzVector::default(); MAX_LOOP];
                rr[0] = loop_momenta[0].map(|x| Dual7::from_real(x));
                rr[1] = loop_momenta[1].map(|x| Dual7::from_real(x));

                for i in 0..3 {
                    rr[0][i + 1][i + 1] = 1.0;
                    rr[1][i + 1][i + 4] = 1.0;
                }

                return self.deform_generic(&rr);
            }
            "hardcoded_double_triangle_deformation" => {
                let k1 = Vector3::new(loop_momenta[0].x, loop_momenta[0].y, loop_momenta[0].z);
                let k2 = Vector3::new(loop_momenta[1].x, loop_momenta[1].y, loop_momenta[1].z);

                let mut k1_dual7 = k1.map(|x| Dual7::from_real(x));
                let mut k2_dual7 = k2.map(|x| Dual7::from_real(x));

                // Initialise the correct dual components of the input to 1.0
                for i in 0..6 {
                    if i < 3 {
                        k1_dual7[i][i + 1] = 1.0;
                    } else {
                        k2_dual7[i - 3][i + 1] = 1.0;
                    }
                }

                // hard-coded for now. The deformation should be a struct on its own
                // and *not* part of the loop line struct!
                let mut p_dual7 = Vector3::new(0., 0., 0.).map(|x| Dual7::from_real(x));
                let mut p0_dual7 = Dual7::from_real(1.0);
                let lambda1 = Dual7::from_real(1.0);
                let lambda2 = Dual7::from_real(1.0);
                // One should typically choose a dimensionful width sigma for the gaussian since the
                // dimensionality of the exponent coded below is GeV^2.

                // For now, hard-code it to 20.
                let sigma = Dual7::from_real(20.0);
                let sigma_squared = sigma * sigma;

                let exp_257 = k1_dual7.dot(&k1_dual7).sqrt()
                    + (k1_dual7 + p_dual7).dot(&(k1_dual7 + p_dual7)).sqrt()
                    - p0_dual7;
                let f257 = (-((exp_257 * exp_257) / sigma_squared)).exp();
                let exp_136 = k2_dual7.dot(&k2_dual7).sqrt()
                    + (k2_dual7 - p_dual7).dot(&(k2_dual7 - p_dual7)).sqrt()
                    - p0_dual7;
                let f136 = (-((exp_136 * exp_136) / sigma_squared)).exp();
                let exp_f4 = (k1_dual7 + p_dual7).dot(&(k1_dual7 + p_dual7)).sqrt()
                    + k2_dual7.dot(&k2_dual7).sqrt()
                    + (k1_dual7 + k2_dual7).dot(&(k1_dual7 + k2_dual7)).sqrt()
                    - p0_dual7;
                let f4 = (-((exp_f4 * exp_f4) / sigma_squared)).exp();

                let kappa1 = ((k1_dual7 / k1_dual7.dot(&k1_dual7).sqrt()
                    + (k1_dual7 + p_dual7)
                        / (k1_dual7 + p_dual7).dot(&(k1_dual7 + p_dual7)).sqrt())
                    * f257
                    + ((k1_dual7 + p_dual7)
                        / (k1_dual7 + p_dual7).dot(&(k1_dual7 + p_dual7)).sqrt()
                        + (k1_dual7 + k2_dual7)
                            / (k1_dual7 + k2_dual7).dot(&(k1_dual7 + k2_dual7)).sqrt())
                        * f4)
                    * lambda1;

                let kappa2 = ((k2_dual7 / k2_dual7.dot(&k2_dual7).sqrt()
                    + (k2_dual7 - p_dual7)
                        / (k2_dual7 - p_dual7).dot(&(k2_dual7 - p_dual7)).sqrt())
                    * f136
                    + (k2_dual7 / k2_dual7.dot(&k2_dual7).sqrt()
                        + (k1_dual7 + k2_dual7)
                            / (k1_dual7 + k2_dual7).dot(&(k1_dual7 + k2_dual7)).sqrt())
                        * f4)
                    * lambda2;

                let mut jac_mat = [[Complex::new(0., 0.); 6]; 6];
                for i in 0..6 {
                    jac_mat[i][i] += Complex::new(1., 0.);
                    for j in 0..6 {
                        if i < 3 {
                            jac_mat[i][j] += Complex::new(0.0, k1_dual7[i][j + 1]);
                        } else {
                            jac_mat[i][j] += Complex::new(0.0, k2_dual7[i - 3][j + 1]);
                        }
                    }
                }

                let kappa1_vec3 = kappa1.map(|x| x.real());
                let kappa2_vec3 = kappa2.map(|x| x.real());

                let jac = utils::determinant6x6(&jac_mat);

                r[0] = LorentzVector::from_args(0., kappa1_vec3[0], kappa1_vec3[1], kappa1_vec3[2]);
                r[1] = LorentzVector::from_args(0., kappa2_vec3[0], kappa2_vec3[1], kappa2_vec3[2]);

                (r, jac)
            }

            x => panic!("Requested deformation type unknown: {:?}", x),
        }
    }

    fn deform_generic<U: dual_num::Dim + dual_num::DimName>(
        &self,
        loop_momenta: &[LorentzVector<DualN<f64, U>>],
    ) -> ([LorentzVector<f64>; MAX_LOOP], Complex)
    where
        dual_num::DefaultAllocator: dual_num::Allocator<f64, U>,
        dual_num::Owned<f64, U>: Copy,
    {
        let mut cut_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut deform_dirs = [LorentzVector::default(); MAX_LOOP];
        let mut kappas = [LorentzVector::default(); MAX_LOOP];

        // surface equation:
        // m*|sum_i^L a_i q_i^cut + p_vec| + sum_i^L a_i*r_i * |q_i^cut| - p^0 = 0
        // where m=delta_sign a = sig_ll_in_cb, r = residue_sign
        for (
            cut_structure_index,
            cut,
            onshell_ll_index,
            onshell_prop_index,
            _delta_sign,
            sig_ll_in_cb,
            surface_shift,
        ) in self.ellipsoids.iter()
        {
            // compute v_i = q_i^cut / sqrt(|q_i^cut|^2 + m^2)
            // construct the normalized 3-momenta that flow through the cut propagators
            // the 0th component is to be ignored
            let mut cut_counter = 0;
            for (c, ll) in cut.iter().zip(self.loop_lines.iter()) {
                if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = c {
                    let mut mom = LorentzVector::default();
                    for (l, c) in loop_momenta.iter().zip(ll.signature.iter()) {
                        mom = mom + l * From::from(*c as f64);
                    }
                    mom = mom + ll.propagators[*i].q.map(|x| From::from(x));
                    mom = mom / (mom.spatial_squared() + ll.propagators[*i].m_squared).sqrt();
                    cut_dirs[cut_counter] = mom;
                    cut_counter += 1;
                }
            }

            // compute w=(sum_i^L a_i q_i^cut + p_vec) / sqrt(|sum_i^L a_i q_i^cut + p_vec|^2 + m^2)
            // determine sum_i^L a_i q_i^cut using the loop-momentum basis
            let mut surface_dir = LorentzVector::default();
            for (l, &sig) in loop_momenta
                .iter()
                .zip(self.loop_lines[*onshell_ll_index].signature.iter())
            {
                surface_dir = surface_dir + l * From::from(sig as f64);
            }

            let surface_prop = &self.loop_lines[*onshell_ll_index].propagators[*onshell_prop_index];
            surface_dir = surface_dir + surface_prop.q.map(|x| From::from(x));

            let surface_dir_norm =
                surface_dir / (surface_dir.spatial_squared() + surface_prop.m_squared).sqrt();

            // compute n_i = (v_i + a_i * w) * abs(a_i)
            // n_i is the normal vector of the ith entry of the cut momentum basis
            // this n_i is 0 when the surface does not depend on the ith cut loop line
            // we update cut_dirs from v_i to n_i
            for (cut_dir, &sign) in cut_dirs[..self.n_loops].iter_mut().zip(sig_ll_in_cb.iter()) {
                *cut_dir = (*cut_dir + surface_dir_norm * From::from(sign as f64))
                    * From::from(sign.abs() as f64);
            }

            // convert from cut momentum basis to loop momentum basis
            for (i, deform_dir) in deform_dirs[..self.n_loops].iter_mut().enumerate() {
                *deform_dir = LorentzVector::default();
                for (&sign, cut_dir) in self.cb_to_lmb_mat[*cut_structure_index]
                    [i * self.n_loops..(i + 1) * self.n_loops]
                    .iter()
                    .zip(cut_dirs[..self.n_loops].iter())
                {
                    *deform_dir = *deform_dir + cut_dir * From::from(sign as f64);
                }
            }

            // multiply the direction deform_dirs by the exponential function
            // for the exponential function we need to evaluate the propagator of the ellipsoid
            // and for that we need to determine the cut energies

            // compute the cut momenta with the affine term subtracted
            let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
            let mut index = 0;
            for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
                if *ll_cut != Cut::NoCut {
                    cut_momenta[index] = ll.get_cut_momentum(loop_momenta, ll_cut);
                    index += 1;
                }
            }

            // evaluate the inverse propagator of the surface
            // the momentum map from the cut momenta is in signs and the shift is known as well
            // compute the energies for the loop momenta
            let aij = 10.;
            let mut mom = LorentzVector::new();
            for (&s, c) in sig_ll_in_cb.iter().zip(cut_momenta.iter()) {
                mom = mom + c * From::from(s as f64);
            }

            let inv =
                (mom + surface_shift.map(|x| From::from(x))).square() - surface_prop.m_squared;
            let dampening = (-inv * inv / (aij * self.e_cm_squared.powi(2))).exp();

            for (kappa, dir) in kappas[..self.n_loops].iter_mut().zip(deform_dirs.iter()) {
                *kappa = *kappa + dir * From::from(dampening);
            }
        }

        // now determine the global lambda scaling by going through each propagator
        // for each cut and taking the minimum of their respective lambdas
        let mut lambda_sq = DualN::from_real(1.);
        for (cuts, cb_to_lmb_mat) in self.ltd_cut_options.iter().zip(self.cb_to_lmb_mat.iter()) {
            for cut in cuts {
                // compute the real and imaginary part and the mass of the cut momentum
                let mut cut_momenta = [LorentzVector::default(); MAX_LOOP];
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
                            kappa_cut = kappa_cut + kappa * From::from(sign as f64);
                        }
                        kappa_cuts[index] = kappa_cut;

                        if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = ll_cut {
                            mass_cuts[index] = ll.propagators[*i].m_squared;
                        }
                        index += 1;
                    }
                }

                for onshell_ll in self.loop_lines.iter() {
                    // construct the complex part of the loop line momentum
                    let mut kappa_onshell = LorentzVector::default();
                    for (kappa, &c) in kappas.iter().zip(onshell_ll.signature.iter()) {
                        kappa_onshell = kappa_onshell + kappa * From::from(c as f64);
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
                    kappa_onshell.t = DualN::from_real(0.);
                    for (((s, kc), cm), mass) in onshell_signs
                        .iter()
                        .zip(kappa_cuts.iter())
                        .zip(cut_momenta.iter())
                        .zip(mass_cuts.iter())
                    {
                        kappa_onshell.t +=
                            kc.spatial_dot(cm) / (cm.spatial_squared() + *mass) * 0.5 * (*s as f64);
                    }

                    let k0sq_inv = DualN::from_real(1.0) / kappa_onshell.square();

                    // if the kappa is 0, there is no need for rescaling
                    if !k0sq_inv.is_finite() {
                        continue;
                    }

                    // construct the real part of the loop line momentum
                    let mut mom = LorentzVector::default();
                    for (l, c) in cut_momenta.iter().zip(sig_ll_in_cb.iter()) {
                        mom = mom + l * From::from(*c as f64);
                    }

                    for onshell_prop in onshell_ll.propagators.iter() {
                        let onshell_prop_mom = mom + onshell_prop.q.map(|x| From::from(x));

                        let x = (kappa_onshell.dot(&onshell_prop_mom) * k0sq_inv).powi(2);
                        let y = (onshell_prop_mom.square() - onshell_prop.m_squared) * k0sq_inv;
                        let prop_lambda_sq = Topology::compute_lambda_factor(x, y);

                        if prop_lambda_sq < lambda_sq {
                            lambda_sq = prop_lambda_sq;
                        }
                    }
                }
            }
        }

        // TODO: add a lambda for the expansion parameter
        let lambda = lambda_sq.sqrt();
        for k in kappas[..self.n_loops].iter_mut() {
            *k = *k * lambda;
            k.t = DualN::from_real(0.); // make sure we do not have a left-over deformation
        }

        let mut jac_mat = ArrayVec::new();
        jac_mat.extend((0..9 * self.n_loops * self.n_loops).map(|_| Complex::default()));
        for i in 0..3 * self.n_loops {
            jac_mat[i * 3 * self.n_loops + i] += Complex::new(1., 0.);
            for j in 0..3 * self.n_loops {
                // first index: loop momentum, second: xyz, third: dual
                jac_mat[i * 3 * self.n_loops + j] +=
                    Complex::new(0.0, kappas[i / 3][i % 3 + 1][j + 1]);
            }
        }

        let jac = utils::determinant(&jac_mat);
        // take the real part
        let mut r = [LorentzVector::default(); MAX_LOOP];
        for (rr, k) in r.iter_mut().zip(&kappas[..self.n_loops]) {
            *rr = k.map(|x| x.real());
        }
        (r, jac)
    }

    pub fn deform(&self, k: &[LorentzVector<f64>]) -> ([LorentzVector<f64>; MAX_LOOP], Complex) {
        match self.n_loops {
            1 => self.deform_one_loop(k),
            2 => self.deform_two_loops(k),
            _ => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                (r, Complex::new(1., 0.))
            }
        }
    }

    #[inline]
    pub fn evaluate_cut(
        &self,
        k_def: &mut ArrayVec<[LorentzVector<Complex>; MAX_LOOP]>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
    ) -> Complex {
        // compute the cut energy for each loop line
        let mut cut_energy = [Complex::new(0., 0.); MAX_LOOP];
        let mut index = 0;
        for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
            if *ll_cut != Cut::NoCut {
                cut_energy[index] = ll.get_loop_energy(&k_def, ll_cut);
                index += 1;
            }
        }

        // compute the energies for the loop momenta
        for (i, l) in k_def.iter_mut().enumerate() {
            l.t = Complex::new(0., 0.);
            for (c, e) in mat[i * self.n_loops..(i + 1) * self.n_loops]
                .iter()
                .zip(&cut_energy[..self.n_loops])
            {
                l.t += e * *c as f64;
            }
        }

        let mut r = Complex::new(1., 0.);
        for (ll_cut, ll) in cut.iter().zip(self.loop_lines.iter()) {
            r *= ll.evaluate(&k_def, ll_cut);
        }
        r
    }

    #[inline]
    pub fn evaluate(&self, x: &[f64]) -> Complex {
        // parameterize
        let mut k = [LorentzVector::default(); MAX_LOOP];
        let mut jac_para = 1.;
        for i in 0..self.n_loops {
            let (l_space, jac) = self.parameterize(&x[i * 3..(i + 1) * 3]);
            k[i] = LorentzVector::from_args(0., l_space[0], l_space[1], l_space[2]);
            jac_para *= jac;
        }

        // deform
        let (kappas, jac_def) = self.deform(&k);
        let mut k_def: ArrayVec<[LorentzVector<Complex>; MAX_LOOP]> = (0..self.n_loops)
            .map(|i| k[i].map(|x| Complex::new(x, 0.)) + kappas[i].map(|x| Complex::new(0., x)))
            .collect();

        let mut result = Complex::new(0., 0.);
        for (cuts, mat) in self.ltd_cut_options.iter().zip(self.cb_to_lmb_mat.iter()) {
            // for each cut coming from the same cut structure
            for cut in cuts {
                result += self.evaluate_cut(&mut k_def, cut, mat);
            }
        }

        // add counterterm
        result += self.counterterm();

        result *= Complex::new(0., -2. * PI).powf(self.n_loops as f64); // factor of delta cut
        result *= Complex::new(0., -1. / (2. * PI).powi(4)).powf(self.n_loops as f64); // loop momentum factor
        result *= jac_para * jac_def;

        if !result.is_finite() {
            println!(
                "Bad result: {} with x={:?}, k_def={:?}, jac_para={}, jac_def={}",
                result, x, k_def, jac_para, jac_def
            );
        }

        result
    }
}
