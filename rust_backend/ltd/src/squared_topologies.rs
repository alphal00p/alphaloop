use arrayvec::ArrayVec;
use f128::f128;
use float;
use integrand::IntegrandImplementation;
use itertools::Itertools;
use num::Complex;
use num_traits::FromPrimitive;
use num_traits::NumCast;
use num_traits::ToPrimitive;
use num_traits::{Float, FloatConst};
use num_traits::{One, Zero};
use observables::EventManager;
use serde::Deserialize;
use std::fs::File;
use std::mem;
use topologies::{LTDCache, LTDNumerator, Topology};
use utils;
use utils::Signum;
use vector::{LorentzVector, RealNumberLike};
use {DeformationStrategy, FloatLike, Settings, MAX_LOOP};

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCut {
    pub name: String,
    pub sign: i8,
    pub level: usize,
    pub signature: (Vec<i8>, Vec<i8>),
    pub m_squared: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCutLimits {
    pub diagrams: Vec<Topology>,
    pub conjugate_deformation: Vec<bool>,
    pub symmetry_factor: f64,
    #[serde(default)]
    pub numerator_tensor_coefficients: Vec<(f64, f64)>,
    #[serde(skip_deserializing)]
    pub numerator: LTDNumerator,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCuts {
    pub cuts: Vec<CutkoskyCut>,
    pub n_bubbles: usize,
    pub uv_limits: Vec<CutkoskyCutLimits>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopology {
    pub name: String,
    pub n_loops: usize,
    pub n_incoming_momenta: usize,
    pub e_cm_squared: f64,
    pub overall_numerator: f64,
    pub external_momenta: Vec<LorentzVector<f64>>,
    pub cutkosky_cuts: Vec<CutkoskyCuts>,
    #[serde(skip_deserializing)]
    pub settings: Settings,
    #[serde(skip_deserializing)]
    pub rotation_matrix: [[float; 3]; 3],
}

#[derive(Default)]
pub struct SquaredTopologyCache {
    float_cache: Vec<Vec<Vec<LTDCache<float>>>>,
    quad_cache: Vec<Vec<Vec<LTDCache<f128>>>>,
}

pub trait CachePrecisionSelector<T: FloatLike> {
    fn get(&mut self) -> &mut Vec<Vec<Vec<LTDCache<T>>>>;
}

impl CachePrecisionSelector<float> for SquaredTopologyCache {
    #[inline]
    fn get(&mut self) -> &mut Vec<Vec<Vec<LTDCache<float>>>> {
        &mut self.float_cache
    }
}

impl CachePrecisionSelector<f128> for SquaredTopologyCache {
    #[inline]
    fn get(&mut self) -> &mut Vec<Vec<Vec<LTDCache<f128>>>> {
        &mut self.quad_cache
    }
}

impl SquaredTopology {
    pub fn from_file(filename: &str, settings: &Settings) -> SquaredTopology {
        let f = File::open(filename).expect("Could not open squared topology file");

        let mut squared_topo: SquaredTopology = serde_yaml::from_reader(f).unwrap();
        squared_topo.settings = settings.clone();
        for cutkosky_cuts in &mut squared_topo.cutkosky_cuts {
            for uv_limit in &mut cutkosky_cuts.uv_limits {
                uv_limit.numerator = if uv_limit.numerator_tensor_coefficients.len() == 0 {
                    LTDNumerator::one(squared_topo.n_loops + cutkosky_cuts.n_bubbles)
                } else {
                    LTDNumerator::new(
                        squared_topo.n_loops + cutkosky_cuts.n_bubbles,
                        &uv_limit
                            .numerator_tensor_coefficients
                            .iter()
                            .map(|x| Complex::new(x.0, x.1))
                            .collect::<Vec<_>>(),
                    )
                };

                for d in &mut uv_limit.diagrams {
                    d.settings = settings.clone();
                    d.process(false);
                }
            }
        }

        debug_assert_eq!(
            squared_topo.external_momenta.len(),
            squared_topo.n_incoming_momenta * 2,
            "The number of external momenta is wrong."
        );

        squared_topo.rotation_matrix = [
            [float::one(), float::zero(), float::zero()],
            [float::zero(), float::one(), float::zero()],
            [float::zero(), float::zero(), float::one()],
        ];

        squared_topo
    }

    fn evaluate_signature<T: RealNumberLike + FromPrimitive>(
        signature: &(Vec<i8>, Vec<i8>),
        external_momenta: &[LorentzVector<T>],
        loop_momenta: &[LorentzVector<T>],
    ) -> LorentzVector<T> {
        let mut cut_momentum = LorentzVector::default();
        for (&sign, mom) in signature.0.iter().zip_eq(loop_momenta) {
            if sign != 0 {
                // note: we allow for the sign to be any small integer
                cut_momentum += mom * T::from_i8(sign).unwrap();
            }
        }

        for (&sign, mom) in signature.1.iter().zip_eq(external_momenta) {
            if sign != 0 {
                cut_momentum += mom.multiply_sign(sign);
            }
        }
        cut_momentum
    }

    pub fn create_caches<T: FloatLike>(&self) -> Vec<Vec<Vec<LTDCache<T>>>> {
        let mut caches = vec![];
        for cutkosky_cuts in &self.cutkosky_cuts {
            let mut lim_cache = vec![];
            for uv_limit in &cutkosky_cuts.uv_limits {
                let mut diag_cache = vec![];
                for d in &uv_limit.diagrams {
                    diag_cache.push(LTDCache::<T>::new(d));
                }
                lim_cache.push(diag_cache);
            }
            caches.push(lim_cache);
        }
        caches
    }

    #[inline]
    /// Kahlen function.
    fn lambda<T: FloatLike>(s: T, sqr_mass_a: T, sqr_mass_b: T) -> T {
        s * s + sqr_mass_a * sqr_mass_a + sqr_mass_b * sqr_mass_b
            - Into::<T>::into(2.) * s * sqr_mass_a
            - Into::<T>::into(2.) * sqr_mass_b * sqr_mass_a
            - Into::<T>::into(2.) * s * sqr_mass_b
    }

    /// Solve the momentum conservation delta using Newton's method. In the case of massless propagators and an external momentum with 0 spatial part,
    /// it will take one step to find the solution. This function returns the scaling parameter and its Jacobian.
    pub fn find_scaling<T: FloatLike>(
        cutkosky_cuts: &CutkoskyCuts,
        external_momenta: &[LorentzVector<T>],
        loop_momenta: &[LorentzVector<T>],
        incoming_energy: T,
        debug_level: usize,
    ) -> Option<[(T, T); 2]> {
        // determine an overestimate of the t that solves the energy constraint
        // then -t and t should give two different solutions or do not converge
        let mut t_start = T::zero();
        let mut sum_k = T::zero();
        for cut in &cutkosky_cuts.cuts {
            if cut.level != 0 {
                continue;
            }

            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            let shift = utils::evaluate_signature(&cut.signature.1, external_momenta);
            let k_norm_sq = k.spatial_squared();
            t_start += Float::abs(k.spatial_dot(&shift)) / k_norm_sq;
            sum_k += k_norm_sq.sqrt();
        }

        t_start += incoming_energy / sum_k;

        // find the two solutions
        let mut solutions = [(T::zero(), T::zero()); 2];
        for (i, &mut mut t) in [-t_start, t_start].iter_mut().enumerate() {
            for it in 0..20 {
                let mut f = incoming_energy;
                let mut df = T::zero();

                for cut in &cutkosky_cuts.cuts {
                    if cut.level != 0 {
                        continue;
                    }
                    let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
                    let shift = utils::evaluate_signature(&cut.signature.1, external_momenta);
                    let energy =
                        ((k * t + shift).spatial_squared() + Into::<T>::into(cut.m_squared)).sqrt();
                    f -= energy;
                    df -= (t * k.spatial_squared() + k.spatial_dot(&shift)) / energy;
                }

                if debug_level > 4 {
                    println!("  | t{} finder: f={}, df={}, t={}", i, f, df, t);
                }

                if Float::abs(f) < Into::<T>::into(1e-14) * incoming_energy {
                    if debug_level > 2 {
                        println!("  | t{} = {}", i, t);
                    }

                    solutions[i] = (t, Float::abs(df).inv());
                    break;
                }

                if it == 19 {
                    if debug_level > 2 {
                        println!(
                            "  | no convergence after {} iterations: f={}, df={}, t={}",
                            i, f, df, t
                        );
                    }
                    return None;
                }

                t = t - f / df;
            }
        }

        if Float::abs(solutions[0].0 - solutions[1].0) < Into::<T>::into(1e-12) * incoming_energy {
            panic!(
                "Found the same scaling solution twice: {} for t={} and t={} for k={:?}, ext={:?}",
                solutions[0].0, -t_start, t_start, loop_momenta, external_momenta
            );
        }

        Some(solutions)
    }

    pub fn evaluate<'a, T: FloatLike>(
        &mut self,
        x: &'a [f64],
        cache: &mut [Vec<Vec<LTDCache<T>>>],
        mut event_manager: Option<&mut EventManager>,
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
            let (l_space, jac) = Topology::parameterize(
                &x[i * 3..(i + 1) * 3],
                self.e_cm_squared,
                i,
                &self.settings,
            );

            // there could be some rounding here
            let rot = &self.rotation_matrix;
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

        let result = self.evaluate_mom(&k[..self.n_loops], cache, &mut event_manager) * jac_para;

        if let Some(em) = &mut event_manager {
            if em.track_events {
                for e in em.event_buffer.iter_mut() {
                    e.integrand *= jac_para.to_f64().unwrap();
                }
            }
        }

        // NOTE: there is no unique k_def anymore. it depends on the cut
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]> = ArrayVec::default();
        for l in &k[..self.n_loops] {
            k_def.push(l.map(|x| Complex::new(x, T::zero())));
        }

        (x, k_def, jac_para, Complex::one(), result)
    }

    pub fn evaluate_mom<T: FloatLike>(
        &mut self,
        loop_momenta: &[LorentzVector<T>],
        caches: &mut [Vec<Vec<LTDCache<T>>>],
        event_manager: &mut Option<&mut EventManager>,
    ) -> Complex<T> {
        debug_assert_eq!(
            loop_momenta.len(),
            self.n_loops,
            "The number of loop momenta is wrong."
        );

        let mut external_momenta: ArrayVec<[LorentzVector<T>; MAX_LOOP]> = self
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4]; // FIXME: bound may be too small
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_LOOP + 4];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_LOOP + 4];
        let mut result = Complex::zero();
        for cut_index in 0..self.cutkosky_cuts.len() {
            let cutkosky_cuts = &mut self.cutkosky_cuts[cut_index];

            let scaling_solutions = if self.settings.cross_section.do_rescaling {
                let incoming_energy = external_momenta[..self.n_incoming_momenta]
                    .iter()
                    .map(|m| m.t)
                    .sum();
                SquaredTopology::find_scaling(
                    cutkosky_cuts,
                    &external_momenta[..self.external_momenta.len()],
                    loop_momenta,
                    incoming_energy,
                    self.settings.general.debug,
                )
            } else {
                // we do not scale, so give one positive solution that is just 1
                Some([(-T::one(), T::one()), (T::one(), T::one())])
            };

            if scaling_solutions.is_none() {
                if self.settings.general.debug >= 1 {
                    println!(
                        "Phase space point has no solutions for cut {:?}:",
                        cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                    );
                }
                continue;
            }

            if self.settings.general.debug >= 1 {
                println!(
                    "Cut {}:",
                    cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                );
            }

            for &(scaling, scaling_jac) in &scaling_solutions.unwrap() {
                if scaling < T::zero() {
                    // we ignore all negative solutions
                    continue;
                }

                result += self.evaluate_cut(
                    loop_momenta,
                    &mut cut_momenta,
                    &mut external_momenta,
                    &mut rescaled_loop_momenta,
                    &mut subgraph_loop_momenta,
                    &mut k_def[..self.n_loops + 1],
                    &mut caches[cut_index],
                    event_manager,
                    cut_index,
                    scaling,
                    scaling_jac,
                );
            }
        }

        if let Some(em) = event_manager {
            if em.track_events {
                for e in em.event_buffer.iter_mut() {
                    e.integrand *= self.overall_numerator;
                }
            }
        }

        result *= Into::<T>::into(self.overall_numerator);

        if self.settings.general.debug >= 1 {
            println!("Final result = {:e}", result);
        }

        result
    }

    pub fn evaluate_cut<T: FloatLike>(
        &mut self,
        loop_momenta: &[LorentzVector<T>],
        cut_momenta: &mut [LorentzVector<T>],
        external_momenta: &mut ArrayVec<[LorentzVector<T>; MAX_LOOP]>,
        rescaled_loop_momenta: &mut [LorentzVector<T>],
        subgraph_loop_momenta: &mut [LorentzVector<T>],
        k_def: &mut [LorentzVector<Complex<T>>],
        cache: &mut [Vec<LTDCache<T>>],
        event_manager: &mut Option<&mut EventManager>,
        cut_index: usize,
        scaling: T,
        scaling_jac: T,
    ) -> Complex<T> {
        let cutkosky_cuts = &mut self.cutkosky_cuts[cut_index];

        let mut cut_energies_summed = T::zero();
        let mut scaling_result = Complex::one();

        // evaluate the cuts with the proper scaling
        for (cut_mom, cut) in cut_momenta[..cutkosky_cuts.cuts.len()]
            .iter_mut()
            .zip(cutkosky_cuts.cuts.iter())
        {
            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            let shift = utils::evaluate_signature(
                &cut.signature.1,
                &external_momenta[..self.external_momenta.len()],
            );
            *cut_mom = k * scaling + shift;
            let energy = (cut_mom.spatial_squared() + Into::<T>::into(cut.m_squared)).sqrt();
            cut_mom.t = energy.multiply_sign(cut.sign);
            scaling_result *= num::Complex::new(T::zero(), -<T as FloatConst>::PI() / energy); // add (-2 pi i)/(2E) for every cut

            if cut.level == 0 {
                // only Cutkosky cuts constribute to the energy delta
                cut_energies_summed += energy;
            }
        }

        if self.settings.general.debug >= 1 {
            println!("  | 1/Es = {}", scaling_result);
            println!("  | q0 = {}", cut_energies_summed);
            println!("  | scaling = {}", scaling);
            println!("  | scaling_jac = {}", scaling_jac);
        }

        if self.settings.cross_section.do_rescaling {
            // h is any function that integrates to 1
            let h = (-Float::powi(
                scaling / Into::<T>::into(self.settings.cross_section.rescaling_function_spread),
                2,
            ))
            .exp();
            let h_norm = <T as FloatConst>::PI().sqrt() / Into::<T>::into(2.0)
                * Into::<T>::into(self.settings.cross_section.rescaling_function_spread);

            scaling_result *=
                scaling_jac * Float::powi(scaling, self.n_loops as i32 * 3) * h / h_norm;

            // rescale the loop momenta
            for (rlm, lm) in rescaled_loop_momenta[..self.n_loops]
                .iter_mut()
                .zip_eq(loop_momenta)
            {
                *rlm = lm * scaling;
            }
        } else {
            // set 0th component of external momenta to the sum of the cuts
            // NOTE: we assume 1->1 here and that it is outgoing
            external_momenta[0].t = cut_energies_summed;
            external_momenta[1].t = cut_energies_summed;
        }

        let e_cm_sq = if self.n_incoming_momenta == 2 {
            (external_momenta[0] + external_momenta[1]).square()
        } else {
            external_momenta[0].square()
        };

        if let Some(em) = event_manager {
            if !em.add_event(
                &self.external_momenta,
                &cut_momenta[..cutkosky_cuts.cuts.len()],
                &cutkosky_cuts.cuts,
                &self.rotation_matrix,
            ) {
                // the event was cut
                return Complex::zero();
            }
        }

        // now apply the same procedure for all uv limits
        let mut diag_and_num_contributions = Complex::zero();
        for (cut_uv_limit, diag_cache) in cutkosky_cuts.uv_limits.iter_mut().zip(cache) {
            // set the shifts, which are expressed in the cut basis
            for subgraph in &mut cut_uv_limit.diagrams {
                for ll in &mut subgraph.loop_lines {
                    for p in &mut ll.propagators {
                        p.q = SquaredTopology::evaluate_signature(
                            &p.parametric_shift,
                            &external_momenta[..self.external_momenta.len()],
                            &cut_momenta[..cutkosky_cuts.cuts.len()],
                        )
                        .map(|x| x.to_f64().unwrap());
                    }
                }
                subgraph.e_cm_squared = e_cm_sq.to_f64().unwrap();

                subgraph.update_ellipsoids();

                if self.settings.general.deformation_strategy == DeformationStrategy::Fixed {
                    subgraph.fixed_deformation =
                        subgraph.determine_ellipsoid_overlap_structure(true);

                    if self.settings.general.debug > 0 {
                        // check if the overlap structure makes sense
                        subgraph.check_fixed_deformation();
                    }
                }
            }

            // for the evaluation of the numerator we need complex loop momenta of the supergraph.
            // the order is: cut momenta, momenta left graph, momenta right graph
            // FIXME: this requires that the last level 0 cutkosky cut is at the end of the list
            for (kd, cut_mom) in k_def
                .iter_mut()
                .zip(&cut_momenta[..cutkosky_cuts.cuts.len() - 1])
            {
                *kd = cut_mom.to_complex(true);
            }

            // compute the deformation vectors
            let mut k_def_index = cutkosky_cuts.cuts.len() - 1;
            for ((subgraph, subgraph_cache), &conjugate_deformation) in cut_uv_limit
                .diagrams
                .iter_mut()
                .zip_eq(diag_cache.iter_mut())
                .zip_eq(cut_uv_limit.conjugate_deformation.iter())
            {
                // do the loop momentum map, which is expressed in the loop momentum basis
                // the time component should not matter here
                for (slm, lmm) in subgraph_loop_momenta[..subgraph.n_loops]
                    .iter_mut()
                    .zip_eq(&subgraph.loop_momentum_map)
                {
                    *slm = SquaredTopology::evaluate_signature(
                        lmm,
                        &external_momenta[..self.external_momenta.len()],
                        if self.settings.cross_section.do_rescaling {
                            &rescaled_loop_momenta[..self.n_loops]
                        } else {
                            loop_momenta
                        },
                    );
                }

                let (kappas, jac_def) =
                    subgraph.deform(&subgraph_loop_momenta[..subgraph.n_loops], subgraph_cache);

                for (lm, kappa) in subgraph_loop_momenta[..subgraph.n_loops]
                    .iter()
                    .zip(&kappas)
                {
                    k_def[k_def_index] = if conjugate_deformation {
                        // take the complex conjugate of the deformation
                        lm.map(|x| Complex::new(x, T::zero()))
                            - kappa.map(|x| Complex::new(T::zero(), x))
                    } else {
                        lm.map(|x| Complex::new(x, T::zero()))
                            + kappa.map(|x| Complex::new(T::zero(), x))
                    };
                    k_def_index += 1;
                }

                if conjugate_deformation {
                    scaling_result *= jac_def.conj();
                } else {
                    scaling_result *= jac_def;
                }

                if subgraph
                    .compute_complex_cut_energies(
                        &k_def[k_def_index - subgraph.n_loops..k_def_index],
                        subgraph_cache,
                    )
                    .is_err()
                {
                    panic!("NaN on cut energy");
                }
            }

            // evaluate the cut numerator with the spatial parts of all loop momenta and
            // the energy parts of the Cutkosky cuts.
            // Store the reduced numerator in the left graph cache for now
            // TODO: this is impractical
            cut_uv_limit.numerator.evaluate_reduced_in_lb(
                &k_def,
                cutkosky_cuts.cuts.len() - 1,
                &mut diag_cache[0],
                0,
            );

            let left_cache = &mut diag_cache[0];
            mem::swap(
                &mut left_cache.reduced_coefficient_lb_supergraph[0],
                &mut left_cache.reduced_coefficient_lb[0],
            );

            let mut supergraph_coeff = mem::replace(
                &mut diag_cache[0].reduced_coefficient_lb_supergraph[0],
                vec![],
            );

            // evaluate the left and right graphs for every monomial in the numerator
            for (coeff, powers) in supergraph_coeff
                .iter()
                .zip(&cut_uv_limit.numerator.reduced_coefficient_index_to_powers)
            {
                if coeff.is_zero() {
                    continue;
                }

                // only consider the coefficients that have no powers in the cutkosky cuts
                // TODO: make a more efficient way of skipping the other contributions
                if powers[..cutkosky_cuts.cuts.len() - 1]
                    .iter()
                    .any(|p| *p != 0)
                {
                    continue;
                }

                let mut num_result = *coeff;

                let mut def_mom_index = cutkosky_cuts.cuts.len() - 1;
                for (subgraph, subgraph_cache) in cut_uv_limit
                    .diagrams
                    .iter_mut()
                    .zip_eq(diag_cache.iter_mut())
                {
                    // build the subgraph numerator
                    // TODO: build the reduced numerator directly
                    // TODO: cache, such that we never compute a subgraph with the same numerator twice
                    for n in &mut subgraph.numerator.coefficients {
                        *n = Complex::zero();
                    }

                    if subgraph.n_loops == 0 {
                        subgraph.numerator.coefficients[0] = Complex::one();
                    } else {
                        // now set the coefficient to 1 for the current monomial in the subgraph
                        // by mapping the powers in the reduced numerator back to the power pattern
                        // of the complete numerator
                        // TODO: use binary search
                        let mut subgraph_powers = [0; 10]; // TODO: make max rank a constant
                        let mut power_index = 0;
                        for (lmi, p) in powers[def_mom_index..def_mom_index + subgraph.n_loops]
                            .iter()
                            .enumerate()
                        {
                            for _ in 0..*p {
                                subgraph_powers[power_index] = lmi * 4;
                                power_index += 1;
                            }
                        }
                        if power_index == 0 {
                            subgraph.numerator.coefficients[0] = Complex::one();
                        } else {
                            subgraph.numerator.coefficients[subgraph
                                .numerator
                                .sorted_linear
                                .iter()
                                .position(|x| x[..] == subgraph_powers[..power_index])
                                .unwrap()
                                + 1] += Complex::one();
                        }
                    }

                    let mut res = if subgraph.loop_lines.len() > 0 {
                        subgraph
                            .evaluate_all_dual_integrands::<T>(
                                if subgraph.n_loops != 0 {
                                    &mut k_def[def_mom_index..def_mom_index + subgraph.n_loops]
                                } else {
                                    &mut []
                                },
                                subgraph_cache,
                            )
                            .unwrap()
                    } else {
                        // if the graph has no propagators, it is one and not zero
                        Complex::one()
                    };

                    res *= utils::powi(
                        num::Complex::new(
                            T::zero(),
                            Into::<T>::into(-2.) * <T as FloatConst>::PI(),
                        ),
                        subgraph.n_loops,
                    ); // factor of delta cut

                    num_result *= res;

                    if self.settings.general.debug >= 1 {
                        println!(
                            "  | res {} ({}l) = {:e}",
                            subgraph.name, subgraph.n_loops, res
                        );
                    }

                    def_mom_index += subgraph.n_loops;
                }
                diag_and_num_contributions += num_result;
            }

            mem::swap(
                &mut supergraph_coeff,
                &mut diag_cache[0].reduced_coefficient_lb_supergraph[0],
            );
        }

        scaling_result *= diag_and_num_contributions;

        scaling_result *= utils::powi(
            num::Complex::new(
                Into::<T>::into(1.)
                    / <T as Float>::powi(Into::<T>::into(2.) * <T as FloatConst>::PI(), 4),
                T::zero(),
            ),
            self.n_loops,
        );

        // multiply the flux factor
        scaling_result /= if self.n_incoming_momenta == 2 {
            Into::<T>::into(2.)
                * (SquaredTopology::lambda(
                    (external_momenta[0] + external_momenta[1]).square(),
                    external_momenta[0].square(),
                    external_momenta[1].square(),
                ))
                .sqrt()
        } else {
            let e_cm = external_momenta[0].square().sqrt();
            e_cm * Into::<T>::into(2.0)
        };

        if self.settings.cross_section.picobarns {
            // return a weight in picobarns (from GeV^-2)
            scaling_result *= Into::<T>::into(0.389379304e9);
        }

        // divide by the symmetry factor of the final state
        //scaling_result /= Into::<T>::into(cutkosky_cuts.symmetry_factor); // TODO: no more sym factor!

        if let Some(em) = event_manager {
            // set the event result
            if em.track_events {
                em.event_buffer.last_mut().unwrap().integrand = Complex::new(
                    scaling_result.re.to_f64().unwrap(),
                    scaling_result.im.to_f64().unwrap(),
                );
            }
        }

        if self.settings.general.debug >= 1 {
            println!("  | scaling res = {:e}", scaling_result);
        }

        scaling_result
    }
}

impl IntegrandImplementation for SquaredTopology {
    type Cache = SquaredTopologyCache;

    /// Create a rotated version of this squared topology. The axis needs to be normalized.
    fn rotate(&self, angle: float, axis: (float, float, float)) -> SquaredTopology {
        let cos_t = angle.cos();
        let sin_t = angle.sin();
        let cos_t_bar = float::one() - angle.cos();

        let rot_matrix: [[float; 3]; 3] = [
            [
                cos_t + axis.0 * axis.0 * cos_t_bar,
                axis.0 * axis.1 * cos_t_bar - axis.2 * sin_t,
                axis.0 * axis.2 * cos_t_bar + axis.1 * sin_t,
            ],
            [
                axis.0 * axis.1 * cos_t_bar + axis.2 * sin_t,
                cos_t + axis.1 * axis.1 * cos_t_bar,
                axis.1 * axis.2 * cos_t_bar - axis.0 * sin_t,
            ],
            [
                axis.0 * axis.2 * cos_t_bar - axis.1 * sin_t,
                axis.1 * axis.2 * cos_t_bar + axis.0 * sin_t,
                cos_t + axis.2 * axis.2 * cos_t_bar,
            ],
        ];

        let mut rotated_topology = self.clone();
        rotated_topology.rotation_matrix = rot_matrix.clone();

        for e in &mut rotated_topology.external_momenta {
            let old_x = float::from_f64(e.x).unwrap();
            let old_y = float::from_f64(e.y).unwrap();
            let old_z = float::from_f64(e.z).unwrap();
            e.x = (rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z)
                .to_f64()
                .unwrap();
            e.y = (rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z)
                .to_f64()
                .unwrap();
            e.z = (rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z)
                .to_f64()
                .unwrap();
        }

        for cut in rotated_topology.cutkosky_cuts.iter_mut() {
            for uv_limit in cut.uv_limits.iter_mut() {
                uv_limit.numerator = uv_limit.numerator.rotate(rot_matrix);
            }
        }

        rotated_topology
    }

    #[inline]
    fn evaluate_float<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut SquaredTopologyCache,
        event_manager: Option<&mut EventManager>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<float>>; MAX_LOOP]>,
        float,
        Complex<float>,
        Complex<float>,
    ) {
        self.evaluate(x, cache.get(), event_manager)
    }

    #[inline]
    fn evaluate_f128<'a>(
        &mut self,
        x: &'a [f64],
        cache: &mut SquaredTopologyCache,
        event_manager: Option<&mut EventManager>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<f128>>; MAX_LOOP]>,
        f128,
        Complex<f128>,
        Complex<f128>,
    ) {
        self.evaluate(x, cache.get(), event_manager)
    }

    #[inline]
    fn create_cache(&self) -> SquaredTopologyCache {
        SquaredTopologyCache {
            float_cache: self.create_caches(),
            quad_cache: self.create_caches(),
        }
    }
}
