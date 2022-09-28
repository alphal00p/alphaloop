use crate::overlap::DeformationOverlap;
use crate::squared_topologies::DeformationSubgraph;
use crate::{
    AdditiveMode, ExpansionCheckStrategy, FloatLike, IRHandling, OverallDeformationScaling,
    PoleCheckStrategy, MAX_LOOP,
};
use hyperdual::Hyperdual;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use num::Complex;
use num_traits::{Float, One, Pow, Zero};
use utils::Signum;

use crate::utils;

macro_rules! deform {
    ($self:ident, $loop_momenta:ident, $cache:ident, $fn:ident, $output:ident, $($loops:expr),*) => (
        match $self.n_loops {
            0 => {}
            $($loops => {
                let mut r = [LorentzVector::default(); MAX_LOOP];
                for j in 0..$loops {
                    r[j] = $loop_momenta[j].map(|x| Hyperdual::<T, {$loops * 3 + 1}>::from_real(x));

                    for i in 0..3 {
                        r[j][i + 1][i + 1 + j * 3] = T::one();
                    }
                }

                let (kappas, lambda) = $self.$fn::<T, {$loops * 3 + 1}>(&mut r[..$loops], $cache);

                for j in 0..$loops {
                    for k in 0..3 {
                        linearize(&kappas[j][k + 1], &mut $output.0[j*3*{$loops * 3 + 1} + k*{$loops * 3 + 1}..])
                    }
                }

                linearize(&lambda, &mut $output.1);
            }
            )+
            _ => panic!("No supported dual: compile with `higher_orders` feature?"),
        }
    )
}

fn linearize<T: FloatLike, const N: usize>(dual: &Hyperdual<T, N>, output: &mut [T]) {
    for (o, d) in output.iter_mut().zip(dual.iter()) {
        *o = *d;
    }
}

#[derive(Debug, Clone)]
/// A cache for objects needed for the deformation
pub struct CutInfo<T: FloatLike, const N: usize> {
    pub id: usize,
    pub momentum: LorentzVector<Hyperdual<T, N>>,
    pub real_energy: Hyperdual<T, N>,
    pub spatial_and_mass_sq: Hyperdual<T, N>,
    pub spatial_and_uv_mass_sq: Hyperdual<T, N>,
    pub shift: LorentzVector<Hyperdual<T, N>>,
    pub kappa: LorentzVector<Hyperdual<T, N>>,
    pub kappa_sq: Hyperdual<T, N>,
    pub kappa_dot_mom: Hyperdual<T, N>,
    pub m_sq: Hyperdual<T, N>,
    pub a: Hyperdual<T, N>,
    pub b: Hyperdual<T, N>,
    pub c: Hyperdual<T, N>,
}

impl<T: FloatLike, const N: usize> Default for CutInfo<T, N> {
    fn default() -> CutInfo<T, N> {
        CutInfo {
            id: 0,
            momentum: LorentzVector::default(),
            real_energy: Hyperdual::default(),
            spatial_and_mass_sq: Hyperdual::default(),
            spatial_and_uv_mass_sq: Hyperdual::default(),
            shift: LorentzVector::default(),
            kappa: LorentzVector::default(),
            kappa_sq: Hyperdual::default(),
            kappa_dot_mom: Hyperdual::default(),
            m_sq: Hyperdual::default(),
            a: Hyperdual::default(),
            b: Hyperdual::default(),
            c: Hyperdual::default(),
        }
    }
}

#[derive(Debug, Clone)]
/// A cache for objects needed during LTD computation
pub struct DeformationCacheI<T: FloatLike, const N: usize> {
    pub ellipsoid_eval: Vec<Option<Hyperdual<T, N>>>,
    pub ellipsoid_normal_norm_eval: Vec<Option<Hyperdual<T, N>>>,
    pub deform_dirs: Vec<LorentzVector<Hyperdual<T, N>>>,
    pub cut_energies: Vec<Hyperdual<T, N>>,
    pub cut_info: Vec<CutInfo<T, N>>,
}

impl<T: FloatLike, const N: usize> Default for DeformationCacheI<T, N> {
    fn default() -> DeformationCacheI<T, N> {
        DeformationCacheI {
            ellipsoid_eval: vec![],
            ellipsoid_normal_norm_eval: vec![],
            deform_dirs: vec![],
            cut_energies: vec![],
            cut_info: vec![],
        }
    }
}

impl<T: FloatLike, const N: usize> DeformationCacheI<T, N> {
    fn new(
        num_loops: usize,
        num_surfaces: usize,
        num_propagators: usize,
    ) -> DeformationCacheI<T, N> {
        DeformationCacheI {
            ellipsoid_eval: vec![None; num_surfaces],
            ellipsoid_normal_norm_eval: vec![None; num_surfaces],
            deform_dirs: vec![LorentzVector::default(); num_surfaces * num_loops],
            cut_energies: vec![Hyperdual::default(); num_propagators],
            cut_info: vec![CutInfo::default(); num_propagators],
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct DeformationCache<T: FloatLike> {
    one_loop: DeformationCacheI<T, 4>,
    two_loop: DeformationCacheI<T, 7>,
    three_loop: DeformationCacheI<T, 10>,
    #[cfg(feature = "higher_loops")]
    four_loop: DeformationCacheI<T, 13>,
    #[cfg(feature = "higher_loops")]
    five_loop: DeformationCacheI<T, 16>,
    #[cfg(feature = "higher_loops")]
    six_loop: DeformationCacheI<T, 19>,
    pub complex_cut_energies: Vec<Complex<T>>,
    pub complex_prop_spatial: Vec<Complex<T>>,
    pub complex_ellipsoids: Vec<Vec<Complex<T>>>,
    pub overall_lambda: T, // used to log the minimum
    pub propagators_eval: Vec<Complex<T>>,
    pub propagator_powers: Vec<usize>,
    pub overlap_structure: Vec<DeformationOverlap>,
    pub deformation_jacobian_matrix: Vec<Complex<T>>,
}

impl<T: FloatLike> DeformationCache<T> {
    pub fn new(topo: &DeformationSubgraph) -> DeformationCache<T> {
        let num_propagators = topo.propagators.len();
        DeformationCache {
            one_loop: DeformationCacheI::<T, 4>::new(
                topo.n_loops,
                topo.thresholds.len(),
                num_propagators,
            ),
            two_loop: DeformationCacheI::<T, 7>::new(
                topo.n_loops,
                topo.thresholds.len(),
                num_propagators,
            ),
            three_loop: DeformationCacheI::<T, 10>::new(
                topo.n_loops,
                topo.thresholds.len(),
                num_propagators,
            ),
            #[cfg(feature = "higher_loops")]
            four_loop: DeformationCacheI::<T, 13>::new(
                topo.n_loops,
                topo.thresholds.len(),
                num_propagators,
            ),
            #[cfg(feature = "higher_loops")]
            five_loop: DeformationCacheI::<T, 16>::new(
                topo.n_loops,
                topo.thresholds.len(),
                num_propagators,
            ),
            #[cfg(feature = "higher_loops")]
            six_loop: DeformationCacheI::<T, 19>::new(
                topo.n_loops,
                topo.thresholds.len(),
                num_propagators,
            ),
            complex_cut_energies: vec![Complex::default(); num_propagators],
            complex_prop_spatial: vec![Complex::default(); num_propagators],
            complex_ellipsoids: vec![vec![Complex::default(); num_propagators]; num_propagators],
            overall_lambda: T::zero(),
            propagators_eval: vec![Complex::zero(); num_propagators],
            propagator_powers: vec![1; num_propagators],
            overlap_structure: vec![],
            deformation_jacobian_matrix: vec![],
        }
    }
}

pub trait CacheSelector<T: FloatLike, const N: usize> {
    fn get_cache(&self) -> &DeformationCacheI<T, N>;
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, N>;
}

impl<T: FloatLike> CacheSelector<T, 4> for DeformationCache<T> {
    #[inline]
    fn get_cache(&self) -> &DeformationCacheI<T, 4> {
        &self.one_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, 4> {
        &mut self.one_loop
    }
}

impl<T: FloatLike> CacheSelector<T, 7> for DeformationCache<T> {
    #[inline]
    fn get_cache(&self) -> &DeformationCacheI<T, 7> {
        &self.two_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, 7> {
        &mut self.two_loop
    }
}

impl<T: FloatLike> CacheSelector<T, 10> for DeformationCache<T> {
    #[inline]
    fn get_cache(&self) -> &DeformationCacheI<T, 10> {
        &self.three_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, 10> {
        &mut self.three_loop
    }
}

#[cfg(feature = "higher_loops")]
impl<T: FloatLike> CacheSelector<T, 13> for DeformationCache<T> {
    #[inline]
    fn get_cache(&self) -> &DeformationCacheI<T, 13> {
        &self.four_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, 13> {
        &mut self.four_loop
    }
}

#[cfg(feature = "higher_loops")]
impl<T: FloatLike> CacheSelector<T, 16> for DeformationCache<T> {
    #[inline]
    fn get_cache(&self) -> &DeformationCacheI<T, 16> {
        &self.five_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, 16> {
        &mut self.five_loop
    }
}

#[cfg(feature = "higher_loops")]
impl<T: FloatLike> CacheSelector<T, 19> for DeformationCache<T> {
    #[inline]
    fn get_cache(&self) -> &DeformationCacheI<T, 19> {
        &self.six_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut DeformationCacheI<T, 19> {
        &mut self.six_loop
    }
}

#[derive(Default)]
pub struct LTDCacheAllPrecisions {
    float_cache: DeformationCache<f64>,
    quad_cache: DeformationCache<f128::f128>,
}

impl LTDCacheAllPrecisions {
    pub fn new(topo: &DeformationSubgraph) -> LTDCacheAllPrecisions {
        LTDCacheAllPrecisions {
            float_cache: DeformationCache::new(topo),
            quad_cache: DeformationCache::new(topo),
        }
    }
}

pub trait CachePrecisionSelector<T: FloatLike> {
    fn get(&mut self) -> &mut DeformationCache<T>;
}

impl CachePrecisionSelector<f64> for LTDCacheAllPrecisions {
    #[inline]
    fn get(&mut self) -> &mut DeformationCache<f64> {
        &mut self.float_cache
    }
}

impl CachePrecisionSelector<f128::f128> for LTDCacheAllPrecisions {
    #[inline]
    fn get(&mut self) -> &mut DeformationCache<f128::f128> {
        &mut self.quad_cache
    }
}

impl DeformationSubgraph {
    #[inline]
    pub fn get_expansion_threshold(&self) -> f64 {
        match self.settings.deformation.scaling.expansion_check_strategy {
            ExpansionCheckStrategy::Ratio => {
                self.settings.deformation.scaling.expansion_threshold.abs()
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
        cache: &mut DeformationCache<T>,
    ) -> Hyperdual<T, N>
    where
        DeformationCache<T>: CacheSelector<T, N>,
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

        for f in &self.propagators {
            let mut kappa_cut = LorentzVector::default();
            for (kappa, &sign) in kappas.iter().zip_eq(f.amp_sig.iter()) {
                kappa_cut += kappa.multiply_sign(sign);
            }

            let info = &mut cache.cut_info[f.id];
            info.kappa = kappa_cut;
            info.kappa_sq = kappa_cut.spatial_squared_impr();
            info.kappa_dot_mom = kappa_cut.spatial_dot_impr(&info.momentum);

            if self.settings.deformation.scaling.pole_check_strategy != PoleCheckStrategy::None {
                info.a = (info.kappa_sq * info.real_energy.powi(2) - info.kappa_dot_mom.powi(2))
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
                let c =
                    Hyperdual::<T, N>::from_real(Into::<T>::into(self.get_expansion_threshold()));

                let lambda_exp_sq = match self.settings.deformation.scaling.expansion_check_strategy
                {
                    ExpansionCheckStrategy::FirstLambdaOrder => {
                        let lambda_exp = c * info.spatial_and_mass_sq / info.kappa_dot_mom.abs(); // note: not holomorphic
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
                                * T::convert_from(&self.e_cm_squared)
                                * T::convert_from(&self.e_cm_squared));
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

        // TODO: convert surfaces to lmb
        if self.settings.deformation.scaling.pole_check_strategy != PoleCheckStrategy::None {
            // process all existing and non-existing ellipsoids and pinches
            for s in &self.thresholds {
                let mut a = Hyperdual::<T, N>::from_real(T::zero());
                let mut b = Hyperdual::<T, N>::from_real(T::zero());
                let mut c = Hyperdual::<T, N>::from_real(T::zero());
                for f in &s.foci {
                    let info = &cache.cut_info[*f];
                    a += info.a;
                    b += info.b;
                    c += info.c + info.shift.t;
                }

                let prop_lambda_sq = match self.settings.deformation.scaling.pole_check_strategy {
                    // only works at one-loop
                    PoleCheckStrategy::Exact => {
                        assert_eq!(self.n_loops, 1);

                        let cut_info = &cache.cut_info[s.foci[0]];
                        let on_shell_info = &cache.cut_info[s.foci[1]];

                        // now determine the exact solution for one-loop surfaces
                        let p0 = on_shell_info.shift.t - cut_info.shift.t;
                        let a = cut_info.spatial_and_mass_sq
                            - on_shell_info.spatial_and_mass_sq
                            - p0.powi(2);

                        let b = kappas[0]
                            .spatial_dot(&(cut_info.momentum - on_shell_info.momentum))
                            * Into::<T>::into(2.);

                        let aa =
                            kappas[0].spatial_squared() * p0.powi(2) * Into::<T>::into(4.) - b * b;
                        let bb = (a * b
                            - kappas[0].spatial_dot(&on_shell_info.momentum)
                                * p0.powi(2)
                                * Into::<T>::into(4.))
                            * Into::<T>::into(2.);
                        let cc = a * a
                            - on_shell_info.spatial_and_mass_sq * p0.powi(2) * Into::<T>::into(4.);

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
                        let prop_lambda_sq = DeformationSubgraph::compute_lambda_factor(x, y);
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

                        let t_c =
                            t_out + b / Into::<T>::into(self.settings.deformation.scaling.theta_c);

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

        if sigma.is_zero() {
            lambda_sq.sqrt()
        } else {
            (smooth_min_num / smooth_min_den).sqrt()
        }
    }

    fn deform_fixed<T: FloatLike, const N: usize>(
        &self,
        loop_momenta: &[LorentzVector<Hyperdual<T, N>>],
        cache: &mut DeformationCache<T>,
    ) -> [LorentzVector<Hyperdual<T, N>>; MAX_LOOP]
    where
        DeformationCache<T>: CacheSelector<T, N>,
    {
        let mut kappas = [LorentzVector::default(); MAX_LOOP];
        let mut kappa_source = [LorentzVector::default(); MAX_LOOP];
        let mij_min_sq = Into::<T>::into(self.compute_min_mij().powi(2));
        let mut mij_sq =
            Into::<T>::into(self.settings.deformation.fixed.m_ij.abs().powi(2)) * mij_min_sq;

        let overlap_structure = std::mem::take(&mut cache.overlap_structure);

        if overlap_structure.is_empty() {
            return kappas;
        }

        let cache_i = cache.get_cache_mut();

        // evaluate all excluded ellipsoids and pinches
        for (i, surf) in self.thresholds.iter().enumerate() {
            if !surf.exists {
                continue;
            }

            let mut cut_energy = Hyperdual::from_real(Into::<T>::into(surf.shift.t));
            let mut cut_dirs = [LorentzVector::default(); MAX_LOOP];

            for f in &surf.foci {
                cut_energy += cache_i.cut_info[*f].real_energy;

                // TODO: check
                for (cut_dir, c) in cut_dirs.iter_mut().zip(&self.propagators[*f].amp_sig) {
                    *cut_dir += cache_i.cut_info[*f].momentum
                        / cache_i.cut_info[*f].real_energy.multiply_sign(*c);
                }
            }

            let mut norm_der = Hyperdual::<T, N>::zero();
            for c in &cut_dirs[..self.n_loops] {
                norm_der += c.spatial_squared_impr();
            }

            cache_i.ellipsoid_eval[i] = Some(cut_energy);
            cache_i.ellipsoid_normal_norm_eval[i] = Some(norm_der);
        }

        for k in &mut kappa_source {
            *k = LorentzVector::default();
        }

        for (overlap_index, d) in overlap_structure.iter().enumerate() {
            if self.settings.general.debug > 2 {
                println!("Deformation contribution from overlap {}:", overlap_index);
            }

            let mut s = Hyperdual::<T, N>::one();
            let mut softmin_num = Hyperdual::<T, N>::zero();
            let mut softmin_den = Hyperdual::<T, N>::zero();

            // dampen every excluded threshold
            for &surf_id in &d.excluded_surface_indices {
                let t = self.thresholds.iter().find(|t| t.id == surf_id).unwrap();
                if t.id < self.settings.deformation.fixed.m_ijs.len() {
                    mij_sq =
                        Into::<T>::into(self.settings.deformation.fixed.m_ijs[t.id].abs().powi(2))
                            * mij_min_sq;
                }
                // t is the weighing factor that is 0 if we are on both cut_i and cut_j
                // at the same time and goes to 1 otherwise

                // The surface shift should never be zero at this stage.
                let t = cache_i.ellipsoid_eval[surf_id].unwrap().powi(2)
                    / Into::<T>::into(t.shift.t.powi(2));

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
                    let sup = t / (t + mij_sq);
                    s *= sup;
                    sup
                };

                if self.settings.general.debug > 2 {
                    println!(
                        "  | surf {}: t={:e}, suppression={:e}",
                        surf_id,
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
                    kappa_source[ii] += -dir / normalization * s;
                }
            } else {
                for ii in 0..self.n_loops {
                    let dir = loop_momenta[ii] - d.deformation_sources[ii].cast();
                    // the kappa returned by this function is expected to be dimensionless
                    kappa_source[ii] += -dir * s
                        / Hyperdual::<T, N>::from_real(T::convert_from(&self.e_cm_squared.sqrt()));
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
                    1. / overlap_structure.len() as f64,
                ));
            }
        }

        let mut lambda_sq = Hyperdual::<T, N>::one();

        // Setup a veto region at the intersection of pinched and non-pinched E-surfaces.
        if self.settings.deformation.fixed.ir_handling_strategy != IRHandling::None {
            let mut min_ellipse = Hyperdual::<T, N>::from_real(Into::<T>::into(1e99));
            let mut min_pinch = Hyperdual::<T, N>::from_real(Into::<T>::into(1e99));
            for s in &self.thresholds {
                if s.pinched && self.settings.deformation.fixed.ir_beta_pinch < 0. {
                    continue;
                } else if !s.pinched && self.settings.deformation.fixed.ir_beta_ellipse < 0. {
                    continue;
                }

                // TODO: the unwrap_or below is necessary for when a deformation is active but
                // exactly zero because no E-surfaces exist! This case must be better handled!
                let e = cache_i.ellipsoid_eval[s.id]
                    .unwrap_or(Hyperdual::<T, N>::from_real(Into::<T>::into(1e99)))
                    .powi(2);
                let n = Into::<T>::into(s.shift.t.powi(2));

                let t = (e
                    / (Into::<T>::into(self.settings.deformation.fixed.ir_k_com)
                        * T::convert_from(&self.e_cm_squared)
                        + Into::<T>::into(self.settings.deformation.fixed.ir_k_shift) * n))
                    .pow(Into::<T>::into(self.settings.deformation.fixed.ir_alpha));

                if !s.pinched && min_ellipse > t {
                    min_ellipse = t;
                } else if s.pinched && min_pinch > t {
                    min_pinch = t;
                }
            }

            let mut min_e = Hyperdual::<T, N>::from_real(Into::<T>::into(1e99));
            if self.settings.deformation.fixed.ir_beta_energy >= 0. {
                for f in &self.propagators {
                    let normalised_e = cache_i.cut_info[f.id].spatial_and_mass_sq
                        / T::convert_from(&self.e_cm_squared);
                    if min_e > normalised_e {
                        min_e = normalised_e;
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
                            lambda_sq = Hyperdual::<T, N>::from_real(T::nan()); // TODO: improve into a proper escalation
                        }
                        IRHandling::DismissDeformation => {
                            lambda_sq = Hyperdual::<T, N>::zero();
                        }
                        IRHandling::None => {
                            unreachable!();
                        }
                    }
                }
            }
        }

        let lambda = lambda_sq.sqrt();

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

        if self.settings.general.debug > 2 {
            for ii in 0..self.n_loops {
                println!("  | kappa{}={:e}", ii + 1, kappas[ii].real());
            }
        }

        cache.overlap_structure = overlap_structure;

        kappas
    }

    fn set_deformation_params<T: FloatLike, const N: usize>(
        &self,
        loop_momenta: &[LorentzVector<Hyperdual<T, N>>],
        cache: &mut DeformationCache<T>,
    ) where
        DeformationCache<T>: CacheSelector<T, N>,
    {
        // compute all cut energies
        let cache_i = cache.get_cache_mut();

        // reset surface evaluations
        for e in &mut cache_i.ellipsoid_eval {
            *e = None;
        }

        // update the cut info
        for p in &self.propagators {
            let mut mom = LorentzVector::default();
            for (l, &c) in izip!(loop_momenta.iter(), p.amp_sig.iter()) {
                mom += l.multiply_sign(c);
            }

            let q: LorentzVector<Hyperdual<T, N>> = p.shift.cast();
            let m_sq =
                Hyperdual::<T, N>::from_real(Into::<T>::into(p.mass) * Into::<T>::into(p.mass));

            let mom_sq = (mom + q).spatial_squared_impr();
            let cm = mom_sq + m_sq;
            let energy = cm.sqrt();
            cache_i.cut_energies[p.id] = energy;

            let info = &mut cache_i.cut_info[p.id];
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

    fn deform_generic<T: FloatLike, const N: usize>(
        &self,
        loop_momenta: &[LorentzVector<Hyperdual<T, N>>],
        cache: &mut DeformationCache<T>,
    ) -> ([LorentzVector<Hyperdual<T, N>>; MAX_LOOP], Hyperdual<T, N>)
    where
        DeformationCache<T>: CacheSelector<T, N>,
    {
        self.set_deformation_params(loop_momenta, cache);

        let mut kappas = self.deform_fixed(loop_momenta, cache);

        // make sure the kappa has the right dimension by multiplying in the scale
        let scale = Hyperdual::<T, N>::from_real(
            T::convert_from(&self.e_cm_squared).sqrt()
                * T::from_f64(self.settings.deformation.overall_scaling_constant).unwrap(),
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
            Hyperdual::from_real(Into::<T>::into(self.settings.deformation.scaling.lambda).abs())
        };

        cache.overall_lambda = lambda.real();

        if self.settings.general.debug > 2 {
            println!("lambda={:e}", lambda.real());
        }

        (kappas, lambda)
    }

    /// Compute the deformation and the partial derivatives in all the loop momenta
    /// and write them linearized in the output
    pub fn deform<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<T>],
        cache: &mut DeformationCache<T>,
        mut output: (&mut [T], &mut [T]),
    ) {
        #[cfg(not(feature = "higher_loops"))]
        deform!(self, loop_momenta, cache, deform_generic, output, 1, 2, 3);

        #[cfg(feature = "higher_loops")]
        deform!(
            self,
            loop_momenta,
            cache,
            deform_generic,
            output,
            1,
            2,
            3,
            4,
            5,
            6
        );
    }
}
