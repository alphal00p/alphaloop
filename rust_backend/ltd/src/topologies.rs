use amplitude::Amplitude;
use dashboard::{StatusUpdate, StatusUpdateSender};
use dual_num::{DualN, Scalar, U10, U13, U16, U19, U4, U7};
use float;
use fnv::FnvHashMap;
use num::Complex;
use num_traits::{Float, Signed, Zero};
use partial_fractioning::PFCache;
use partial_fractioning::PartialFractioning;
use scs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::mem;
use std::ptr;
use std::slice;
use std::time::Instant;
use utils;
use utils::Signum;
use vector::{LorentzVector, RealNumberLike};
use {FloatLike, Settings, MAX_LOOP};

#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceType {
    Ellipsoid,
    Hyperboloid,
    Pinch,
}

/// Ellipsoid and hyperboloid surfaces
#[derive(Debug, Clone)]
pub struct Surface {
    pub unique_ellipsoid_id: Option<usize>,
    pub group: usize,
    pub exists: bool,
    pub surface_type: SurfaceType,
    pub cut_structure_index: usize,
    pub cut_option_index: usize,
    pub cut: Vec<Cut>,
    pub onshell_ll_index: usize,
    pub onshell_prop_index: usize,
    pub delta_sign: i8,
    pub sig_ll_in_cb: Vec<i8>,
    pub signs: Vec<i8>,
    pub shift: LorentzVector<float>,
    pub id: Vec<((usize, usize), i8, i8)>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Propagators {
    #[serde(skip_deserializing)]
    pub id: usize, // global id
    pub m_squared: f64,
    pub q: LorentzVector<f64>,
    #[serde(default)]
    pub parametric_shift: (Vec<i8>, Vec<i8>),
    #[serde(default)]
    pub signature: Vec<i8>,
    #[serde(default = "set_one")]
    pub power: usize,
}

fn set_one() -> usize {
    1
}

#[derive(Debug, Clone, Deserialize)]
pub struct LoopLine {
    pub start_node: usize,
    pub end_node: usize,
    pub signature: Vec<i8>,
    pub propagators: Vec<Propagators>,
}

pub struct CutList<'a>(pub &'a Vec<Cut>);

impl<'a> fmt::Display for CutList<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}", self.0[0])?;
        for x in self.0.iter().skip(1) {
            write!(f, ",{}", x)?;
        }
        write!(f, ")")
    }
}

#[derive(Debug, Copy, Clone, Deserialize, PartialEq)]
pub enum Cut {
    NoCut,
    PositiveCut(usize),
    NegativeCut(usize),
}

impl fmt::Display for Cut {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Cut::NoCut => write!(f, "_"),
            Cut::PositiveCut(i) => write!(f, "+{}", i),
            Cut::NegativeCut(i) => write!(f, "-{}", i),
        }
    }
}

#[derive(Debug, Clone)]
/// A cache for objects needed during LTD computation
pub struct CutInfo<T: Scalar + Signed + RealNumberLike, U: dual_num::Dim + dual_num::DimName>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    pub id: usize,
    pub momentum: LorentzVector<DualN<T, U>>,
    pub real_energy: DualN<T, U>,
    pub spatial_and_mass_sq: DualN<T, U>,
    pub shift: LorentzVector<DualN<T, U>>,
    pub kappa: LorentzVector<DualN<T, U>>,
    pub kappa_sq: DualN<T, U>,
    pub kappa_dot_mom: DualN<T, U>,
    pub mass: DualN<T, U>,
    pub a: DualN<T, U>,
    pub b: DualN<T, U>,
    pub c: DualN<T, U>,
}

impl<T: Scalar + Signed + RealNumberLike, U: dual_num::Dim + dual_num::DimName> Default
    for CutInfo<T, U>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    fn default() -> CutInfo<T, U> {
        CutInfo {
            id: 0,
            momentum: LorentzVector::default(),
            real_energy: DualN::default(),
            spatial_and_mass_sq: DualN::default(),
            shift: LorentzVector::default(),
            kappa: LorentzVector::default(),
            kappa_sq: DualN::default(),
            kappa_dot_mom: DualN::default(),
            mass: DualN::default(),
            a: DualN::default(),
            b: DualN::default(),
            c: DualN::default(),
        }
    }
}

#[derive(Debug, Clone)]
/// A cache for objects needed during LTD computation
pub struct LTDCacheI<T: Scalar + Signed + RealNumberLike, U: dual_num::Dim + dual_num::DimName>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    pub ellipsoid_eval: Vec<Option<DualN<T, U>>>,
    pub deform_dirs: Vec<LorentzVector<DualN<T, U>>>,
    pub non_empty_cuts: Vec<(usize, usize)>,
    pub deformation_jacobian: Vec<Complex<T>>,
    pub cut_energies: Vec<DualN<T, U>>,
    pub cut_info: Vec<CutInfo<T, U>>,
    pub computed_cut_ll: Vec<usize>,
}

impl<T: Scalar + Signed + RealNumberLike, U: dual_num::Dim + dual_num::DimName> Default
    for LTDCacheI<T, U>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    fn default() -> LTDCacheI<T, U> {
        LTDCacheI {
            ellipsoid_eval: vec![],
            deform_dirs: vec![],
            non_empty_cuts: vec![],
            deformation_jacobian: vec![],
            cut_energies: vec![],
            cut_info: vec![],
            computed_cut_ll: vec![],
        }
    }
}

impl<T: Scalar + Signed + RealNumberLike, U: dual_num::Dim + dual_num::DimName> LTDCacheI<T, U>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    fn new(num_loops: usize, num_surfaces: usize, num_propagators: usize) -> LTDCacheI<T, U> {
        LTDCacheI {
            ellipsoid_eval: vec![None; num_surfaces],
            deform_dirs: vec![LorentzVector::default(); num_surfaces * num_loops],
            non_empty_cuts: vec![(0, 0); num_surfaces],
            deformation_jacobian: vec![Complex::default(); 9 * num_loops * num_loops],
            cut_energies: vec![DualN::default(); num_propagators],
            cut_info: vec![CutInfo::default(); num_propagators],
            computed_cut_ll: vec![0; num_propagators],
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct LTDCache<T: Scalar + Signed + RealNumberLike> {
    one_loop: LTDCacheI<T, U4>,
    two_loop: LTDCacheI<T, U7>,
    three_loop: LTDCacheI<T, U10>,
    four_loop: LTDCacheI<T, U13>,
    five_loop: LTDCacheI<T, U16>,
    six_loop: LTDCacheI<T, U19>,
    pub complex_cut_energies: Vec<Complex<T>>,
    pub complex_prop_spatial: Vec<Complex<T>>,
    pub complex_loop_line_eval: Vec<Vec<[Complex<T>; 2]>>,
    pub complex_ellipsoids: Vec<Vec<Complex<T>>>,
    pub overall_lambda: T, // used to log the minimum
    pub numerator_momentum_cache: Vec<Complex<T>>,
    pub reduced_coefficient_lb: Vec<Vec<Complex<T>>>,
    pub reduced_coefficient_lb_supergraph: Vec<Vec<Complex<T>>>,
    pub reduced_coefficient_cb: Vec<Complex<T>>,
    pub pf_cache: PFCache,
    pub propagators: FnvHashMap<(usize, usize), Complex<T>>, // TODO: remove hashmap
    pub propagators_eval: Vec<Complex<T>>,
    pub propagator_powers: Vec<usize>,
}

impl<T: Scalar + Signed + RealNumberLike> LTDCache<T> {
    pub fn new(topo: &Topology) -> LTDCache<T> {
        let num_propagators = topo.loop_lines.iter().map(|x| x.propagators.len()).sum();
        LTDCache {
            one_loop: LTDCacheI::<T, U4>::new(topo.n_loops, topo.surfaces.len(), num_propagators),
            two_loop: LTDCacheI::<T, U7>::new(topo.n_loops, topo.surfaces.len(), num_propagators),
            three_loop: LTDCacheI::<T, U10>::new(
                topo.n_loops,
                topo.surfaces.len(),
                num_propagators,
            ),
            four_loop: LTDCacheI::<T, U13>::new(topo.n_loops, topo.surfaces.len(), num_propagators),
            five_loop: LTDCacheI::<T, U16>::new(topo.n_loops, topo.surfaces.len(), num_propagators),
            six_loop: LTDCacheI::<T, U19>::new(topo.n_loops, topo.surfaces.len(), num_propagators),
            complex_cut_energies: vec![Complex::default(); num_propagators],
            complex_prop_spatial: vec![Complex::default(); num_propagators],
            complex_loop_line_eval: topo
                .loop_lines
                .iter()
                .map(|ll| vec![[Complex::default(); 2]; ll.propagators.len()])
                .collect(),
            complex_ellipsoids: vec![vec![Complex::default(); num_propagators]; num_propagators],
            overall_lambda: T::zero(),
            numerator_momentum_cache: vec![],
            reduced_coefficient_lb: vec![
                vec![Complex::default(); topo.n_loops];
                if topo.settings.general.use_amplitude {
                    topo.amplitude.diagrams.len()
                } else {
                    1
                }
            ],
            reduced_coefficient_lb_supergraph: vec![vec![Complex::default(); topo.n_loops]],
            reduced_coefficient_cb: vec![Complex::default(); topo.n_loops],
            pf_cache: PFCache::new(num_propagators), // NOTE: 1-Loop
            propagators: HashMap::default(),
            propagators_eval: vec![Complex::zero(); num_propagators],
            propagator_powers: vec![1; num_propagators],
        }
    }
}

pub trait CacheSelector<T: Scalar + Signed + RealNumberLike, U: dual_num::Dim + dual_num::DimName>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    fn get_cache(&self) -> &LTDCacheI<T, U>
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy;

    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U>
    where
        dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
        dual_num::Owned<T, U>: Copy;
}

impl<T: Scalar + Signed + RealNumberLike> CacheSelector<T, U4> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U4> {
        &self.one_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U4> {
        &mut self.one_loop
    }
}

impl<T: Scalar + Signed + RealNumberLike> CacheSelector<T, U7> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U7> {
        &self.two_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U7> {
        &mut self.two_loop
    }
}

impl<T: Scalar + Signed + RealNumberLike> CacheSelector<T, U10> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U10> {
        &self.three_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U10> {
        &mut self.three_loop
    }
}

impl<T: Scalar + Signed + RealNumberLike> CacheSelector<T, U13> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U13> {
        &self.four_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U13> {
        &mut self.four_loop
    }
}

impl<T: Scalar + Signed + RealNumberLike> CacheSelector<T, U16> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U16> {
        &self.five_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U16> {
        &mut self.five_loop
    }
}

impl<T: Scalar + Signed + RealNumberLike> CacheSelector<T, U19> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U19> {
        &self.six_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U19> {
        &mut self.six_loop
    }
}

#[derive(Default)]
pub struct LTDCacheAllPrecisions {
    float_cache: LTDCache<float>,
    quad_cache: LTDCache<f128::f128>,
}

impl LTDCacheAllPrecisions {
    pub fn new(topo: &Topology) -> LTDCacheAllPrecisions {
        LTDCacheAllPrecisions {
            float_cache: LTDCache::new(topo),
            quad_cache: LTDCache::new(topo),
        }
    }
}

pub trait CachePrecisionSelector<T: FloatLike> {
    fn get(&mut self) -> &mut LTDCache<T>;
}

impl CachePrecisionSelector<float> for LTDCacheAllPrecisions {
    #[inline]
    fn get(&mut self) -> &mut LTDCache<float> {
        &mut self.float_cache
    }
}

impl CachePrecisionSelector<f128::f128> for LTDCacheAllPrecisions {
    #[inline]
    fn get(&mut self) -> &mut LTDCache<f128::f128> {
        &mut self.quad_cache
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct FixedDeformationOverlap {
    pub deformation_sources: Vec<LorentzVector<f64>>,
    pub excluded_surface_ids: Vec<Vec<((usize, usize), i8, i8)>>,
    pub overlap: Option<Vec<usize>>,
    #[serde(default)]
    pub radius: f64,
    #[serde(skip_deserializing)]
    pub excluded_surface_indices: Vec<usize>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct FixedDeformationLimit {
    pub deformation_per_overlap: Vec<FixedDeformationOverlap>,
    pub excluded_propagators: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct ConstantDeformation {
    pub alpha: Vec<LorentzVector<f64>>,
    pub beta: Vec<LorentzVector<f64>>,
}

#[derive(Debug, Clone, Default)]
pub struct LTDNumerator {
    pub coefficients: Vec<Complex<f64>>,
    pub n_loops: usize,
    pub max_rank: usize,
    pub reduced_size: usize,
    pub sorted_linear: Vec<Vec<usize>>,
    pub coefficient_index_map: Vec<(usize, usize)>,
    pub coefficient_index_to_powers: Vec<[u8; MAX_LOOP]>,
    pub reduced_coefficient_index_to_powers: Vec<[u8; MAX_LOOP]>,
    pub reduced_blocks: Vec<usize>,
    pub coefficients_modified: bool,
    pub non_empty_coeff_map_to_reduced_numerator: Vec<Vec<(usize, usize, Complex<f64>)>>,
    pub coeff_map_to_reduced_numerator: Vec<Vec<usize>>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Topology {
    pub name: String,
    pub n_loops: usize,
    pub analytical_result_real: Option<f64>,
    pub analytical_result_imag: Option<f64>,
    pub maximum_ratio_expansion_threshold: f64,
    #[serde(default)]
    pub global_seed: Option<[u8; 32]>,
    #[serde(default)]
    pub e_cm_squared: f64,
    #[serde(default)]
    pub on_shell_flag: usize,
    pub external_kinematics: Vec<LorentzVector<f64>>,
    pub loop_lines: Vec<LoopLine>,
    #[serde(default)]
    pub propagator_id_to_ll_id: Vec<(usize, usize)>,
    pub ltd_cut_structure: Vec<Vec<i8>>,
    #[serde(default)]
    pub ltd_cut_options: Vec<Vec<Vec<Cut>>>, // cartesian product of cut structures
    #[serde(default)]
    pub cb_to_lmb_mat: Vec<Vec<i8>>, // a map from cut momenta to topology loop momenta
    #[serde(default)]
    pub settings: Settings,
    #[serde(default, skip_deserializing)]
    pub surfaces: Vec<Surface>,
    #[serde(default, skip_deserializing)]
    pub all_ellipsoid_surfaces: Vec<Surface>,
    #[serde(skip_deserializing)]
    pub rotation_matrix: [[float; 3]; 3],
    #[serde(default)]
    pub fixed_deformation: Vec<FixedDeformationLimit>,
    pub constant_deformation: Option<ConstantDeformation>,
    #[serde(skip_deserializing)]
    pub all_excluded_surfaces: Vec<bool>,
    #[serde(skip_deserializing)]
    pub numerator: LTDNumerator,
    #[serde(skip_deserializing)]
    pub partial_fractioning: PartialFractioning,
    #[serde(default)]
    pub numerator_tensor_coefficients: Vec<(f64, f64)>,
    #[serde(default)]
    pub amplitude: Amplitude,
    #[serde(default)]
    pub loop_momentum_map: Vec<(Vec<i8>, Vec<i8>)>,
    #[serde(skip_deserializing)]
    pub socp_problem: SOCPProblem,
}

impl Topology {
    pub fn from_file(filename: &str, settings: &Settings) -> HashMap<String, Topology> {
        let f = File::open(filename).expect("Could not open topology file");

        let mut topologies: Vec<Topology> = serde_yaml::from_reader(f).unwrap();

        for t in &mut topologies {
            t.settings = settings.clone();
        }

        topologies
            .into_iter()
            .map(|t| (t.name.clone(), t))
            .collect()
    }

    pub fn print_info(&self, status_update_sender: &mut StatusUpdateSender) {
        let num_cuts: usize = self.ltd_cut_options.iter().map(|c| c.len()).sum();
        if self.settings.general.debug > 0 {
            println!("Number of cuts: {}", num_cuts);
        }
        let mut n_unique_e_surface = 0;
        let mut n_pinches = 0;
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            if surf_index == surf.group && surf.exists {
                if surf.surface_type == SurfaceType::Ellipsoid {
                    n_unique_e_surface += 1;
                } else if surf.surface_type == SurfaceType::Pinch {
                    n_pinches += 1;
                }
            }
        }
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Number of unique existing E-surfaces: {}",
                n_unique_e_surface
            )))
            .unwrap();
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Number of unique existing pinches: {}",
                n_pinches
            )))
            .unwrap();
        if self.fixed_deformation.len() > 0 {
            let maximal_overlap_structure: Vec<i32> = self
                .fixed_deformation
                .iter()
                .find(|x| x.excluded_propagators.len() == 0)
                .unwrap()
                .deformation_per_overlap
                .iter()
                .map(|c| (n_unique_e_surface - (c.excluded_surface_indices.len() as i32)))
                .collect();
            let radii: Vec<f64> = self
                .fixed_deformation
                .iter()
                .flat_map(|fd| fd.deformation_per_overlap.iter().map(|fdo| fdo.radius))
                .collect();
            let n_sources: usize = self
                .fixed_deformation
                .iter()
                .flat_map(|fd| fd.deformation_per_overlap.iter().map(|_| 1))
                .sum();

            if self.settings.general.debug > 0 {
                println!(
                    "Number of E-surfaces part of each maximal overlap: {:?}",
                    maximal_overlap_structure
                );
                println!("Total number of sources: {}", n_sources);
                println!(
                    "Min radius: {}",
                    radii
                        .iter()
                        .fold(std::f64::INFINITY, |acc, x| f64::min(acc, *x))
                );
                println!(
                    "Max radius: {}",
                    radii
                        .iter()
                        .fold(std::f64::NEG_INFINITY, |acc, x| f64::max(acc, *x))
                );
            }
        }

        if self.settings.general.debug > 0 {
            println!(
                "M_ij considered: {}",
                self.settings.deformation.fixed.m_ij.abs() * self.compute_min_mij()
            );
        }
        match self.analytical_result_real {
            Some(_) => status_update_sender
                .send(StatusUpdate::Message(format!(
                    "Analytic result: {:e}",
                    num::Complex::<f64>::new(
                        self.analytical_result_real.unwrap(),
                        self.analytical_result_imag.unwrap()
                    )
                )))
                .unwrap(),
            _ => {}
        }
    }

    pub fn determine_ellipsoid_overlap_structure(
        &mut self,
        update_excluded_surfaces: bool,
    ) -> Vec<FixedDeformationLimit> {
        let ellipsoid_ids: Vec<usize> = self
            .surfaces
            .iter()
            .enumerate()
            .filter_map(|(i, s)| {
                if s.exists && s.surface_type == SurfaceType::Ellipsoid && s.group == i {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        if ellipsoid_ids.is_empty() {
            if self.settings.general.debug > 1 {
                println!("No ellipsoids");
            }
            return vec![];
        }

        let overlap = self.find_overlap_structure(&ellipsoid_ids, &[]);
        if self.settings.general.debug > 1 {
            println!("overlap={:#?}", overlap);
        }

        // TODO: also determine subspace deformation vectors
        let fixed_deformation = vec![FixedDeformationLimit {
            deformation_per_overlap: overlap,
            excluded_propagators: vec![],
        }];

        if update_excluded_surfaces {
            // update the excluded surfaces
            for ex in &mut self.all_excluded_surfaces {
                *ex = false;
            }

            for d_lim in &fixed_deformation {
                for d in &d_lim.deformation_per_overlap {
                    for surf_id in &d.excluded_surface_indices {
                        self.all_excluded_surfaces[*surf_id] = true;
                    }
                }
            }
        }

        fixed_deformation
    }

    /// Construct the SCS problem with all the ellipsoids in the problem.
    fn construct_socp_problem(&mut self, ellipsoid_ids: &[usize], minimize: bool) {
        let mut p = &mut self.socp_problem;

        let multiplier = if minimize { 6 * self.n_loops } else { 1 };
        p.radius_computation = minimize;

        p.cone.l = (ellipsoid_ids.len() * multiplier) as scs::scs_int;

        let mut b_index = 0;

        for m in 0..multiplier {
            for &e_id in ellipsoid_ids {
                p.b[b_index] = -self.surfaces[e_id]
                    .shift
                    .t
                    .multiply_sign(self.surfaces[e_id].delta_sign);
                b_index += 1;
            }
        }

        let mut lm_tag = [false; MAX_LOOP];

        // construct the focus part of the shift
        // TODO: do not generate the focus more than once if k[m/3] does not appear in the focus
        let mut focus_count = 0;
        for m in 0..multiplier {
            for (ll_index, ll) in self.loop_lines.iter().enumerate() {
                for (prop_index, prop) in ll.propagators.iter().enumerate() {
                    let mut included = false;
                    'el_loop: for &e_id in ellipsoid_ids {
                        for ((surf_ll_index, surf_prop_index), _, _) in &self.surfaces[e_id].id {
                            if ll_index == *surf_ll_index && prop_index == *surf_prop_index {
                                for (tag, &sig) in lm_tag.iter_mut().zip(&ll.signature) {
                                    *tag |= sig != 0;
                                }
                                included = true;
                                break 'el_loop;
                            }
                        }
                    }

                    if included {
                        p.focus_list[focus_count] = (ll_index, prop_index, m);
                        p.b[b_index] = 0.;
                        b_index += 1;
                        if prop.m_squared != 0. {
                            p.b[b_index] = prop.m_squared.sqrt();
                            b_index += 1;
                            p.q[focus_count] = 5;
                        } else {
                            p.q[focus_count] = 4;
                        };
                        p.b[b_index] = prop.q.x;
                        b_index += 1;
                        p.b[b_index] = prop.q.y;
                        b_index += 1;
                        p.b[b_index] = prop.q.z;
                        b_index += 1;
                        focus_count += 1;
                    }
                }
            }
        }

        let mut var_count = 0;
        for (vm, &x) in p.var_map.iter_mut().zip(lm_tag.iter()) {
            if x {
                *vm = var_count;
                var_count += 1;
            }
        }

        let width = if minimize {
            3 * var_count + focus_count + 1
        } else {
            3 * var_count + focus_count
        };

        p.cone.qsize = focus_count as scs::scs_int;
        p.a.m = b_index as scs::scs_int;
        p.a.n = width as scs::scs_int;
        p.data.m = p.a.m;
        p.data.n = p.a.n;

        for aa in &mut p.a_dense[..width * b_index] {
            *aa = 0.;
        }

        for cc in &mut p.c[..width] {
            *cc = 0.;
        }

        let focus_start = if minimize {
            p.c[3 * var_count] += -1.; // maximize the radius, TODO: sign?
            3 * var_count + 1
        } else {
            3 * var_count
        };

        // write the ellipsoid constraints: f1 + f2 < c, etc
        let mut row_counter = 0; // the start of the focus information
        for m in 0..multiplier {
            for &e_id in ellipsoid_ids {
                for ((ll_index, prop_index), _, _) in self.surfaces[e_id].id.iter() {
                    let foc_index = p
                        .focus_list
                        .iter()
                        .position(|x| *x == (*ll_index, *prop_index, m))
                        .unwrap();
                    p.a_dense[row_counter * width + focus_start + foc_index] = 1.;
                }
                row_counter += 1;
            }
        }

        // now for every focus f, add the norm constraints |k_x+p_x,m| < f
        for (foc_index, (ll_index, prop_index, radius_index)) in
            p.focus_list[..focus_count].iter().enumerate()
        {
            let prop = &self.loop_lines[*ll_index].propagators[*prop_index];

            // set the linear equation
            p.a_dense[row_counter * width + focus_start + foc_index] = -1.0;
            row_counter += 1;
            if prop.m_squared != 0. {
                // the mass is an empty line
                row_counter += 1;
            }

            let (rad_lm, rad_dir, rad_sign) =
                (radius_index / 6, radius_index % 3, radius_index % 2);

            for dir_index in 0..3 {
                for (var_index, &v) in prop.signature.iter().enumerate() {
                    if v != 0 {
                        p.a_dense[row_counter * width + 3 * p.var_map[var_index] + dir_index] =
                            Into::<float>::into(-v);

                        if minimize && rad_lm == var_index && rad_dir == dir_index {
                            p.a_dense[row_counter * width + 3 * var_count] =
                                if rad_sign == 0 { -1. } else { 1.0 };
                        }
                    }
                }
                row_counter += 1;
            }
        }

        // now write a in compressed column format
        let mut non_zero_index = 0;
        p.p[0] = 0;
        for col in 0..width {
            for row in 0..b_index {
                if p.a_dense[row * width + col] != 0. {
                    p.x[non_zero_index] = p.a_dense[row * width + col];
                    p.i[non_zero_index] = row as scs::scs_int;
                    non_zero_index += 1;
                }
            }
            p.p[col + 1] = non_zero_index as scs::scs_int;
        }
    }

    /// Evaluate an E-surface, where eval < 0 is the inside.
    fn evaluate_surface<T: FloatLike>(&self, loop_momenta: &[LorentzVector<T>], s: &Surface) -> T {
        let mut eval = Into::<T>::into(s.shift.t).multiply_sign(s.delta_sign);
        for ((ll, pp), _, _) in &s.id {
            let prop = &self.loop_lines[*ll].propagators[*pp];
            // get the loop momentum
            let mut mom = prop.q.cast();
            for (m, sig) in loop_momenta[..self.n_loops].iter().zip(&prop.signature) {
                mom += m.multiply_sign(*sig);
            }
            eval += (mom.spatial_squared() + Into::<T>::into(prop.m_squared)).sqrt();
        }
        eval
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
                for appears in &self.socp_problem.pair_appears_in_overlap_structurex
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
                        if !self.socp_problem.pair_overlap_matrix
                            [*x as usize * num_ellipsoids + new_num]
                        {
                            continue 'next_num_loop;
                        }

                        if !self.socp_problem.pair_appears_in_overlap_structurex
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
                                if self.socp_problem.pair_overlap_matrix[k * num_ellipsoids + j]
                                    && !self.socp_problem.pair_appears_in_overlap_structurex
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
        all_ellipsoids: &[usize],
    ) -> FixedDeformationOverlap {
        // TODO: recycle from pre-allocated overlap structures
        let mut var_count = 0;
        let mut sources = vec![LorentzVector::default(); self.n_loops];
        for &i in &self.socp_problem.var_map[..self.n_loops] {
            if i < self.n_loops {
                sources[i] = LorentzVector::from_args(
                    0.,
                    self.socp_problem.sol_x[i * 3],
                    self.socp_problem.sol_x[i * 3 + 1],
                    self.socp_problem.sol_x[i * 3 + 2],
                );
                var_count += 1;
            }
        }

        let mut excluded_surface_indices = Vec::with_capacity(all_ellipsoids.len());
        for e in all_ellipsoids {
            if !ellipsoid_list.contains(e) {
                excluded_surface_indices.push(*e);
            }
        }

        let mut radius = 0.;

        if self.socp_problem.radius_computation {
            radius = self.socp_problem.sol_x[var_count * 3];
        } else {
            // give a meausure for distance
            for (foc_count, foc_val) in self.socp_problem.c
                [3 * var_count..self.socp_problem.a.n as usize]
                .iter()
                .zip(&self.socp_problem.sol_x[3 * var_count..self.socp_problem.a.n as usize])
            {
                radius += foc_count * foc_val;
            }
        }

        FixedDeformationOverlap {
            deformation_sources: sources,
            excluded_surface_ids: vec![],
            excluded_surface_indices: excluded_surface_indices,
            overlap: Some(ellipsoid_list.to_vec()),
            radius,
        }
    }

    /// Find a point on a surface along a line and return the scale:
    /// `E(dir * t + k) = 0`, returning `t` closest to `t_start`.
    fn get_scale_to_reach_surface<T: FloatLike>(
        &self,
        s: &Surface,
        k: &[LorentzVector<T>],
        dir: &[LorentzVector<T>],
        t_start: T,
    ) -> Option<T> {
        let mut t = t_start;
        let tolerance = Into::<T>::into(1e-15 * self.e_cm_squared.sqrt());
        for it in 0..30 {
            let mut f = Into::<T>::into(s.shift.t.multiply_sign(s.delta_sign));
            let mut df = T::zero();

            for ((ll_index, prop_index), _, _) in &s.id {
                let shift = self.loop_lines[*ll_index].propagators[*prop_index].q.cast();
                let k = utils::evaluate_signature(&self.loop_lines[*ll_index].signature, k);
                let dir = utils::evaluate_signature(&self.loop_lines[*ll_index].signature, dir);
                let energy = ((k + dir * t + shift).spatial_squared()
                    + Into::<T>::into(
                        self.loop_lines[*ll_index].propagators[*prop_index].m_squared,
                    ))
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

    /// Generate a point inside an E-surface. A point that is guaranteed to be on the inside is
    /// `c_i = -p s_i m_i/(sum_j m_j)`, in the cut basis where `s` is the surface signature in the cut basis.
    fn get_point_inside_surface<T: FloatLike>(&self, s: &Surface) -> [LorentzVector<T>; MAX_LOOP] {
        let mut ll_on_foci = [LorentzVector::default(); MAX_LOOP];
        let mut cut_shift = [LorentzVector::<T>::default(); MAX_LOOP]; // qs of cut
        let mut cut_masses = [T::zero(); MAX_LOOP];
        let mut cut_counter = 0;
        let mut mass_sum = T::zero();
        for (cut, ll) in s.cut.iter().zip(self.loop_lines.iter()) {
            if let Cut::NegativeCut(cut_prop_index) | Cut::PositiveCut(cut_prop_index) = cut {
                // only add the shift if this cut momentum occurs in the surface
                if s.sig_ll_in_cb[cut_counter] != 0 {
                    cut_shift[cut_counter] = ll.propagators[*cut_prop_index].q.cast();
                    cut_masses[cut_counter] =
                        Into::<T>::into(ll.propagators[*cut_prop_index].m_squared.sqrt());
                    mass_sum += cut_masses[cut_counter];
                }
                cut_counter += 1;
            }
        }

        mass_sum += Into::<T>::into(
            self.loop_lines[s.onshell_ll_index].propagators[s.onshell_prop_index]
                .m_squared
                .sqrt(),
        );

        // do the basis transformation
        for (i, lm) in ll_on_foci[..self.n_loops].iter_mut().enumerate() {
            for (((&sign, shift), &mass), &cut_mom_sig) in self.cb_to_lmb_mat[s.cut_structure_index]
                [i * self.n_loops..(i + 1) * self.n_loops]
                .iter()
                .zip(cut_shift[..self.n_loops].iter())
                .zip(cut_masses[..self.n_loops].iter())
                .zip(s.sig_ll_in_cb.iter())
            {
                if sign != 0 {
                    *lm += -shift.multiply_sign(sign);

                    if mass_sum > T::zero() {
                        *lm += -s
                            .shift
                            .cast()
                            .multiply_sign(cut_mom_sig)
                            .multiply_sign(sign)
                            * (mass / mass_sum);
                    }
                }
            }
            lm.t = T::zero();
        }

        debug_assert!(self.evaluate_surface(&ll_on_foci, s) < T::zero());
        ll_on_foci
    }

    /// Get the normal at the point on an E-surface. We always assume that the surface is positive.
    fn get_normal<T: FloatLike>(
        &self,
        surface: &Surface,
        k: &[LorentzVector<T>],
        normalize: bool,
    ) -> [LorentzVector<T>; MAX_LOOP] {
        let mut df = [LorentzVector::default(); MAX_LOOP];

        for ((ll_index, prop_index), _, _) in &surface.id {
            let shift = self.loop_lines[*ll_index].propagators[*prop_index].q.cast();
            let mut mom_part =
                utils::evaluate_signature(&self.loop_lines[*ll_index].signature, k) + shift;
            mom_part.t = T::zero();
            let energy = (mom_part.spatial_squared()
                + Into::<T>::into(self.loop_lines[*ll_index].propagators[*prop_index].m_squared))
            .sqrt();

            if energy == T::zero() || Float::is_nan(energy) {
                continue;
            }

            for (d, s) in df.iter_mut().zip(&self.loop_lines[*ll_index].signature) {
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
    fn overlap_heuristic(&self, surf1: &Surface, surf2: &Surface) -> bool {
        let mut x1 = self.get_point_inside_surface::<f64>(surf1);
        let mut x2 = self.get_point_inside_surface::<f64>(surf2);

        if self.evaluate_surface(&x1, surf2) < 0. || self.evaluate_surface(&x2, surf1) < 0. {
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
            debug_assert!(self.evaluate_surface(&x1, surf1).abs() < 1e-14 * self.e_cm_squared);
            debug_assert!(self.evaluate_surface(&x2, surf2).abs() < 1e-14 * self.e_cm_squared);

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
                    -self.e_cm_squared, // start with an underestimate
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

            debug_assert!(self.evaluate_surface(&x1, surf1) < 0.);
            debug_assert!(self.evaluate_surface(&x2, surf2) < 0.);

            if it == 49 {
                if self.settings.general.debug > 2 {
                    println!(
                        "  | no convergence after {} iterations for surface {} and {}",
                        it, surf1.group, surf2.group
                    );
                }
            }
        }
        false
    }

    fn find_overlap_structure(
        &mut self,
        ellipsoid_list: &[usize],
        _extra_constraints: &[usize],
    ) -> Vec<FixedDeformationOverlap> {
        let mut overlap_structure = vec![];

        let t = Instant::now();

        let mut problem_count = 0;
        let mut pair_overlap_count = 0;

        if self.settings.deformation.fixed.use_heuristic_centers {
            if ellipsoid_list.len() == 1 {
                // put ourselves on a focus
                let ll_on_foci = self.get_point_inside_surface(&self.surfaces[ellipsoid_list[0]]);
                let r: f64 = -self.evaluate_surface(
                    &ll_on_foci[..self.n_loops],
                    &self.surfaces[ellipsoid_list[0]],
                );
                assert!(r > 0.);

                return vec![FixedDeformationOverlap {
                    deformation_sources: ll_on_foci[..self.n_loops].to_vec(),
                    excluded_surface_ids: vec![],
                    excluded_surface_indices: vec![],
                    overlap: Some(ellipsoid_list.to_vec()),
                    radius: r,
                }];
            }

            // the origin is often inside all surfaces, so do a quick test for that
            let mut origin_inside_radius = self.e_cm_squared;
            let origin = [LorentzVector::default(); MAX_LOOP];
            for &e in ellipsoid_list {
                let r: f64 = -self.evaluate_surface(&origin[..self.n_loops], &self.surfaces[e]);
                if r < origin_inside_radius {
                    origin_inside_radius = r;
                }
                if r < 0. {
                    if self.settings.general.debug > 0 {
                        println!("Origin not inside for {:?}: {}", &self.surfaces[e], r);
                        let id = &self.surfaces[e].id;
                        println!(
                            "shift1={}",
                            self.loop_lines[(id[0].0).0].propagators[(id[0].0).1].q
                        );
                        println!(
                            "shift2={}",
                            self.loop_lines[(id[1].0).0].propagators[(id[1].0).1].q
                        );
                    }
                    break;
                }
            }

            if origin_inside_radius > 0. {
                return vec![FixedDeformationOverlap {
                    deformation_sources: origin[..self.n_loops].to_vec(),
                    excluded_surface_ids: vec![],
                    excluded_surface_indices: vec![],
                    overlap: Some(ellipsoid_list.to_vec()),
                    radius: origin_inside_radius,
                }];
            }
        }

        // first check for full overlap before constructing pair information
        /*self.construct_socp_problem(&ellipsoid_list, true, false);
        if self.socp_problem.solve() > 0 {
            overlap_structure
                .push(self.build_deformation_overlap_info(ellipsoid_list, ellipsoid_list));
            self.socp_problem.finish();
            return overlap_structure;
        }*/

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
                for (foc1, _, _) in self.surfaces[ellipsoid_list[i]].id.iter() {
                    for (foc2, _, _) in self.surfaces[ellipsoid_list[j]].id.iter() {
                        if foc1.0 == foc2.0 {
                            loop_line_overlap = true;
                        }
                        if foc1 == foc2 {
                            num_overlap += 1;
                        }
                    }
                }

                let mut focus_overlap = false;
                if !loop_line_overlap
                    || (num_overlap + 1 >= self.surfaces[ellipsoid_list[i]].id.len()
                        && num_overlap + 1 >= self.surfaces[ellipsoid_list[j]].id.len())
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
                let has_overlap_heuristic = self.overlap_heuristic(
                    &self.surfaces[ellipsoid_list[i]],
                    &self.surfaces[ellipsoid_list[j]],
                );
                ball_time += Instant::now().duration_since(ti).as_nanos();

                if self.settings.general.debug < 2 && has_overlap_heuristic {
                    self.socp_problem.pair_overlap_matrix[i * ellipsoid_list.len() + j] = true;
                    self.socp_problem.pair_overlap_matrix[j * ellipsoid_list.len() + i] = true;
                    pair_overlap_count += 1;
                    continue;
                }

                self.construct_socp_problem(&[ellipsoid_list[i], ellipsoid_list[j]], false);

                // perform the ECOS test
                let ti = Instant::now();
                self.socp_problem.initialize_workspace_ecos();
                let r = self.socp_problem.solve_ecos();
                let has_overlap = r == 0 || r == -2 || r == 10;
                self.socp_problem.finish_ecos();
                ecos_time += Instant::now().duration_since(ti).as_nanos();

                if focus_overlap && !has_overlap {
                    panic!(
                        "Foci overlap but SCS says that there is no overlap for {} and {}",
                        ellipsoid_list[i], ellipsoid_list[j]
                    );
                }

                if has_overlap != has_overlap_heuristic {
                    if self.settings.general.debug > 2 {
                        println!(
                            "Heuristic failed: {} vs {}",
                            has_overlap, has_overlap_heuristic
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
                            "Inconcistency between SCS and ECOS: {} vs {} for {:?}",
                            has_overlap_scs,
                            has_overlap,
                            &[ellipsoid_list[i], ellipsoid_list[j]]
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
            while self.next_set(
                ellipsoid_list.len(),
                &mut option[..n],
                &mut different[..n],
                initialize,
            ) {
                initialize = false;

                for (ot, o) in option_translated[..n].iter_mut().zip(&option[..n]) {
                    *ot = ellipsoid_list[*o as usize];
                }

                self.construct_socp_problem(&option_translated[..n], false);

                // perform the ECOS test
                let ti = Instant::now();
                self.socp_problem.initialize_workspace_ecos();
                let r = self.socp_problem.solve_ecos();
                let has_overlap = r == 0 || r == -2 || r == 10;
                self.socp_problem.finish_ecos();
                ecos_time += Instant::now().duration_since(ti).as_nanos();

                // verify with SCS
                if self.settings.general.debug > 1 {
                    self.socp_problem.initialize_workspace_scs();
                    let has_overlap_scs = self.socp_problem.solve_scs() > 0;
                    self.socp_problem.finish_scs();

                    if has_overlap != has_overlap_scs {
                        panic!(
                            "Inconcistency between SCS and ECOS: {} vs {} for {:?}",
                            has_overlap_scs,
                            has_overlap,
                            &option_translated[..n]
                        );
                    }
                }

                if has_overlap {
                    // find centers with maximal radius
                    self.construct_socp_problem(&option_translated[..n], true);
                    self.socp_problem.initialize_workspace_ecos();
                    let r = self.socp_problem.solve_ecos();
                    assert!(
                        r == 0 || r == -2 || r == 10,
                        "ECOS returns failure when it should find a center: {}",
                        r
                    );
                    self.socp_problem.finish_ecos();

                    for (i, &ei) in option[..n].iter().enumerate() {
                        for &ej in option[i + 1..n].iter() {
                            self.socp_problem.pair_appears_in_overlap_structurex
                                [ei as usize * ellipsoid_list.len() + ej as usize] = true;
                            self.socp_problem.pair_appears_in_overlap_structurex
                                [ej as usize * ellipsoid_list.len() + ei as usize] = true;
                        }
                    }

                    overlap_structure.push(
                        self.build_deformation_overlap_info(
                            &option_translated[..n],
                            ellipsoid_list,
                        ),
                    );
                }

                problem_count += 1;
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
    pub radius_computation: bool,
    pub pair_overlap_matrix: Vec<bool>,
    pub pair_appears_in_overlap_structurex: Vec<bool>,
    pub focus_list: Vec<(usize, usize, usize)>,
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
        SOCPProblem::new(50 * 6, 20 * 6)
    }
}

impl Clone for SOCPProblem {
    fn clone(&self) -> SOCPProblem {
        let mut p = SOCPProblem::default();

        // now copy the settings from the old problem
        // all other data is not copied
        *p.settings = *self.settings;
        p
    }
}

impl SOCPProblem {
    pub fn new(max_ellipsoids: usize, max_foci: usize) -> SOCPProblem {
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

        let pair_overlap_matrix = vec![false; max_ellipsoids * max_ellipsoids];
        let pair_appears_in_overlap_structurex = vec![false; max_ellipsoids * max_ellipsoids];

        let option = vec![0isize; max_ellipsoids];
        let option_translated = vec![0; max_ellipsoids];
        let different = vec![false; max_ellipsoids];

        let num_constraints = max_ellipsoids + max_foci * 5;
        let var_length = 3 * MAX_LOOP + max_foci;
        let num_non_empty = max_ellipsoids * MAX_LOOP + max_foci + 5 * MAX_LOOP * max_foci;

        let focus_list = vec![(0, 0, 0); max_foci];
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
