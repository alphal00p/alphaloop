use amplitude::Amplitude;
use dual_num::{DualN, Scalar, U10, U13, U16, U19, U4, U7};
use float;
use num::Complex;
use num_traits::{Signed, Zero};
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
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
    pub overall_lambda: T, // used to log the minimum
    pub numerator_momentum_cache: Vec<Complex<T>>,
    pub reduced_coefficient_lb: Vec<Complex<T>>,
    pub reduced_coefficient_cb: Vec<Complex<T>>,
    pub propagators: HashMap<(usize, usize), Complex<T>>, // TODO: remove hashmap
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
            overall_lambda: T::zero(),
            numerator_momentum_cache: vec![],
            reduced_coefficient_lb: vec![Complex::default(); topo.n_loops],
            reduced_coefficient_cb: vec![Complex::default(); topo.n_loops],
            propagators: HashMap::new(),
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
    #[serde(skip_deserializing)]
    pub ellipsoids_not_in_cuts: Vec<Vec<Vec<usize>>>,
    #[serde(default)]
    pub fixed_deformation: Vec<FixedDeformationLimit>,
    pub constant_deformation: Option<ConstantDeformation>,
    #[serde(skip_deserializing)]
    pub all_excluded_surfaces: Vec<bool>,
    #[serde(skip_deserializing)]
    pub numerator: LTDNumerator,
    #[serde(default)]
    pub numerator_tensor_coefficients: Vec<(f64, f64)>,
    #[serde(default)]
    pub amplitude: Amplitude,
    #[serde(default)]
    pub loop_momentum_map: Vec<(Vec<i8>, Vec<i8>)>,
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

    pub fn print_info(&self) {
        let num_cuts: usize = self.ltd_cut_options.iter().map(|c| c.len()).sum();
        println!("Number of cuts: {}", num_cuts);
        let mut n_unique_e_surface = 0;
        for (surf_index, surf) in self.surfaces.iter().enumerate() {
            if surf_index != surf.group
                || surf.surface_type != SurfaceType::Ellipsoid
                || !surf.exists
            {
                continue;
            }
            n_unique_e_surface += 1;
        }
        println!(
            "Number of unique existing non-pinched E-surfaces: {}",
            n_unique_e_surface
        );
        if self.fixed_deformation.len() > 0 {
            let maximal_overlap_structure: Vec<i32> = self.fixed_deformation[0]
                .deformation_per_overlap
                .iter()
                .map(|c| (n_unique_e_surface - (c.excluded_surface_ids.len() as i32)))
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
        println!(
            "Expansion threshold considered: {}",
            self.get_expansion_threshold()
        );
        println!(
            "M_ij considered: {}",
            self.settings.deformation.fixed.m_ij.abs() * self.compute_min_mij()
        );
        match self.analytical_result_real {
            Some(_) => println!(
                "Analytic result: {:e}",
                num::Complex::<f64>::new(
                    self.analytical_result_real.unwrap(),
                    self.analytical_result_imag.unwrap()
                )
            ),
            _ => println!("Analytic result not available."),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct LTDNumerator {
    pub coefficients: Vec<Complex<f64>>,
    pub n_loops: usize,
    pub max_rank: usize,
    pub reduced_size: usize,
    pub sorted_linear: Vec<Vec<usize>>,
    pub coefficient_index_map: Vec<(usize, usize)>,
    pub coefficient_index_to_powers: Vec<[usize; MAX_LOOP]>,
    pub reduced_coefficient_index_to_powers: Vec<[usize; MAX_LOOP]>,
    pub powers_to_position: HashMap<[usize; MAX_LOOP], usize>,
}
