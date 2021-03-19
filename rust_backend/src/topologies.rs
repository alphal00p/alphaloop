use crate::amplitude::Amplitude;
use crate::dashboard::{StatusUpdate, StatusUpdateSender};
use crate::integrand::IntegrandImplementation;
use crate::integrand::IntegrandSample;
use crate::observables::EventManager;
use crate::partial_fractioning::PFCache;
use crate::partial_fractioning::{PartialFractioning, PartialFractioningMultiLoops};
use crate::utils;
use crate::{float, FloatLike, Settings, MAX_LOOP};
use color_eyre::{Help, Report};
use dual_num::{DualN, U10, U13, U16, U19, U4, U7};
use eyre::WrapErr;
use fnv::FnvHashMap;
use havana::{ContinuousGrid, Grid};
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use mpolynomial::MPolynomial;
use num::Complex;
use num_traits::{Float, FloatConst, FromPrimitive, Inv, One, Signed, ToPrimitive, Zero};
use rand::Rng;
use scs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::mem;
use std::ptr;
use std::slice;
use std::time::Instant;
use utils::Signum;

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

#[derive(Debug, Clone)]
pub struct Focus {
    pub id: usize,
    pub signature: Vec<i8>,
    pub m_squared: f64,
    pub shift: LorentzVector<f64>,
}

#[derive(Debug, Clone)]
pub struct ESurface {
    pub id: usize,
    pub foci: Vec<Focus>,
    pub shift: LorentzVector<f64>,
    pub focus_basis_signature: Vec<i8>,
    pub focus_basis_mat: Vec<f64>,
    pub point_on_the_inside: Vec<LorentzVector<f64>>,
}

impl Topology {
    pub fn print_surface(&self, s: &Surface) {
        for ((ll, pp), e_sign, _) in &s.id {
            print!("{}sqrt((", if *e_sign > 0 { "+" } else { "-" });
            let prop = &self.loop_lines[*ll].propagators[*pp];
            // get the loop momentum
            for (i, sig) in prop.signature.iter().enumerate() {
                if *sig != 0 {
                    print!("{}k{}", if *sig > 0 { "+" } else { "-" }, i);
                }
            }
            print!("+({},{},{}))^2", prop.q.x, prop.q.y, prop.q.z);

            print!("+{})", prop.m_squared);
        }
        println!("{:+} exists={}", s.shift.t, s.exists);
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct Propagators {
    #[serde(skip_deserializing)]
    pub id: usize, // global id
    pub name: String,
    pub m_squared: f64,
    pub q: LorentzVector<f64>,
    #[serde(default)]
    pub parametric_shift: (Vec<i8>, Vec<i8>),
    #[serde(default)]
    pub signature: Vec<i8>,
    #[serde(default = "set_one")]
    pub power: usize,
    #[serde(default)]
    pub uv: bool,
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
pub struct CutInfo<T: FloatLike, U: dual_num::Dim + dual_num::DimName>
where
    dual_num::DefaultAllocator: dual_num::Allocator<T, U>,
    dual_num::Owned<T, U>: Copy,
{
    pub id: usize,
    pub momentum: LorentzVector<DualN<T, U>>,
    pub real_energy: DualN<T, U>,
    pub spatial_and_mass_sq: DualN<T, U>,
    pub spatial_and_uv_mass_sq: DualN<T, U>,
    pub shift: LorentzVector<DualN<T, U>>,
    pub kappa: LorentzVector<DualN<T, U>>,
    pub kappa_sq: DualN<T, U>,
    pub kappa_dot_mom: DualN<T, U>,
    pub mass: DualN<T, U>,
    pub a: DualN<T, U>,
    pub b: DualN<T, U>,
    pub c: DualN<T, U>,
}

impl<T: FloatLike, U: dual_num::Dim + dual_num::DimName> Default for CutInfo<T, U>
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
            spatial_and_uv_mass_sq: DualN::default(),
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
pub struct LTDCacheI<T: FloatLike, U: dual_num::Dim + dual_num::DimName>
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

impl<T: FloatLike, U: dual_num::Dim + dual_num::DimName> Default for LTDCacheI<T, U>
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

impl<T: FloatLike, U: dual_num::Dim + dual_num::DimName> LTDCacheI<T, U>
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
pub struct LTDCache<T: FloatLike> {
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
    pub numerator_cache_outdated: Vec<bool>,
    pub reduced_coefficient_lb_mpoly: MPolynomial<Complex<T>>,
    pub reduced_coefficient_lb: Vec<Vec<Complex<T>>>,
    pub reduced_coefficient_lb_supergraph: Vec<Vec<Complex<T>>>,
    pub reduced_coefficient_cb: Vec<Complex<T>>,
    pub pf_cache: PFCache<T>,
    pub propagators: FnvHashMap<(usize, usize), Complex<T>>, // TODO: remove hashmap
    pub propagators_eval: Vec<Complex<T>>,
    pub propagator_powers: Vec<usize>,
    pub cached_topology_integrand: Vec<(usize, Complex<T>)>,
}

impl<T: FloatLike> LTDCache<T> {
    pub fn new(topo: &Topology) -> LTDCache<T> {
        let num_propagators = topo.loop_lines.iter().map(|x| x.propagators.len()).sum();
        let num_propagators_deg_1l = topo
            .loop_lines
            .iter()
            .filter(|x| !x.signature.iter().all(|x| *x == 0))
            .map(|x| x.propagators.iter().map(|p| p.power).sum::<usize>())
            .sum();
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
            numerator_cache_outdated: vec![],
            cached_topology_integrand: vec![],
            reduced_coefficient_lb_mpoly: MPolynomial::new(topo.n_loops),
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
            pf_cache: PFCache::<T>::new(num_propagators_deg_1l, topo.n_loops),
            propagators: HashMap::default(),
            propagators_eval: vec![Complex::zero(); num_propagators],
            propagator_powers: vec![1; num_propagators],
        }
    }
}

pub trait CacheSelector<T: FloatLike, U: dual_num::Dim + dual_num::DimName>
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

impl<T: FloatLike> CacheSelector<T, U4> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U4> {
        &self.one_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U4> {
        &mut self.one_loop
    }
}

impl<T: FloatLike> CacheSelector<T, U7> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U7> {
        &self.two_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U7> {
        &mut self.two_loop
    }
}

impl<T: FloatLike> CacheSelector<T, U10> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U10> {
        &self.three_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U10> {
        &mut self.three_loop
    }
}

impl<T: FloatLike> CacheSelector<T, U13> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U13> {
        &self.four_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U13> {
        &mut self.four_loop
    }
}

impl<T: FloatLike> CacheSelector<T, U16> for LTDCache<T> {
    #[inline]
    fn get_cache(&self) -> &LTDCacheI<T, U16> {
        &self.five_loop
    }

    #[inline]
    fn get_cache_mut(&mut self) -> &mut LTDCacheI<T, U16> {
        &mut self.five_loop
    }
}

impl<T: FloatLike> CacheSelector<T, U19> for LTDCache<T> {
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
    pub coefficient_index_to_powers: Vec<Vec<u8>>,
    pub reduced_coefficient_index_to_powers: Vec<Vec<u8>>,
    pub reduced_blocks: Vec<usize>,
    pub coefficients_modified: bool,
    pub non_empty_coeff_map_to_reduced_numerator: Vec<Vec<(usize, usize, Complex<f64>)>>,
    pub coeff_map_to_reduced_numerator: Vec<Vec<usize>>,
}
#[derive(Debug, Clone, Default)]
pub struct ReducedLTDNumerator<T: FloatLike> {
    pub coefficients: Vec<Complex<T>>,
    pub n_loops: usize,
    pub max_rank: usize,
    pub size: usize,
    //pub coefficient_index_map: Vec<(usize, usize)>,
    //pub coefficient_index_to_powers: Vec<[u8; MAX_LOOP]>,
    pub reduced_coefficient_index_to_powers: Vec<Vec<u8>>,
    pub reduced_blocks: Vec<usize>,
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
    #[serde(skip_deserializing)]
    pub partial_fractioning_multiloops: PartialFractioningMultiLoops,
    #[serde(default)]
    pub numerator_tensor_coefficients_sparse: Vec<(Vec<usize>, (f64, f64))>,
    #[serde(default)]
    pub numerator_tensor_coefficients: Vec<(f64, f64)>,
    #[serde(default)]
    pub amplitude: Amplitude,
    #[serde(default)]
    pub loop_momentum_map: Vec<(Vec<i8>, Vec<i8>)>,
    #[serde(skip_deserializing)]
    pub socp_problem: SOCPProblem,
    #[serde(skip_deserializing)]
    pub subspaces: Vec<(usize, usize, Vec<(usize, usize)>)>,
}

impl Topology {
    pub fn from_file(
        filename: &str,
        settings: &Settings,
    ) -> Result<HashMap<String, Topology>, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open topology file {}", filename))
            .suggestion("Does the path exist?")?;

        let mut topologies: Vec<Topology> = serde_yaml::from_reader(f)
            .wrap_err("Could not parse topology file")
            .suggestion("Is it a correct yaml file")?;

        for t in &mut topologies {
            t.settings = settings.clone();
        }

        Ok(topologies
            .into_iter()
            .map(|t| (t.name.clone(), t))
            .collect())
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

    /// Create a rotated version of this topology. The axis needs to be normalized.
    fn rotate(&self, angle: float, axis: (float, float, float)) -> Topology {
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
        rotated_topology.name = rotated_topology.name + "_rot";
        rotated_topology.rotation_matrix = rot_matrix.clone();

        for e in &mut rotated_topology.external_kinematics {
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

        for ll in &mut rotated_topology.loop_lines {
            for p in &mut ll.propagators {
                let old_x = float::from_f64(p.q.x).unwrap();
                let old_y = float::from_f64(p.q.y).unwrap();
                let old_z = float::from_f64(p.q.z).unwrap();
                p.q.x = (rot_matrix[0][0] * old_x
                    + rot_matrix[0][1] * old_y
                    + rot_matrix[0][2] * old_z)
                    .to_f64()
                    .unwrap();
                p.q.y = (rot_matrix[1][0] * old_x
                    + rot_matrix[1][1] * old_y
                    + rot_matrix[1][2] * old_z)
                    .to_f64()
                    .unwrap();
                p.q.z = (rot_matrix[2][0] * old_x
                    + rot_matrix[2][1] * old_y
                    + rot_matrix[2][2] * old_z)
                    .to_f64()
                    .unwrap();
            }
        }

        for surf in &mut rotated_topology.surfaces {
            let old_x = surf.shift.x;
            let old_y = surf.shift.y;
            let old_z = surf.shift.z;
            surf.shift.x =
                rot_matrix[0][0] * old_x + rot_matrix[0][1] * old_y + rot_matrix[0][2] * old_z;
            surf.shift.y =
                rot_matrix[1][0] * old_x + rot_matrix[1][1] * old_y + rot_matrix[1][2] * old_z;
            surf.shift.z =
                rot_matrix[2][0] * old_x + rot_matrix[2][1] * old_y + rot_matrix[2][2] * old_z;
        }

        // now rotate the fixed deformation vectors
        for d_lim in &mut rotated_topology.fixed_deformation {
            for d in &mut d_lim.deformation_per_overlap {
                for source in &mut d.deformation_sources {
                    let old_x = source.x;
                    let old_y = source.y;
                    let old_z = source.z;
                    source.x = rot_matrix[0][0] * old_x
                        + rot_matrix[0][1] * old_y
                        + rot_matrix[0][2] * old_z;
                    source.y = rot_matrix[1][0] * old_x
                        + rot_matrix[1][1] * old_y
                        + rot_matrix[1][2] * old_z;
                    source.z = rot_matrix[2][0] * old_x
                        + rot_matrix[2][1] * old_y
                        + rot_matrix[2][2] * old_z;
                }
            }
        }
        // now rotate the numerators
        rotated_topology.numerator = rotated_topology.numerator.rotate(rot_matrix);
        if rotated_topology.settings.general.use_amplitude {
            for diag in rotated_topology.amplitude.diagrams.iter_mut() {
                diag.numerator = diag.numerator.rotate(rot_matrix);
            }
        }

        rotated_topology
    }
}

impl IntegrandImplementation for Topology {
    type Cache = LTDCacheAllPrecisions;

    fn create_stability_check(&self, num_checks: usize) -> Vec<Topology> {
        let mut rng = rand::thread_rng();
        let mut topologies = vec![];

        for _ in 0..num_checks {
            let angle =
                float::from_f64(rng.gen::<f64>() * 2.).unwrap() * <float as FloatConst>::PI();
            let mut rv = (
                float::from_f64(rng.gen()).unwrap(),
                float::from_f64(rng.gen()).unwrap(),
                float::from_f64(rng.gen()).unwrap(),
            ); // rotation axis
            let inv_norm = (rv.0 * rv.0 + rv.1 * rv.1 + rv.2 * rv.2).sqrt().inv();
            rv = (rv.0 * inv_norm, rv.1 * inv_norm, rv.2 * inv_norm);

            topologies.push(self.rotate(angle, rv));
        }

        topologies
    }

    fn get_target(&self) -> Option<Complex<f64>> {
        if self.analytical_result_real.is_some() || self.analytical_result_imag.is_some() {
            return Some(Complex::new(
                self.analytical_result_real.unwrap_or(0.),
                self.analytical_result_imag.unwrap_or(0.),
            ));
        }
        None
    }

    fn set_partial_fractioning(&mut self, enable: bool) {
        if enable {
            self.settings.general.partial_fractioning_threshold = 1e-99;
        } else {
            self.settings.general.partial_fractioning_threshold = -1.;
        }
    }

    fn create_grid(&self) -> Grid {
        Grid::ContinuousGrid(ContinuousGrid::new(
            self.n_loops,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
        ))
    }

    #[inline]
    fn evaluate_float<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut LTDCacheAllPrecisions,
        _events: Option<&mut EventManager>,
    ) -> Complex<float> {
        self.evaluate(x, cache.get())
    }

    #[inline]
    fn evaluate_f128<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut LTDCacheAllPrecisions,
        _events: Option<&mut EventManager>,
    ) -> Complex<f128::f128> {
        self.evaluate(x, cache.get())
    }

    #[inline]
    fn create_cache(&self) -> LTDCacheAllPrecisions {
        LTDCacheAllPrecisions::new(self)
    }

    fn set_precision(&mut self, _prec: usize) {
        // ignored
    }
}

impl Topology {
    pub fn determine_ellipsoid_overlap_structure(
        &mut self,
        update_excluded_surfaces: bool,
    ) -> Vec<FixedDeformationLimit> {
        if self.settings.general.debug > 1 {
            println!("Determining overlap structure for {}", self.name);
        }
        if self
            .surfaces
            .iter()
            .enumerate()
            .filter_map(|(i, s)| {
                if self.settings.general.debug > 1 {
                    self.print_surface(&s);
                }

                if s.exists && s.surface_type == SurfaceType::Ellipsoid && s.group == i {
                    Some(i)
                } else {
                    None
                }
            })
            .count()
            == 0
        {
            if self.settings.general.debug > 1 {
                println!("  | no ellipsoids");
            }
            return vec![];
        }

        // TODO: prevent allocating so many times
        let mut fixed_deformation = vec![];
        let mut subspaces = mem::replace(&mut self.subspaces, vec![]);
        let mut esurfaces = Vec::with_capacity(self.surfaces.len());
        let mut non_existing_e_surfaces = Vec::with_capacity(10);
        let mut constant_e_surfaces = Vec::with_capacity(10);

        for subspace in &subspaces {
            esurfaces.clear();
            non_existing_e_surfaces.clear();
            constant_e_surfaces.clear();

            if self.settings.general.debug > 0 {
                println!("  | doing subspace {:?}", subspace);
            }
            if self.settings.general.debug > 1 {
                println!(
                    "  | subspace basis: {:#?}",
                    self.ltd_cut_options[subspace.0][subspace.1]
                );
            }

            // rewrite all E-surfaces in the subspace basis, fill in the constraint,
            // and check if the E-surfaces still exist
            for i in 0..self.surfaces.len() {
                let s = &self.surfaces[i];
                if !s.exists || s.surface_type != SurfaceType::Ellipsoid || s.group != i {
                    continue;
                }

                let mut foci = Vec::with_capacity(s.id.len());
                let mut energy_shift = s.shift.t.multiply_sign(s.delta_sign);
                // construct a new matrix that transforms the momenta in the ellipsoid to the ellipsoid basis
                // where only one focus has a shift
                for ((foc_ll, foc_p), _, _) in &s.id {
                    let prop = &self.loop_lines[*foc_ll].propagators[*foc_p];
                    let mut shift = prop.q;

                    // transform the signature using the transposed cb_to_lmb matrix
                    let mut new_signature = vec![0; self.n_loops];
                    for (i, new_sig) in new_signature.iter_mut().enumerate() {
                        for (j, old_sig) in prop.signature.iter().enumerate() {
                            *new_sig +=
                                self.cb_to_lmb_mat[subspace.0][j * self.n_loops + i] * old_sig;
                        }
                    }

                    // incorporate the shift from the basis transformation
                    // all fixed components of the subspace basis should be removed from the signature
                    let mut subspace_basis_index = 0;
                    for (ll, c) in self.ltd_cut_options[subspace.0][subspace.1]
                        .iter()
                        .enumerate()
                    {
                        if let Cut::PositiveCut(prop_index) | Cut::NegativeCut(prop_index) = *c {
                            shift -= self.loop_lines[ll].propagators[prop_index]
                                .q
                                .multiply_sign(new_signature[subspace_basis_index]);

                            if subspace.2.contains(&(ll, prop_index)) {
                                new_signature[subspace_basis_index] = 0;
                            }

                            subspace_basis_index += 1;
                        }
                    }

                    if new_signature.iter().all(|i| *i == 0) {
                        energy_shift += (shift.spatial_squared() + prop.m_squared).sqrt();
                    } else {
                        foci.push(Focus {
                            id: prop.id,
                            signature: new_signature,
                            shift,
                            m_squared: self.loop_lines[*foc_ll].propagators[*foc_p].m_squared,
                        });
                    }
                }

                if foci.len() == 0 {
                    // ellipsoid is just a number!
                    if energy_shift < 0. {
                        // the constraint is on the inside of the ellipsoid, so we consider it existing
                        constant_e_surfaces.push(i);
                    } else {
                        non_existing_e_surfaces.push(i);
                    }
                    continue;
                }

                let subspace_signature_matrix: Vec<f64> = foci[1..]
                    .iter()
                    .map(|f| f.signature.iter().map(|c| *c as f64))
                    .flatten()
                    .collect();

                // take the right inverse of the reduced signature matrix
                let surf_basis_len = foci.len() - 1;
                let b = na::DMatrix::from_row_slice(
                    surf_basis_len,
                    self.n_loops,
                    &subspace_signature_matrix,
                );

                // take the right inverse
                let right_inverse = match (b.clone() * b.transpose()).try_inverse() {
                    Some(inv) => b.transpose() * inv,
                    None => {
                        if self.settings.general.debug > 1 {
                            println!("  | surface {} factorizes", i);
                        }
                        // if the right inverse cannot be computed, it means that the e-surface in the subspace
                        // factorizes into the sum of two or more e-surfaces:
                        // |k|+|l|+|m|+|k-l+m|+c under constraint k-l=p goes to
                        // |k|+|k-p|+|m|+|m+p+|+c
                        // in this case, we use an SOCP to test for its existence and get a point on the inside.
                        let mut surface = ESurface {
                            id: i,
                            foci,
                            shift: LorentzVector::from_args(energy_shift, 0., 0., 0.),
                            focus_basis_signature: vec![], // TODO: use Option
                            focus_basis_mat: vec![],
                            point_on_the_inside: vec![],
                        };

                        self.construct_socp_problem(&[0], std::slice::from_ref(&surface), true);
                        self.socp_problem.initialize_workspace_ecos();
                        let r = self.socp_problem.solve_ecos();
                        let exists = r == 0 || r == -2 || r == 10;
                        self.socp_problem.finish_ecos();

                        if exists {
                            let o = self.build_deformation_overlap_info(
                                &[0],
                                std::slice::from_ref(&surface),
                            );

                            surface.point_on_the_inside = o.deformation_sources;
                            esurfaces.push(surface);
                        } else {
                            non_existing_e_surfaces.push(i);
                        }

                        continue;
                    }
                };

                // note: this transpose is ONLY there because the iteration is column-wise instead of row-wise
                let mat: Vec<f64> = right_inverse.transpose().iter().cloned().collect();

                // now check if the surface still exists
                let mass_sum: f64 = foci.iter().map(|f| f.m_squared.sqrt()).sum();

                // determine the surface sign in the focus basis
                let mut signature = vec![0; surf_basis_len];
                for (i, s) in signature.iter_mut().enumerate() {
                    let mut sig_f = 0.;
                    for (j, depsig) in foci[0].signature.iter().enumerate() {
                        // use transpose
                        sig_f += *depsig as f64 * mat[j * surf_basis_len + i];
                    }

                    assert!(
                        (sig_f - sig_f.round()).abs() < 1e-8,
                        "Signature is not -1,0,1: {}",
                        sig_f
                    );
                    *s = sig_f as i8;
                }

                // determine the 4-vector surface shift using the matrix
                let mut surface_shift = foci[0].shift;
                for (&sign, f) in signature.iter().zip_eq(&foci[1..]) {
                    surface_shift -= f.shift.multiply_sign(sign);
                }

                surface_shift.t = energy_shift;

                if surface_shift.square() - mass_sum.powi(2) >= 1e-8 * self.e_cm_squared
                    && surface_shift.t < 0.
                {
                    if self.settings.general.debug > 0 {
                        let r = surface_shift.square() - mass_sum.powi(2);
                        if r < 1e-6 * self.e_cm_squared {
                            println!("Small ellipsoid {} detected: {}", i, r);
                        }
                    }

                    let mut surface = ESurface {
                        id: i,
                        foci,
                        shift: surface_shift,
                        focus_basis_signature: signature,
                        focus_basis_mat: mat,
                        point_on_the_inside: vec![],
                    };

                    let p = self.get_point_inside_surface::<f64>(&surface);
                    surface.point_on_the_inside.extend(&p[..self.n_loops]);
                    esurfaces.push(surface);
                } else {
                    non_existing_e_surfaces.push(i);
                }
            }

            let mut overlap = self.find_overlap_structure(&esurfaces);

            // map the source back to the loop momentum basis
            for o in &mut overlap {
                let mut sources_lmb = vec![LorentzVector::default(); self.n_loops];

                for (source_lmb, r) in sources_lmb
                    .iter_mut()
                    .zip_eq(self.cb_to_lmb_mat[subspace.0].chunks_exact(self.n_loops))
                {
                    for (i, (source, c)) in o.deformation_sources.iter().zip_eq(r).enumerate() {
                        // find the shift
                        let llprop = self.ltd_cut_options[subspace.0][subspace.1]
                            .iter()
                            .enumerate()
                            .filter_map(|(i, v)| match *v {
                                Cut::PositiveCut(j) => Some((i, j)),
                                Cut::NegativeCut(j) => Some((i, j)),
                                Cut::NoCut => None,
                            })
                            .nth(i)
                            .unwrap();
                        let shift = &self.loop_lines[llprop.0].propagators[llprop.1].q;

                        *source_lmb += (source - shift).multiply_sign(*c);
                        source_lmb.t = 0.;
                    }
                }

                o.deformation_sources.clear();
                o.deformation_sources.extend(&sources_lmb[..self.n_loops]);

                // add all non-existing E-surfaces to the exclusion
                o.excluded_surface_indices.extend(&non_existing_e_surfaces);
                o.overlap.as_mut().map(|ov| ov.extend(&constant_e_surfaces));
            }

            fixed_deformation.push(FixedDeformationLimit {
                deformation_per_overlap: overlap,
                excluded_propagators: subspace.2.clone(),
            });
        }

        mem::swap(&mut subspaces, &mut self.subspaces);

        if self.settings.general.debug > 1 {
            println!("  | deformation={:#?}", fixed_deformation);
        }

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

        if self.settings.general.debug > 0 && self.fixed_deformation.len() > 0 {
            self.compare_fixed_deformation(&mut fixed_deformation);
        }

        fixed_deformation
    }

    /// Compare a fixed deformation with the one from yaml.
    fn compare_fixed_deformation(&mut self, fixed_deformation: &mut [FixedDeformationLimit]) {
        // check the deformation with the one from the yaml file
        assert_eq!(fixed_deformation.len(), self.fixed_deformation.len());

        for yf in &mut self.fixed_deformation {
            for yf_ov in &mut yf.deformation_per_overlap {
                yf_ov.excluded_surface_indices.sort();
            }
        }

        for f in fixed_deformation.iter_mut() {
            for f_ov in &mut f.deformation_per_overlap {
                f_ov.excluded_surface_indices.sort();
            }
        }

        for f in fixed_deformation.iter() {
            let yf = self
                .fixed_deformation
                .iter_mut()
                .filter(|yf| f.excluded_propagators == yf.excluded_propagators)
                .next()
                .unwrap_or_else(|| {
                    panic!(
                        "Could not find overlap with propagator exclusion {:?}",
                        f.excluded_propagators
                    )
                });

            assert_eq!(
                f.deformation_per_overlap.len(),
                yf.deformation_per_overlap.len()
            );
            for f_ov in &f.deformation_per_overlap {
                let yf_ov = yf
                    .deformation_per_overlap
                    .iter_mut()
                    .filter(|yf_ov| f_ov.excluded_surface_indices == yf_ov.excluded_surface_indices)
                    .next()
                    .unwrap_or_else(|| {
                        panic!(
                            "Could not find source with exclusion {:?} for subspace {:?}",
                            f_ov.excluded_surface_indices, f.excluded_propagators,
                        )
                    });

                // check if the sources are approximately the same
                for (f_s, yf_s) in f_ov
                    .deformation_sources
                    .iter()
                    .zip_eq(yf_ov.deformation_sources.iter())
                {
                    if (f_s - yf_s).spatial_squared() > 1e-5 {
                        println!(
                    "Warning: sources do not match precisely: rust={} vs yaml={} for subspace {:?}",
                    f_s, yf_s, f.excluded_propagators,
                );
                    }
                }
            }
        }
    }

    /// Construct the SOCP problem.
    fn construct_socp_problem(
        &mut self,
        ellipsoid_ids: &[usize],
        ellipsoids: &[ESurface],
        minimize: bool,
    ) {
        let mut p = &mut self.socp_problem;

        // TODO: change the n_loops to the number of specific variables
        let multiplier = if minimize { 6 * self.n_loops } else { 1 };
        p.radius_computation = minimize;

        p.cone.l = (ellipsoid_ids.len() * multiplier) as scs::scs_int;

        let mut b_index = 0;

        for _ in 0..multiplier {
            for &e_id in ellipsoid_ids {
                p.b[b_index] = -ellipsoids[e_id].shift.t;
                b_index += 1;
            }
        }

        let mut lm_tag = [false; MAX_LOOP];

        // construct the focus part of the shift
        // TODO: do not generate the focus more than once if k[m/3] does not appear in the focus!!!
        let mut focus_count = 0;
        for m in 0..multiplier {
            for ll in &self.loop_lines {
                for prop in &ll.propagators {
                    // check if the focus occurs in an ellipsoid
                    'el_loop: for &e_id in ellipsoid_ids {
                        for (foc_index, f) in ellipsoids[e_id].foci.iter().enumerate() {
                            if f.id == prop.id {
                                for (tag, &sig) in lm_tag.iter_mut().zip(&f.signature) {
                                    *tag |= sig != 0;
                                }

                                p.b[b_index] = 0.;
                                b_index += 1;
                                if f.m_squared != 0. {
                                    p.b[b_index] = f.m_squared.sqrt();
                                    b_index += 1;
                                    p.q[focus_count] = 5;
                                } else {
                                    p.q[focus_count] = 4;
                                };
                                p.b[b_index] = f.shift.x;
                                b_index += 1;
                                p.b[b_index] = f.shift.y;
                                b_index += 1;
                                p.b[b_index] = f.shift.z;
                                b_index += 1;

                                p.focus_list[focus_count] = (prop.id, e_id, foc_index, m);
                                focus_count += 1;
                                break 'el_loop;
                            }
                        }
                    }
                }
            }
        }

        let mut var_count = 0;
        for (vm, &x) in p.var_map.iter_mut().zip(lm_tag.iter()) {
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
            // minimize the radius. ECOS very rarely gives a proof of infeasibility if we maximize instead
            p.c[3 * var_count] = 1.;
            3 * var_count + 1
        } else {
            3 * var_count
        };

        // write the ellipsoid constraints: f1 + f2 < c, etc
        let mut row_counter = 0; // the start of the focus information
        for m in 0..multiplier {
            for &e_id in ellipsoid_ids {
                for f in ellipsoids[e_id].foci.iter() {
                    let foc_index = p
                        .focus_list
                        .iter()
                        .position(|x| x.0 == f.id && x.3 == m)
                        .unwrap();
                    p.a_dense[row_counter * width + focus_start + foc_index] = 1.;
                }
                row_counter += 1;
            }
        }

        // now for every focus f, add the norm constraints |k_x+p_x,m| < f
        for (foc_index, (_, foc_e_surf, foc_e_surf_foc_index, radius_index)) in
            p.focus_list[..focus_count].iter().enumerate()
        {
            let focus = &ellipsoids[*foc_e_surf].foci[*foc_e_surf_foc_index];

            // set the linear equation
            p.a_dense[row_counter * width + focus_start + foc_index] = -1.0;
            row_counter += 1;
            if focus.m_squared != 0. {
                // the mass is an empty line
                row_counter += 1;
            }

            let (rad_lm, rad_dir, rad_sign) =
                (radius_index / 6, radius_index % 3, radius_index % 2);

            for dir_index in 0..3 {
                for (var_index, &v) in focus.signature.iter().enumerate() {
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
    fn evaluate_surface<T: FloatLike>(&self, loop_momenta: &[LorentzVector<T>], s: &ESurface) -> T {
        let mut eval = Into::<T>::into(s.shift.t);
        for f in &s.foci {
            // get the loop momentum
            let mut mom = f.shift.cast();
            for (m, sig) in loop_momenta[..self.n_loops].iter().zip(&f.signature) {
                mom += m.multiply_sign(*sig);
            }
            eval += (mom.spatial_squared() + Into::<T>::into(f.m_squared)).sqrt();
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
        all_ellipsoids: &[ESurface],
    ) -> FixedDeformationOverlap {
        // TODO: recycle from pre-allocated overlap structures
        let mut var_count = 0;
        let mut sources = vec![LorentzVector::default(); self.n_loops];

        // map the source back from the loop momenta present in the ellipsoids to the full space
        for (source, &i) in sources
            .iter_mut()
            .zip_eq(&self.socp_problem.var_map[..self.n_loops])
        {
            if i < self.n_loops {
                *source = LorentzVector::from_args(
                    0.,
                    self.socp_problem.sol_x[i * 3],
                    self.socp_problem.sol_x[i * 3 + 1],
                    self.socp_problem.sol_x[i * 3 + 2],
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

        if self.socp_problem.radius_computation {
            radius = -self.socp_problem.sol_x[var_count * 3];
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

        let overlap = ellipsoid_list
            .iter()
            .map(|e_index| all_ellipsoids[*e_index].id)
            .collect();

        FixedDeformationOverlap {
            deformation_sources: sources,
            excluded_surface_ids: vec![],
            excluded_surface_indices,
            overlap: Some(overlap),
            radius,
        }
    }

    /// Find a point on a surface along a line and return the scale:
    /// `E(dir * t + k) = 0`, returning `t` closest to `t_start`.
    fn get_scale_to_reach_surface<T: FloatLike>(
        &self,
        s: &ESurface,
        k: &[LorentzVector<T>],
        dir: &[LorentzVector<T>],
        t_start: T,
    ) -> Option<T> {
        let mut t = t_start;
        let tolerance = Into::<T>::into(1e-15 * self.e_cm_squared.sqrt());
        for it in 0..30 {
            let mut f = Into::<T>::into(s.shift.t);
            let mut df = T::zero();

            for focus in &s.foci {
                let shift = focus.shift.cast();
                let k = utils::evaluate_signature(&focus.signature, k);
                let dir = utils::evaluate_signature(&focus.signature, dir);
                let energy = ((k + dir * t + shift).spatial_squared()
                    + Into::<T>::into(focus.m_squared))
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
    fn get_point_inside_surface<T: FloatLike>(&self, s: &ESurface) -> [LorentzVector<T>; MAX_LOOP] {
        // in the focus basis, the first focus has the sum of basis momenta
        let mass_sum: T = s
            .foci
            .iter()
            .map(|f| Into::<T>::into(f.m_squared).sqrt())
            .sum();

        let mut cb_on_foci = [LorentzVector::default(); MAX_LOOP];

        for ((m, f), sig) in cb_on_foci
            .iter_mut()
            .zip(&s.foci[1..])
            .zip(&s.focus_basis_signature)
        {
            if !mass_sum.is_zero() && !f.m_squared.is_zero() {
                *m = -s.shift.cast().multiply_sign(*sig) * Into::<T>::into(f.m_squared).sqrt()
                    / mass_sum;
            }
        }

        // do the basis transformation from the cut basis
        let mut ll_on_foci = [LorentzVector::default(); MAX_LOOP];
        for (lm, r) in ll_on_foci[..self.n_loops]
            .iter_mut()
            .zip_eq(s.focus_basis_mat.chunks_exact(s.foci.len() - 1))
        {
            for ((sign, cm), f) in r
                .iter()
                .zip_eq(&cb_on_foci[..s.foci.len() - 1])
                .zip_eq(&s.foci[1..])
            {
                *lm += (cm - f.shift.cast()) * Into::<T>::into(*sign);
            }
            lm.t = T::zero();
        }

        debug_assert!(self.evaluate_surface(&ll_on_foci, s) < T::zero());
        ll_on_foci
    }

    /// Get the normal at the point on an E-surface. We always assume that the surface is positive.
    fn get_normal<T: FloatLike>(
        &self,
        surface: &ESurface,
        k: &[LorentzVector<T>],
        normalize: bool,
    ) -> [LorentzVector<T>; MAX_LOOP] {
        let mut df = [LorentzVector::default(); MAX_LOOP];

        for f in &surface.foci {
            let shift = f.shift.cast();
            let mut mom_part = utils::evaluate_signature(&f.signature, k) + shift;
            mom_part.t = T::zero();
            let energy = (mom_part.spatial_squared() + Into::<T>::into(f.m_squared)).sqrt();

            if energy == T::zero() || Float::is_nan(energy) {
                continue;
            }

            for (d, s) in df.iter_mut().zip(&f.signature) {
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
    fn overlap_heuristic(&self, surf1: &ESurface, surf2: &ESurface) -> bool {
        let mut x1 = [LorentzVector::default(); MAX_LOOP];
        let mut x2 = [LorentzVector::default(); MAX_LOOP];
        x1[..self.n_loops].copy_from_slice(&surf1.point_on_the_inside);
        x2[..self.n_loops].copy_from_slice(&surf2.point_on_the_inside);

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
                        it, surf1.id, surf2.id
                    );
                }
            }
        }
        false
    }

    fn find_overlap_structure(
        &mut self,
        ellipsoid_list: &[ESurface],
    ) -> Vec<FixedDeformationOverlap> {
        let mut overlap_structure = vec![];

        let t = Instant::now();

        let mut problem_count = 0;
        let mut pair_overlap_count = 0;

        if self.settings.deformation.fixed.use_heuristic_centers {
            if ellipsoid_list.len() == 1 {
                // put ourselves on a focus
                let r: f64 = -self
                    .evaluate_surface(&ellipsoid_list[0].point_on_the_inside, &ellipsoid_list[0]);
                assert!(r > 0.);

                return vec![FixedDeformationOverlap {
                    deformation_sources: ellipsoid_list[0].point_on_the_inside.clone(),
                    excluded_surface_ids: vec![],
                    excluded_surface_indices: vec![],
                    overlap: Some(vec![ellipsoid_list[0].id]),
                    radius: r,
                }];
            }

            // the origin is often inside all surfaces, so do a quick test for that
            let mut origin_inside_radius = self.e_cm_squared;
            let origin = [LorentzVector::default(); MAX_LOOP];
            for e in ellipsoid_list {
                let r: f64 = -self.evaluate_surface(&origin[..self.n_loops], e);
                if r <= 0. || r * r < 1e-10 * self.e_cm_squared {
                    origin_inside_radius = -1.;
                    break;
                }
                if r < origin_inside_radius {
                    origin_inside_radius = r;
                }
            }

            if origin_inside_radius > 0. {
                if self.settings.general.debug > 0 {
                    println!("Origin inside for all E-surface");
                }

                return vec![FixedDeformationOverlap {
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
                for foc1 in &ellipsoid_list[i].foci {
                    for foc2 in &ellipsoid_list[j].foci {
                        if foc1
                            .signature
                            .iter()
                            .zip(foc2.signature.iter())
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

                self.construct_socp_problem(&[i, j], ellipsoid_list, false);

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
            while self.next_set(
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

                        let r = -self.evaluate_surface(
                            &ellipsoid_list[o1].point_on_the_inside,
                            &ellipsoid_list[o2],
                        );

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
                    self.construct_socp_problem(&option_translated[..n], ellipsoid_list, false);

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

                        FixedDeformationOverlap {
                            deformation_sources: ellipsoid_list[point_e_id].point_on_the_inside
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

                            self.construct_socp_problem(
                                &option_translated[..n],
                                ellipsoid_list,
                                true,
                            );

                            self.socp_problem.initialize_workspace_ecos();
                            let r = self.socp_problem.solve_ecos();
                            if !(r == 0 || r == -2 || r == 10) {
                                println!(
                                    "ECOS cannot opimize radius for {}: {}. Falling back to default overlap finding.",
                                    self.name, r
                                );

                                self.socp_problem.finish_ecos();
                                self.construct_socp_problem(
                                    &option_translated[..n],
                                    ellipsoid_list,
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
                            self.construct_socp_problem(
                                &option_translated[..n],
                                ellipsoid_list,
                                false,
                            );
                            self.socp_problem.initialize_workspace_ecos();
                            let r = self.socp_problem.solve_ecos();
                            self.socp_problem.finish_ecos();
                            ecos_time += Instant::now().duration_since(ti).as_nanos();
                            problem_count += 1;
                            assert!(r == 0 || r == -2 || r == 10, "ECOS cannot find overlap");
                        }

                        self.build_deformation_overlap_info(&option_translated[..n], ellipsoid_list)
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
        *p.settings = *self.settings;
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
