use amplitude::Amplitude;
use arrayvec::ArrayVec;
use dual_num::{DualN, Scalar, U10, U13, U16, U19, U4, U7};
use float;
use itertools::Itertools;
use num::Complex;
use num_traits::{Float, FloatConst};
use num_traits::{One, Signed, Zero};
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use utils;
use utils::Signum;
use vector::{LorentzVector, RealNumberLike};
use {FloatLike, PythonNumerator, Settings, MAX_LOOP};

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
    pub propagators: HashMap<(usize, usize), Complex<T>>,
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
            propagators: HashMap::new(),
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
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCuts {
    pub cut_names: Vec<String>,
    pub cut_signs: Vec<i8>,
    pub cut_signature: Vec<(Vec<i8>, Vec<i8>)>,
    pub subgraph_left: Topology,
    pub subgraph_right: Topology,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopology {
    pub n_loops: usize,
    pub cutkosky_cuts: Vec<CutkoskyCuts>,
    #[serde(skip_deserializing)]
    pub settings: Settings,
}

impl SquaredTopology {
    pub fn from_file(filename: &str, settings: &Settings) -> SquaredTopology {
        let f = File::open(filename).expect("Could not open squared topology file");

        let mut squared_topo: SquaredTopology = serde_yaml::from_reader(f).unwrap();
        squared_topo.settings = settings.clone();
        for cutkosky_cuts in &mut squared_topo.cutkosky_cuts {
            cutkosky_cuts.subgraph_left.settings = settings.clone();
            cutkosky_cuts.subgraph_right.settings = settings.clone();
            cutkosky_cuts.subgraph_left.process();
            cutkosky_cuts.subgraph_right.process();
        }

        squared_topo
    }

    fn evaluate_signature<T: RealNumberLike>(
        signature: &(Vec<i8>, Vec<i8>),
        external_momenta: &[LorentzVector<T>],
        loop_momenta: &[LorentzVector<T>],
    ) -> LorentzVector<T> {
        let mut cut_momentum = LorentzVector::default();
        for (&sign, mom) in signature.0.iter().zip_eq(loop_momenta) {
            if sign != 0 {
                cut_momentum += mom.multiply_sign(sign);
            }
        }

        for (&sign, mom) in signature.1.iter().zip_eq(external_momenta) {
            if sign != 0 {
                cut_momentum += mom.multiply_sign(sign);
            }
        }
        cut_momentum
    }

    pub fn create_caches<T: FloatLike + Into<f64>>(&self) -> Vec<Vec<LTDCache<T>>> {
        let mut caches = vec![];
        for cutkosky_cuts in &self.cutkosky_cuts {
            caches.push(vec![
                LTDCache::<T>::new(&cutkosky_cuts.subgraph_left),
                LTDCache::<T>::new(&cutkosky_cuts.subgraph_right),
            ]);
        }
        caches
    }

    pub fn evaluate<T: FloatLike + Into<f64>>(
        &mut self,
        external_momenta: &mut [LorentzVector<T>],
        loop_momenta: &[LorentzVector<T>],
        caches: &mut [Vec<LTDCache<T>>],
    ) -> Complex<T> {
        let mut cut_momenta = [LorentzVector::default(); MAX_LOOP + 4]; // FIXME: bound may be too small
        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_LOOP];
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]> =
            (0..MAX_LOOP).map(|_| LorentzVector::default()).collect();
        let mut result = Complex::zero();
        for (cutkosky_cuts, cache) in self.cutkosky_cuts.iter_mut().zip_eq(caches) {
            let mut cut_result = Complex::one();
            let mut q0 = T::zero();

            // evaluate the cuts
            for ((cut, cut_sig), cut_sign) in cut_momenta[..cutkosky_cuts.cut_signature.len()]
                .iter_mut()
                .zip(cutkosky_cuts.cut_signature.iter())
                .zip(cutkosky_cuts.cut_signs.iter())
            {
                *cut = SquaredTopology::evaluate_signature(cut_sig, external_momenta, loop_momenta);
                let energy = cut.spatial_distance();
                cut.t = energy.multiply_sign(*cut_sign);
                cut_result *= num::Complex::new(T::zero(), -<T as FloatConst>::PI() / energy); // add (2 pi i)/(2E) for every cut
                q0 += energy;
            }

            if self.settings.general.debug >= 1 {
                println!("Cut {:?}:", cutkosky_cuts.cut_names);
                println!("  | 1/Es = {}", cut_result);
                println!("  | q0 = {}", q0);
            }

            // set 0th component of external momenta to the sum of the cuts
            external_momenta[0].t = q0;
            external_momenta[1].t = q0; // FIXME: we assume 1->1 here and that it is outgoing

            // set the shifts, which are expressed in the cut basis
            let mut subgraphs = [
                &mut cutkosky_cuts.subgraph_left,
                &mut cutkosky_cuts.subgraph_right,
            ];
            for subgraph in &mut subgraphs {
                for ll in &mut subgraph.loop_lines {
                    for p in &mut ll.propagators {
                        p.q = SquaredTopology::evaluate_signature(
                            &p.parametric_shift,
                            external_momenta,
                            &cut_momenta[..cutkosky_cuts.cut_signature.len()],
                        )
                        .map(|x| x.into());
                    }
                }
            }

            //subgraphs[0].process();
            //subgraphs[1].process();
            subgraphs[0].e_cm_squared = Into::<f64>::into(q0 * q0);
            subgraphs[1].e_cm_squared = Into::<f64>::into(q0 * q0);

            // evaluate
            for (subgraph, subgraph_cache) in subgraphs.iter_mut().zip_eq(cache.iter_mut()) {
                // do the loop momentum map, which is expressed in the loop momentum basis
                // the time component should not matter here
                for (slm, lmm) in subgraph_loop_momenta[..subgraph.n_loops]
                    .iter_mut()
                    .zip_eq(&subgraph.loop_momentum_map)
                {
                    *slm =
                        SquaredTopology::evaluate_signature(lmm, external_momenta, &loop_momenta);
                }

                // for now, we do not deform
                for (kd, k) in k_def[..subgraph.n_loops]
                    .iter_mut()
                    .zip_eq(&subgraph_loop_momenta[..subgraph.n_loops])
                {
                    *kd = k.map(|x| Complex::new(x, T::zero()));
                }

                if subgraph
                    .compute_complex_cut_energies(&k_def[..subgraph.n_loops], subgraph_cache)
                    .is_err()
                {
                    panic!("NaN on cut energy");
                }

                let (mut res, kd) = if subgraph.loop_lines.len() > 0 {
                    subgraph
                        .evaluate_all_dual_integrands::<T, PythonNumerator>(
                            &subgraph_loop_momenta[..subgraph.n_loops],
                            k_def,
                            subgraph_cache,
                            &None,
                        )
                        .unwrap()
                } else {
                    // if the graph has no propagators, it is one and not zero
                    (Complex::one(), k_def)
                };

                res *= utils::powi(
                    num::Complex::new(T::zero(), Into::<T>::into(-2.) * <T as FloatConst>::PI()),
                    subgraph.n_loops,
                ); // factor of delta cut

                k_def = kd;
                cut_result *= res;

                if self.settings.general.debug >= 1 {
                    println!(
                        "  | res {} ({}l) = {:e}",
                        subgraph.name, subgraph.n_loops, res
                    );
                }
            }

            cut_result *= utils::powi(
                num::Complex::new(
                    Into::<T>::into(1.)
                        / <T as Float>::powi(Into::<T>::into(2.) * <T as FloatConst>::PI(), 4),
                    T::zero(),
                ),
                self.n_loops,
            );

            if self.settings.general.debug >= 1 {
                println!("  | res = {:e}", cut_result);
            }

            result += cut_result;
        }

        if self.settings.general.debug >= 1 {
            println!("Final result = {:e}", result);
        }

        result
    }
}
