use dual_num::{DualN, Scalar, U10, U13, U4, U7};
use float;
use num::Complex;
use num_traits::Signed;
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use vector::{LorentzVector, RealNumberLike};
use Settings;

/// Ellipsoid and hyperboloid surfaces
#[derive(Debug, Clone)]
pub struct Surface {
    pub group: usize,
    pub ellipsoid: bool,
    pub cut_structure_index: usize,
    pub cut_option_index: usize,
    pub cut: Vec<Cut>,
    pub onshell_ll_index: usize,
    pub onshell_prop_index: usize,
    pub delta_sign: i8,
    pub sig_ll_in_cb: Vec<i8>,
    pub signs: Vec<i8>,
    pub shift: LorentzVector<float>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Propagators {
    #[serde(skip_deserializing)]
    pub id: usize, // global id
    pub m_squared: f64,
    pub q: LorentzVector<f64>,
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
    pub ellipsoid_eval: Vec<DualN<T, U>>,
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
            ellipsoid_eval: vec![DualN::default(); num_surfaces],
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
    pub complex_cut_energies: Vec<Complex<T>>,
    pub complex_prop_spatial: Vec<Complex<T>>,
    pub complex_loop_line_eval: Vec<Vec<[Complex<T>; 2]>>,
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
            complex_cut_energies: vec![Complex::default(); num_propagators],
            complex_prop_spatial: vec![Complex::default(); num_propagators],
            complex_loop_line_eval: topo
                .loop_lines
                .iter()
                .map(|ll| vec![[Complex::default(); 2]; ll.propagators.len()])
                .collect(),
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

#[derive(Debug, Clone, Deserialize)]
pub struct Topology {
    pub name: String,
    pub n_loops: usize,
    pub analytical_result_real: Option<f64>,
    pub analytical_result_imag: Option<f64>,
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
    #[serde(skip_deserializing)]
    pub rotation_matrix: [[float; 3]; 3],
    #[serde(skip_deserializing)]
    pub ellipsoids_not_in_cuts: Vec<Vec<Vec<usize>>>,
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

    /*pub fn build_loop_line(
        loop_momenta: &[(usize, bool)],
        ext: Vec<LorentzVector<f64>>,
        mass: Vec<f64>,
        compute_qs: bool,
    ) -> LoopLine {
        let qs = if compute_qs {
            // convention of Dario
            let mut qs = vec![ext[0]];
            if ext.len() > 1 {
                for e in &ext[1..ext.len() - 1] {
                    let n = qs.last().unwrap() + e;
                    qs.push(n);
                }
            }
            qs.push(LorentzVector::default());
            qs
        } else {
            // externals are already summed
            ext
        };

        let q_and_mass: Vec<_> = qs.into_iter().zip(mass.into_iter()).collect();

        LoopLine::new(loop_momenta, q_and_mass)
    }

    pub fn get_hardcoded_topology(topology: &str) -> Topology {
        match topology {
            "P1" => {
                // does not need deformation
                let p1 = LorentzVector::from_args(5.23923, -4.18858, 0.74966, -3.05669);
                let p2 = LorentzVector::from_args(6.99881, -2.93659, 5.03338, 3.87619);
                let m = 7.73358;
                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll = build_loop_line(&[(0, true)], vec![p1, p2, -p1 - p2], vec![m, m, m], true);
                (1, e_cm, vec![ll])
            }
            "P5" => {
                // does not need deformation
                let p1 = LorentzVector::from_args(31.54872, -322.40325, 300.53015, -385.58013);
                let p2 = LorentzVector::from_args(103.90430, 202.00974, -451.27794, -435.12848);
                let p3 = LorentzVector::from_args(294.76653, 252.88958, 447.09194, 311.71630);
                let m = 4.68481;
                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll = build_loop_line(
                    &[(0, true)],
                    vec![p1, p2, p3, -p1 - p2 - p3],
                    vec![m; 4],
                    true,
                );
                (1, e_cm, vec![ll])
            }
            "P3" => {
                //P3 in https://arxiv.org/pdf/1510.00187.pdf
                let p1 = LorentzVector::from_args(10.51284, 6.89159, -7.40660, -2.85795);
                let p2 = LorentzVector::from_args(6.45709, 2.46635, 5.84093, 1.22257);
                let m = 0.52559;
                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll = build_loop_line(&[(0, true)], vec![p1, p2, -p1 - p2], vec![m; 3], true);
                (1, e_cm, vec![ll])
            }
            "P4" => {
                //P4 in https://arxiv.org/pdf/1510.00187.pdf
                let p1 = LorentzVector::from_args(95.77004, 31.32025, -34.08106, -9.38565);
                let p2 = LorentzVector::from_args(94.54738, -53.84229, 67.11107, 45.56763);
                let m1 = 83.02643;
                let m2 = 76.12873;
                let m3 = 55.00359;
                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll =
                    build_loop_line(&[(0, true)], vec![p1, p2, -p1 - p2], vec![m1, m2, m3], true);
                (1, e_cm, vec![ll])
            }
            "P7" => {
                //Use Massive Momenta
                //Analytic Result:  -2.38766 e-10 - i 3.03080 e-10
                let p1 = LorentzVector::from_args(62.80274, -49.71968, -5.53340, -79.44048);
                let p2 = LorentzVector::from_args(48.59375, -1.65847, 34.91140, 71.89564);
                let p3 = LorentzVector::from_args(76.75934, -19.14334, -17.10279, 30.22959);
                let m = 9.82998;

                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll = build_loop_line(
                    &[(0, true)],
                    vec![p1, p2, p3, -p1 - p2 - p3],
                    vec![m, m, m, m],
                    true,
                );
                (1, e_cm, vec![ll])
            }
            "onshell-box" => {
                //Use Massless Momenta and Massless Propagators
                //Analytic Result (real): 6.08501e-3
                let p1 = LorentzVector::from_args(0.5, 0.5, 0.0, 0.0);
                let p2 = LorentzVector::from_args(0.5, -0.5, 0.0, 0.0);
                let p3 = LorentzVector::from_args(-0.5, 0.0, 0.5, 0.0);
                let m = 0.;
                println!(
                    "p1^2: {}, p2^2: {}, p3^3: {}",
                    p1.square(),
                    p2.square(),
                    p3.square()
                );
                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll = build_loop_line(
                    &[(0, true)],
                    vec![p1, p2, p3, -p1 - p2 - p3],
                    vec![m, m, m, m],
                    true,
                );
                (1, e_cm, vec![ll])
            }
            "massless-box" => {
                //All but p1 are massless and Massless Propagators
                //Analytic Result (real): ?
                let on_shell_flag = 16;
                println!("select PS: {:0.4b}", on_shell_flag);
                let (p1, p2, p3, _) = match on_shell_flag {
                    0b0000 => (
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            0.2516167274762227,
                            0.05186897514731689,
                            -0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            -0.2516167274762227,
                            -0.05186897514731689,
                            0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            -0.2422804773385922,
                            0.1018206163990297,
                            -0.1317288330743322,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            0.2422804773385922,
                            -0.1018206163990297,
                            0.1317288330743322,
                        ),
                    ),
                    0b1000 => (
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            0.2516167274762227,
                            0.05186897514731689,
                            -0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            -0.2516167274762227,
                            -0.05186897514731689,
                            0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            0.6111111111111111,
                            -0.3205069452819676,
                            0.1346960146655688,
                            -0.1742608664057115,
                        ),
                        LorentzVector::from_args(
                            0.3888888888888889,
                            0.3205069452819676,
                            -0.1346960146655688,
                            0.1742608664057115,
                        ),
                    ),
                    0b0100 => (
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            0.2516167274762227,
                            0.05186897514731689,
                            -0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            -0.2516167274762227,
                            -0.05186897514731689,
                            0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            -0.3662936517508201,
                            0.1539383024749358,
                            -0.1991552758922418,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            0.3662936517508201,
                            -0.1539383024749358,
                            0.1991552758922418,
                        ),
                    ),
                    0b1100 => (
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            0.2516167274762227,
                            0.05186897514731689,
                            -0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            -0.2516167274762227,
                            -0.05186897514731689,
                            0.02012167665288664,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            -0.4120803582196727,
                            0.1731805902843028,
                            -0.2240496853787720,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            0.4120803582196727,
                            -0.1731805902843028,
                            0.2240496853787720,
                        ),
                    ),
                    0b0010 => (
                        LorentzVector::from_args(
                            -0.6250000000000000,
                            0.3661561216082539,
                            0.07548044584409114,
                            -0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            -0.3750000000000000,
                            -0.3661561216082539,
                            -0.07548044584409114,
                            0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            -0.2422804773385922,
                            0.1018206163990297,
                            -0.1317288330743322,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            0.2422804773385922,
                            -0.1018206163990297,
                            0.1317288330743322,
                        ),
                    ),
                    0b1010 => (
                        LorentzVector::from_args(
                            -0.6250000000000000,
                            0.3661561216082539,
                            0.07548044584409114,
                            -0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            -0.3750000000000000,
                            -0.3661561216082539,
                            -0.07548044584409114,
                            0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            0.6111111111111111,
                            -0.3205069452819676,
                            0.1346960146655688,
                            -0.1742608664057115,
                        ),
                        LorentzVector::from_args(
                            0.3888888888888889,
                            0.3205069452819676,
                            -0.1346960146655688,
                            0.1742608664057115,
                        ),
                    ),
                    0b0110 => (
                        LorentzVector::from_args(
                            -0.6250000000000000,
                            0.3661561216082539,
                            0.07548044584409114,
                            -0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            -0.3750000000000000,
                            -0.3661561216082539,
                            -0.07548044584409114,
                            0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            -0.3662936517508201,
                            0.1539383024749358,
                            -0.1991552758922418,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            0.3662936517508201,
                            -0.1539383024749358,
                            0.1991552758922418,
                        ),
                    ),
                    0b1110 => (
                        LorentzVector::from_args(
                            -0.6250000000000000,
                            0.3661561216082539,
                            0.07548044584409114,
                            -0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            -0.3750000000000000,
                            -0.3661561216082539,
                            -0.07548044584409114,
                            0.02928134054272110,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            -0.4120803582196727,
                            0.1731805902843028,
                            -0.2240496853787720,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            0.4120803582196727,
                            -0.1731805902843028,
                            0.2240496853787720,
                        ),
                    ),
                    0b0001 => (
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            0.4271821418762962,
                            0.08806052015143966,
                            -0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            -0.4271821418762962,
                            -0.08806052015143966,
                            0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            -0.2422804773385922,
                            0.1018206163990297,
                            -0.1317288330743322,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            0.2422804773385922,
                            -0.1018206163990297,
                            0.1317288330743322,
                        ),
                    ),
                    0b1001 => (
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            0.4271821418762962,
                            0.08806052015143966,
                            -0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            -0.4271821418762962,
                            -0.08806052015143966,
                            0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            0.6111111111111111,
                            -0.3205069452819676,
                            0.1346960146655688,
                            -0.1742608664057115,
                        ),
                        LorentzVector::from_args(
                            0.3888888888888889,
                            0.3205069452819676,
                            -0.1346960146655688,
                            0.1742608664057115,
                        ),
                    ),
                    0b0101 => (
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            0.4271821418762962,
                            0.08806052015143966,
                            -0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            -0.4271821418762962,
                            -0.08806052015143966,
                            0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            -0.3662936517508201,
                            0.1539383024749358,
                            -0.1991552758922418,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            0.3662936517508201,
                            -0.1539383024749358,
                            0.1991552758922418,
                        ),
                    ),
                    0b1101 => (
                        LorentzVector::from_args(
                            -0.4375000000000000,
                            0.4271821418762962,
                            0.08806052015143966,
                            -0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            -0.5625000000000000,
                            -0.4271821418762962,
                            -0.08806052015143966,
                            0.03416156396650795,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            -0.4120803582196727,
                            0.1731805902843028,
                            -0.2240496853787720,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            0.4120803582196727,
                            -0.1731805902843028,
                            0.2240496853787720,
                        ),
                    ),
                    0b0011 => (
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            0.4882081621443385,
                            0.1006405944587882,
                            -0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            -0.4882081621443385,
                            -0.1006405944587882,
                            0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            -0.2422804773385922,
                            0.1018206163990297,
                            -0.1317288330743322,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            0.2422804773385922,
                            -0.1018206163990297,
                            0.1317288330743322,
                        ),
                    ),
                    0b1011 => (
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            0.4882081621443385,
                            0.1006405944587882,
                            -0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            -0.4882081621443385,
                            -0.1006405944587882,
                            0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            0.6111111111111111,
                            -0.3205069452819676,
                            0.1346960146655688,
                            -0.1742608664057115,
                        ),
                        LorentzVector::from_args(
                            0.3888888888888889,
                            0.3205069452819676,
                            -0.1346960146655688,
                            0.1742608664057115,
                        ),
                    ),
                    0b0111 => (
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            0.4882081621443385,
                            0.1006405944587882,
                            -0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            -0.4882081621443385,
                            -0.1006405944587882,
                            0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            0.4444444444444444,
                            -0.3662936517508201,
                            0.1539383024749358,
                            -0.1991552758922418,
                        ),
                        LorentzVector::from_args(
                            0.5555555555555556,
                            0.3662936517508201,
                            -0.1539383024749358,
                            0.1991552758922418,
                        ),
                    ),
                    0b1111 => (
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            0.4882081621443385,
                            0.1006405944587882,
                            -0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            -0.5000000000000000,
                            -0.4882081621443385,
                            -0.1006405944587882,
                            0.03904178739029480,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            -0.4120803582196727,
                            0.1731805902843028,
                            -0.2240496853787720,
                        ),
                        LorentzVector::from_args(
                            0.5000000000000000,
                            0.4120803582196727,
                            -0.1731805902843028,
                            0.2240496853787720,
                        ),
                    ),
                    _ => unimplemented!("Invalid box PS {:.4b}", on_shell_flag),
                };

                let m = 0.;

                println!(
                    "p1^2: {}, p2^2: {}, p3^3: {}",
                    p1.square(),
                    p2.square(),
                    p3.square()
                );
                let e_cm = (p1 + p2).square().abs().sqrt();
                let ll = build_loop_line(
                    &[(0, true)],
                    vec![p1, p2, p3, -p1 - p2 - p3],
                    vec![m, m, m, m],
                    true,
                );
                (1, e_cm, vec![ll])
            }
            "double-triangle" => {
                let p = LorentzVector::from_args(1.0, 0., 0., 0.);
                let m = 0.;
                let e_cm = p.square().abs().sqrt();

                let ll1 = build_loop_line(
                    &[(0, true)],
                    vec![LorentzVector::default(), -p],
                    vec![m, m],
                    false,
                );
                let ll2 = build_loop_line(
                    &[(1, true)],
                    vec![LorentzVector::default(), p],
                    vec![m, m],
                    false,
                );
                let ll3 = build_loop_line(
                    &[(0, true), (1, true)],
                    vec![LorentzVector::default()],
                    vec![m],
                    false,
                );
                (2, e_cm, vec![ll1, ll2, ll3])
            }
            "triangle-box" => {
                let p1 = LorentzVector::from_args(1.0, 0.0, 0.0, 0.0);
                let p2 = LorentzVector::from_args(-0.59375, -0.0, -0.0, -0.32021721143623744);
                let _p3 = LorentzVector::from_args(-0.40625, -0.0, -0.0, 0.32021721143623744);
                let m = 0.;
                let e_cm = 1.0;

                let ll1 = build_loop_line(
                    &[(0, true)],
                    vec![LorentzVector::default(), -p1],
                    vec![m, m],
                    false,
                );
                let ll2 = build_loop_line(
                    &[(1, true)],
                    vec![LorentzVector::default(), -p2, p1],
                    vec![m, m, m],
                    false,
                );
                let ll3 = build_loop_line(&[(0, true), (1, true)], vec![-p1], vec![m], false);
                (2, e_cm, vec![ll1, ll2, ll3])
            }
            x => unimplemented!("Unknown topology {}", x),
        }
    }*/
}
