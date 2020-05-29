use arrayvec::ArrayVec;
use color_eyre::{Help, Report};
use eyre::WrapErr;
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
use rand::Rng;
use serde::Deserialize;
use std::fs::File;
use std::mem;
use std::ops::Mul;
use topologies::{Cut, LTDCache, LTDNumerator, Topology};
use utils;
use utils::{PolynomialReconstruction, Signum};
use vector::{LorentzVector, RealNumberLike};
use {DeformationStrategy, FloatLike, NormalisingFunction, Settings, MAX_LOOP};

#[cfg(feature = "mg_numerator")]
use dlopen::wrapper::Container;

#[cfg(feature = "mg_numerator")]
mod MGNumeratorMod {
    use dlopen::wrapper::{Container, WrapperApi};
    use libc::{c_char, c_double, c_int};
    use num::Complex;

    #[derive(WrapperApi)]
    pub struct MGNumeratorAPI {
        c_get_numerator: unsafe extern "C" fn(
            p: *const c_double,
            proc_id: *const c_int,
            selected_diagam_left: *const c_int,
            selected_diagram_right: *const c_int,
            ans_re: *mut c_double,
            ans_im: *mut c_double,
        ),
        c_initialise: unsafe extern "C" fn(p: *const c_char),
    }

    pub fn initialise(api_container: &mut Container<MGNumeratorAPI>, param_card_filename: &str) {
        unsafe {
            let mut a = ['\0' as i8; 512];
            for (xa, c) in a.iter_mut().zip(param_card_filename.chars()) {
                *xa = c as i8;
            }
            api_container.c_initialise(a.as_ptr());
        }
    }

    pub fn get_numerator(
        api_container: &mut Container<MGNumeratorAPI>,
        p: &[f64],
        proc_id: usize,
        left_diagram_id: usize,
        right_diagram_id: usize,
    ) -> Complex<f64> {
        let proc_id = proc_id as i32;
        let selected_diagram_left = left_diagram_id as i32;
        let selected_diagram_right = right_diagram_id as i32;

        unsafe {
            let mut ans = Complex::default();
            api_container.c_get_numerator(
                &p[0] as *const f64,
                &proc_id as *const i32,
                &selected_diagram_left as *const i32,
                &selected_diagram_right as *const i32,
                &mut ans.re as *mut f64,
                &mut ans.im as *mut f64,
            );
            ans
        }
    }

    pub fn load() -> Container<MGNumeratorAPI> {
        let mut path = std::env::var("MG_NUMERATOR_PATH")
            .expect("MG_NUMERATOR_PATH needs to be set in the mg_numerator mode.");

        let mut container: Container<MGNumeratorAPI> =
            unsafe { Container::load(&(path.clone() + "lib/libMGnumerators_dynamic.so")) }
                .expect("Could not open library or load symbols");

        path += "/Cards/param_card.dat";
        initialise(&mut container, &path);

        container
    }
}

#[cfg(feature = "mg_numerator")]
mod FORMNumeratorMod {
    use dlopen::wrapper::{Container, WrapperApi};
    use libc::{c_double, c_int};
    use num::Complex;

    #[derive(WrapperApi)]
    pub struct FORMNumeratorAPI {
        evaluate: unsafe extern "C" fn(p: *const c_double, diag: c_int, conf: c_int, out: *mut c_double),
    }

    pub fn get_numerator(
        api_container: &mut Container<FORMNumeratorAPI>,
        p: &[f64],
        diag: usize,
        conf: usize,
        poly: &mut[f64]
    ) {
        unsafe {
            api_container.evaluate(&p[0] as *const f64, diag as i32, conf as i32, &mut poly[0] as *mut f64);
        }
    }

    pub fn load() -> Container<FORMNumeratorAPI> {
        let path = std::env::var("MG_NUMERATOR_PATH")
            .expect("MG_NUMERATOR_PATH needs to be set in the mg_numerator mode.");

        let container: Container<FORMNumeratorAPI> =
            unsafe { Container::load(&(path.clone() + "lib/libFORM_numerators.so")) }
                .expect("Could not open library or load symbols");

        container
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCut {
    pub name: String,
    pub sign: i8,
    pub level: usize,
    pub particle_id: i8,
    pub signature: (Vec<i8>, Vec<i8>),
    pub m_squared: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCutLimits {
    pub diagrams: Vec<Topology>,
    pub conjugate_deformation: Vec<bool>,
    #[serde(default)]
    pub numerator_tensor_coefficients_sparse: Vec<(Vec<usize>, (f64, f64))>,
    #[serde(default)]
    pub numerator_tensor_coefficients: Vec<(f64, f64)>,
    #[serde(skip_deserializing)]
    pub numerator: LTDNumerator,
    pub cb_to_lmb: Option<Vec<i8>>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub energy_polynomial_matrix: Vec<Vec<f64>>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub energy_polynomial_map: Vec<usize>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub energy_polynomial_samples: Vec<Vec<f64>>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub energy_polynomial_coefficients: Vec<Complex<f64>>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCuts {
    pub cuts: Vec<CutkoskyCut>,
    pub n_bubbles: usize,
    pub uv_limits: Vec<CutkoskyCutLimits>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MGNumeratorCallSignature {
    pub left_diagram_id: usize,
    pub proc_id: usize,
    pub right_diagram_id: usize,
}

#[derive(Debug, Clone, Deserialize)]
pub struct FORMNumeratorCallSignature {
    pub id: usize,
}

#[derive(Deserialize)]
pub struct MGNumerator {
    call_signature: Option<MGNumeratorCallSignature>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub mg_numerator: Option<Container<MGNumeratorMod::MGNumeratorAPI>>,
}

impl Clone for MGNumerator {
    fn clone(&self) -> Self {
        MGNumerator {
            call_signature: self.call_signature.clone(),
            #[cfg(feature = "mg_numerator")]
            mg_numerator: if self.call_signature.is_some() {
                Some(MGNumeratorMod::load())
            } else {
                None
            },
        }
    }
}

#[derive(Deserialize)]
pub struct FORMNumerator {
    call_signature: Option<FORMNumeratorCallSignature>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub form_numerator: Option<Container<FORMNumeratorMod::FORMNumeratorAPI>>,
    #[serde(skip_deserializing)]
    #[cfg(feature = "mg_numerator")]
    pub form_numerator_buffer: Vec<f64>,
}

impl Clone for FORMNumerator {
    fn clone(&self) -> Self {
        FORMNumerator {
            call_signature: self.call_signature.clone(),
            #[cfg(feature = "mg_numerator")]
            form_numerator: if self.call_signature.is_some() {
                Some(FORMNumeratorMod::load())
            } else {
                None
            },
            form_numerator_buffer: self.form_numerator_buffer.clone(),
        }
    }
}

#[derive(Clone, Deserialize)]
pub struct SquaredTopology {
    pub name: String,
    pub n_loops: usize,
    pub n_incoming_momenta: usize,
    pub e_cm_squared: f64,
    pub overall_numerator: f64,
    pub numerator_in_loop_momentum_basis: bool,
    pub external_momenta: Vec<LorentzVector<f64>>,
    pub cutkosky_cuts: Vec<CutkoskyCuts>,
    #[serde(skip_deserializing)]
    pub settings: Settings,
    #[serde(skip_deserializing)]
    pub rotation_matrix: [[float; 3]; 3],
    pub topo: Topology,
    #[serde(rename = "MG_numerator")]
    pub mg_numerator: MGNumerator,
    #[serde(rename = "FORM_numerator")]
    pub form_numerator: FORMNumerator,
}

#[derive(Clone)]
pub struct SquaredTopologySet {
    pub name: String,
    pub n_loops: usize,
    pub e_cm_squared: f64,
    topologies: Vec<SquaredTopology>,
    pub multiplicity: Vec<usize>,
    pub settings: Settings,
    pub rotation_matrix: [[float; 3]; 3],
    pub multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f64>>)>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetTopology {
    pub name: String,
    pub multiplicity: usize,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetInput {
    pub name: String,
    topologies: Vec<SquaredTopologySetTopology>,
}

impl SquaredTopologySet {
    pub fn from_one(mut squared_topology: SquaredTopology) -> SquaredTopologySet {
        let channels = squared_topology.generate_multi_channeling_channels();

        SquaredTopologySet {
            name: squared_topology.name.clone(),
            n_loops: squared_topology.n_loops,
            settings: squared_topology.settings.clone(),
            e_cm_squared: squared_topology.e_cm_squared,
            rotation_matrix: squared_topology.rotation_matrix.clone(),
            topologies: vec![squared_topology],
            multiplicity: vec![1],
            multi_channeling_channels: channels,
        }
    }

    pub fn from_file(filename: &str, settings: &Settings) -> Result<SquaredTopologySet, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open squared topology set file {}", filename))
            .suggestion("Does the path exist?")?;

        let squared_topology_set_input: SquaredTopologySetInput = serde_yaml::from_reader(f)
            .wrap_err("Could not parse squared topology set file")
            .suggestion("Is it a correct yaml file")?;

        let mut topologies: Vec<SquaredTopology> = vec![];
        let mut multiplicity: Vec<usize> = vec![];
        let mut multi_channeling_channels = vec![];

        for topo in squared_topology_set_input.topologies {
            let filename = std::path::Path::new(&filename)
                .with_file_name(topo.name)
                .with_extension("yaml");
            let mut squared_topology =
                SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                    .wrap_err("Could not load subtopology file")?;

            if !topologies.is_empty() && squared_topology.n_loops != topologies[0].n_loops {
                panic!("Topology sets require all topologies to have the same number of loops");
            }

            multi_channeling_channels.extend(squared_topology.generate_multi_channeling_channels());
            topologies.push(squared_topology);
            multiplicity.push(topo.multiplicity);
        }

        // filter duplicate channels between supergraphs
        let mut unique_multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f64>>)> =
            vec![];

        'mc_loop: for c in multi_channeling_channels {
            for uc in &unique_multi_channeling_channels {
                if uc.0 == c.0 && uc.1 == c.1 {
                    // test for true equality
                    if uc
                        .2
                        .iter()
                        .zip(c.2.iter())
                        .all(|(s1, s2)| (s1 - s2).spatial_squared() < 1.0e-15)
                    {
                        continue 'mc_loop;
                    }
                }
            }
            unique_multi_channeling_channels.push(c);
        }

        let rotation_matrix = [
            [float::one(), float::zero(), float::zero()],
            [float::zero(), float::one(), float::zero()],
            [float::zero(), float::zero(), float::one()],
        ];

        Ok(SquaredTopologySet {
            name: squared_topology_set_input.name,
            n_loops: topologies[0].n_loops,
            e_cm_squared: topologies[0].e_cm_squared,
            topologies,
            rotation_matrix,
            settings: settings.clone(),
            multiplicity,
            multi_channeling_channels: unique_multi_channeling_channels,
        })
    }

    pub fn create_caches<T: FloatLike>(&self) -> Vec<Vec<Vec<Vec<LTDCache<T>>>>> {
        self.topologies.iter().map(|t| t.create_caches()).collect()
    }

    pub fn multi_channeling<'a, T: FloatLike>(
        &mut self,
        x: &'a [f64],
        cache: &mut [Vec<Vec<Vec<LTDCache<T>>>>],
        mut event_manager: Option<&mut EventManager>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        T,
        Complex<T>,
        Complex<T>,
    ) {
        // paramaterize and consider the result in a channel basis
        let n_loops = x.len() / 3;
        let mut k_channel = [LorentzVector::default(); MAX_LOOP];
        for i in 0..n_loops {
            let (l_space, _) = Topology::parameterize::<T>(
                &x[i * 3..(i + 1) * 3],
                self.e_cm_squared,
                i,
                &self.settings,
            );

            let rot = &self.rotation_matrix;
            k_channel[i] = LorentzVector::from_args(
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
        }

        let mut k_lmb = [LorentzVector::default(); MAX_LOOP];
        let mut k_other_channel = [LorentzVector::default(); MAX_LOOP];
        let mut event_counter = 0;
        let mut result = Complex::zero();
        for (channel, _, channel_shift) in &self.multi_channeling_channels {
            // transform to the loop momentum basis
            for (kk, r) in k_lmb[..self.n_loops]
                .iter_mut()
                .zip_eq(channel.chunks(self.n_loops))
            {
                *kk = LorentzVector::default();

                for ((ko, s), shift) in k_channel[..self.n_loops]
                    .iter()
                    .zip_eq(r.iter())
                    .zip_eq(channel_shift)
                {
                    *kk += (ko - shift.cast()).multiply_sign(*s);
                }
            }

            // determine the normalization constant
            let mut normalization = T::zero();
            for (_, other_channel_inv, other_channel_shift) in &self.multi_channeling_channels {
                // transform from the loop momentum basis to the other channel basis
                for ((kk, r), shift) in k_other_channel[..self.n_loops]
                    .iter_mut()
                    .zip_eq(other_channel_inv.chunks(self.n_loops))
                    .zip_eq(other_channel_shift)
                {
                    *kk = shift.cast();
                    for (ko, s) in k_lmb[..self.n_loops].iter().zip_eq(r.iter()) {
                        *kk += ko.multiply_sign(*s);
                    }
                }

                let mut inv_jac_para = T::one();
                for i in 0..n_loops {
                    let rot = &self.rotation_matrix;

                    // undo the rotation
                    let rotated = LorentzVector::from_args(
                        T::zero(),
                        <T as NumCast>::from(rot[0][0]).unwrap() * k_other_channel[i].x
                            + <T as NumCast>::from(rot[1][0]).unwrap() * k_other_channel[i].y
                            + <T as NumCast>::from(rot[2][0]).unwrap() * k_other_channel[i].z,
                        <T as NumCast>::from(rot[0][1]).unwrap() * k_other_channel[i].x
                            + <T as NumCast>::from(rot[1][1]).unwrap() * k_other_channel[i].y
                            + <T as NumCast>::from(rot[2][1]).unwrap() * k_other_channel[i].z,
                        <T as NumCast>::from(rot[0][2]).unwrap() * k_other_channel[i].x
                            + <T as NumCast>::from(rot[1][2]).unwrap() * k_other_channel[i].y
                            + <T as NumCast>::from(rot[2][2]).unwrap() * k_other_channel[i].z,
                    );
                    let (_, jac) =
                        Topology::inv_parametrize(&rotated, self.e_cm_squared, i, &self.settings);

                    inv_jac_para *= jac;
                }

                normalization += inv_jac_para;
            }

            // evaluate the integrand in this channel
            let mut channel_result = Complex::zero();
            for ((t, cache), &m) in self
                .topologies
                .iter_mut()
                .zip_eq(cache.iter_mut())
                .zip_eq(&self.multiplicity)
            {
                channel_result += t.evaluate_mom(&k_lmb[..n_loops], cache, &mut event_manager)
                    / normalization
                    * Into::<T>::into(m as f64);
                if let Some(em) = &mut event_manager {
                    if em.track_events {
                        for e in em.event_buffer[event_counter..].iter_mut() {
                            e.integrand *= m as f64 / normalization.to_f64().unwrap();
                        }
                        event_counter = em.event_buffer.len();
                    }
                }
            }

            result += channel_result;
        }

        // we don't have a meaningful jacobian nor k_def
        let k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]> = (0..self.n_loops)
            .map(|_| LorentzVector::default())
            .collect();
        (x, k_def, T::one(), Complex::one(), result)
    }

    pub fn evaluate<'a, T: FloatLike>(
        &mut self,
        x: &'a [f64],
        cache: &mut [Vec<Vec<Vec<LTDCache<T>>>>],
        mut event_manager: Option<&mut EventManager>,
    ) -> (
        &'a [f64],
        ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]>,
        T,
        Complex<T>,
        Complex<T>,
    ) {
        if self.settings.general.multi_channeling && !self.multi_channeling_channels.is_empty() {
            return self.multi_channeling(x, cache, event_manager);
        }

        // jointly parameterize all squared topologies
        let n_loops = x.len() / 3;
        let mut k = [LorentzVector::default(); MAX_LOOP];
        let mut jac_para = T::one();
        for i in 0..n_loops {
            // set the loop index to i + 1 so that we can also shift k
            let (l_space, jac) = Topology::parameterize(
                &x[i * 3..(i + 1) * 3],
                self.e_cm_squared, // NOTE: taking e_cm from the first graph
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

        let mut result = Complex::zero();
        let mut event_counter = 0;
        for ((t, cache), &m) in self
            .topologies
            .iter_mut()
            .zip_eq(cache.iter_mut())
            .zip_eq(&self.multiplicity)
        {
            result += t.evaluate_mom(&k[..n_loops], cache, &mut event_manager)
                * jac_para
                * Into::<T>::into(m as f64);

            if let Some(em) = &mut event_manager {
                if em.track_events {
                    for e in em.event_buffer[event_counter..].iter_mut() {
                        e.integrand *= jac_para.to_f64().unwrap() * m as f64;
                    }
                    event_counter = em.event_buffer.len();
                }
            }
        }

        // NOTE: there is no unique k_def anymore. it depends on the cut
        let mut k_def: ArrayVec<[LorentzVector<Complex<T>>; MAX_LOOP]> = ArrayVec::default();
        for l in &k[..n_loops] {
            k_def.push(l.map(|x| Complex::new(x, T::zero())));
        }

        (x, k_def, jac_para, Complex::one(), result)
    }
}

/// Cache for squared topology sets.
#[derive(Default)]
pub struct SquaredTopologyCache {
    float_cache: Vec<Vec<Vec<Vec<LTDCache<float>>>>>,
    quad_cache: Vec<Vec<Vec<Vec<LTDCache<f128>>>>>,
}

pub trait CachePrecisionSelector<T: FloatLike> {
    fn get(&mut self) -> &mut Vec<Vec<Vec<Vec<LTDCache<T>>>>>;
}

impl CachePrecisionSelector<float> for SquaredTopologyCache {
    #[inline]
    fn get(&mut self) -> &mut Vec<Vec<Vec<Vec<LTDCache<float>>>>> {
        &mut self.float_cache
    }
}

impl CachePrecisionSelector<f128> for SquaredTopologyCache {
    #[inline]
    fn get(&mut self) -> &mut Vec<Vec<Vec<Vec<LTDCache<f128>>>>> {
        &mut self.quad_cache
    }
}

impl SquaredTopology {
    pub fn from_file(filename: &str, settings: &Settings) -> Result<SquaredTopology, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open squared topology file {}", filename))
            .suggestion("Does the path exist?")?;

        let mut squared_topo: SquaredTopology = serde_yaml::from_reader(f)
            .wrap_err("Could not parse squared topology file")
            .suggestion("Is it a correct yaml file")?;

        squared_topo.settings = settings.clone();
        for cutkosky_cuts in &mut squared_topo.cutkosky_cuts {
            for uv_limit in &mut cutkosky_cuts.uv_limits {
                uv_limit.numerator = if uv_limit.numerator_tensor_coefficients.len() == 0
                    && uv_limit.numerator_tensor_coefficients_sparse.len() == 0
                {
                    LTDNumerator::one(squared_topo.n_loops + cutkosky_cuts.n_bubbles)
                } else {
                    if !uv_limit.numerator_tensor_coefficients_sparse.is_empty() {
                        LTDNumerator::from_sparse(
                            squared_topo.n_loops + cutkosky_cuts.n_bubbles,
                            &uv_limit.numerator_tensor_coefficients_sparse,
                        )
                    } else {
                        LTDNumerator::new(
                            squared_topo.n_loops + cutkosky_cuts.n_bubbles,
                            &uv_limit
                                .numerator_tensor_coefficients
                                .iter()
                                .map(|x| Complex::new(x.0, x.1))
                                .collect::<Vec<_>>(),
                        )
                    }
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

        #[cfg(feature = "mg_numerator")]
        {
            if squared_topo.mg_numerator.call_signature.is_some() {
                squared_topo.mg_numerator.mg_numerator = Some(MGNumeratorMod::load());
                squared_topo.get_energy_polynomial_matrix();
            }

            if squared_topo.form_numerator.call_signature.is_some() {
                squared_topo.form_numerator.form_numerator = Some(FORMNumeratorMod::load());
            }
        }

        Ok(squared_topo)
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

    /// Reconstruct the shape of the energy polynomial for the MG numerator and
    /// construct the transformation matrix that finds the coefficients.
    #[cfg(feature = "mg_numerator")]
    pub fn get_energy_polynomial_matrix(&mut self) {
        let mut rng = rand::thread_rng();

        let e_cm = self.e_cm_squared.sqrt();
        let max_rank: usize = 8; // FIXME: make dynamic
        let mut poly_rec = PolynomialReconstruction::new(max_rank, self.n_loops, e_cm);

        let mut mg_numerator = mem::replace(&mut self.mg_numerator.mg_numerator, None);
        let mut poly_rec_buffer = vec![Complex::zero(); max_rank.pow(self.n_loops as u32)];

        for cutkosky_cuts in &mut self.cutkosky_cuts {
            let mut shifts = [LorentzVector::default(); MAX_LOOP + 4];
            // get the shifts
            for (shift, cut) in shifts.iter_mut().zip(cutkosky_cuts.cuts.iter()) {
                *shift = utils::evaluate_signature(&cut.signature.1, &self.external_momenta);
            }

            for cut_uv_limit in cutkosky_cuts.uv_limits.iter_mut() {
                if let Some(call_signature) = &self.mg_numerator.call_signature {
                    let mut all_external_momenta: ArrayVec<[f64; (MAX_LOOP + 4) * 8]> =
                        ArrayVec::new();

                    let cut_energies: Vec<f64> =
                        (0..self.n_loops).map(|_| rng.gen::<f64>() * e_cm).collect();

                    for e in &self.external_momenta[..self.n_incoming_momenta] {
                        all_external_momenta
                            .try_extend_from_slice(&[e.t, 0., e.x, 0., e.y, 0., e.z, 0.])
                            .unwrap();
                    }

                    for _ in 0..self.n_loops {
                        all_external_momenta
                            .try_extend_from_slice(&[
                                rng.gen::<f64>() * e_cm,
                                0.,
                                rng.gen::<f64>() * e_cm,
                                0.,
                                rng.gen::<f64>() * e_cm,
                                0.,
                                rng.gen::<f64>() * e_cm,
                                0.,
                            ])
                            .unwrap();
                    }

                    let n_incoming_momenta = self.n_incoming_momenta;
                    let n_loops = self.n_loops;
                    let ltd_start_index = cutkosky_cuts.cuts.len() - 1;
                    let mut f = |x: &[f64]| {
                        // set the energies of the LTD momenta
                        if let Some(m) = &cut_uv_limit.cb_to_lmb {
                            for (i, r) in m.chunks(n_loops).enumerate() {
                                all_external_momenta[(n_incoming_momenta + i) * 8] = 0.;
                                all_external_momenta[(n_incoming_momenta + i) * 8 + 1] = 0.;
                                // x is only for the LTD momenta, so we need to
                                for (j, ((&sign, cut_e), shift)) in r
                                    .iter()
                                    .zip_eq(cut_energies.iter())
                                    .zip_eq(&shifts[..n_loops])
                                    .enumerate()
                                {
                                    if sign != 0 {
                                        if j < ltd_start_index {
                                            all_external_momenta[(n_incoming_momenta + i) * 8] +=
                                                (*cut_e - shift.t).multiply_sign(sign)
                                        } else {
                                            all_external_momenta[(n_incoming_momenta + i) * 8] +=
                                                (x[j - ltd_start_index] - shift.t)
                                                    .multiply_sign(sign);
                                        }
                                    }
                                }
                            }
                        }

                        MGNumeratorMod::get_numerator(
                            mg_numerator.as_mut().unwrap(),
                            &all_external_momenta,
                            call_signature.proc_id,
                            call_signature.left_diagram_id,
                            call_signature.right_diagram_id,
                        )
                    };

                    let n_vars = self.n_loops + 1 - cutkosky_cuts.cuts.len();

                    poly_rec_buffer.clear();
                    poly_rec_buffer.resize(max_rank.pow(n_vars as u32), Complex::zero());

                    poly_rec.reconstruct(&mut f, max_rank, n_vars, &mut poly_rec_buffer);

                    if self.settings.general.debug > 1 {
                        println!("Fit polynomial coefficients: {:?}", poly_rec_buffer);
                    }

                    let mut non_zero_coefficients = vec![];
                    let mut max_rank_non_zero = 0;

                    for (i, v) in poly_rec_buffer.iter().enumerate() {
                        // TODO: check if these bounds are appropriate
                        if v.norm_sqr() > self.e_cm_squared * 1e-6 {
                            non_zero_coefficients.push(i);

                            for j in 0..n_vars {
                                let pow = (i / max_rank.pow((n_vars - 1 - j) as u32)) % max_rank;
                                if pow > max_rank_non_zero {
                                    max_rank_non_zero = pow;
                                }
                            }
                        }
                    }

                    if self.settings.general.debug > 0 {
                        println!(
                            "Numerator shape: {:?} with max rank {}",
                            non_zero_coefficients, max_rank_non_zero
                        );
                    }

                    cut_uv_limit.numerator = LTDNumerator::from_sparse(
                        n_loops,
                        &[(vec![0; max_rank_non_zero], (0., 0.))],
                    );

                    // also give the subgraph numerators the proper size
                    // TODO: set the proper rank of only the variables of the subgraph
                    for subgraph in &mut cut_uv_limit.diagrams {
                        if subgraph.n_loops > 0 {
                            subgraph.numerator = LTDNumerator::from_sparse(
                                subgraph.n_loops,
                                &[(vec![0; max_rank_non_zero], (0., 0.))],
                            );
                        }
                    }

                    // collect all non-zero components and construct a new matrix for this particular shape
                    let mut mat = vec![];
                    let mut all_sample_points = vec![];
                    let mut numerator_to_reduced_index_map = vec![];

                    cut_uv_limit.energy_polynomial_coefficients.clear();
                    cut_uv_limit
                        .energy_polynomial_coefficients
                        .resize(non_zero_coefficients.len(), Complex::zero());

                    // construct the map to the reduced index
                    for i in &non_zero_coefficients {
                        let mut power_map = vec![];
                        for j in 0..n_vars {
                            let pow = (i / max_rank.pow((n_vars - 1 - j) as u32)) % max_rank;
                            for _ in 0..pow {
                                power_map.push((cutkosky_cuts.cuts.len() - 1 + j) * 4);
                            }
                        }

                        if power_map.len() == 0 {
                            numerator_to_reduced_index_map.push(0);
                        } else {
                            let index_full = cut_uv_limit.numerator.sorted_linear[..]
                                .binary_search_by(|x| utils::compare_slice(&x[..], &power_map[..]))
                                .unwrap()
                                + 1;

                            let index = cut_uv_limit.numerator.coeff_map_to_reduced_numerator
                                [cutkosky_cuts.cuts.len() - 1][index_full];
                            numerator_to_reduced_index_map.push(index);
                        }
                    }

                    let total_samples = if max_rank_non_zero > 0 {
                        non_zero_coefficients.len() + 2 // take 2 more samples
                    } else {
                        non_zero_coefficients.len()
                    };

                    for _ in 0..total_samples {
                        let sample_points: Vec<f64> =
                            (0..n_vars).map(|_| rng.gen::<f64>() * e_cm).collect();
                        let mut row = vec![];

                        for i in &non_zero_coefficients {
                            let eval: f64 = sample_points
                                .iter()
                                .enumerate()
                                .map(|(j, sp)| {
                                    sp.powi(
                                        ((i / max_rank.pow((n_vars - 1 - j) as u32)) % max_rank)
                                            as i32,
                                    )
                                })
                                .product();
                            row.push(eval);
                        }

                        all_sample_points.push(sample_points);
                        mat.push(row);
                    }

                    let matl: Vec<f64> = mat.iter().flatten().cloned().collect();
                    let m = nalgebra::DMatrix::from_row_slice(
                        total_samples,
                        non_zero_coefficients.len(),
                        &matl,
                    );

                    let inv = (m.transpose().mul(m.clone()))
                        .try_inverse()
                        .unwrap()
                        .mul(m.transpose());

                    // TODO: test stability by sampling points

                    cut_uv_limit.energy_polynomial_matrix = inv
                        .row_iter()
                        .map(|r| r.iter().cloned().collect::<Vec<_>>())
                        .collect();
                    cut_uv_limit.energy_polynomial_map = numerator_to_reduced_index_map;
                    cut_uv_limit.energy_polynomial_samples = all_sample_points;
                }
            }
        }

        mem::swap(&mut mg_numerator, &mut self.mg_numerator.mg_numerator);
    }

    pub fn generate_multi_channeling_channels(
        &mut self,
    ) -> Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f64>>)> {
        // construct the channels from all loop momentum basis (i.e. the LTD cuts)
        self.topo.process(true);

        let mut multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f64>>)> =
            vec![];
        for (cuts, mat) in self
            .topo
            .ltd_cut_options
            .iter()
            .zip_eq(self.topo.cb_to_lmb_mat.iter())
        {
            for cut in cuts {
                let mut lmb_to_cb_mat = vec![];
                let mut cut_shifts = vec![];

                let mut include = true;
                for (ll_cut_sig, ll_cut) in cut.iter().zip_eq(&self.topo.loop_lines) {
                    if let Cut::PositiveCut(i) | Cut::NegativeCut(i) = ll_cut_sig {
                        // optionally leave out massive propagators for efficiency
                        if !self
                            .settings
                            .general
                            .multi_channeling_including_massive_propagators
                            && ll_cut.propagators[*i].m_squared > 0.
                        {
                            include = false;
                            break;
                        }
                        cut_shifts.push(ll_cut.propagators[*i].q);
                        lmb_to_cb_mat.extend(&ll_cut.signature);
                    }
                }

                if include {
                    // check if we have seen the channel before
                    // this is possible since we only look at the spatial parts of the shift
                    for uc in &multi_channeling_channels {
                        if &uc.0 == mat && uc.1 == lmb_to_cb_mat {
                            // test for true equality
                            if uc
                                .2
                                .iter()
                                .zip(cut_shifts.iter())
                                .all(|(s1, s2)| (s1 - s2).spatial_squared() < 1.0e-15)
                            {
                                include = false;
                            }
                        }
                    }

                    if include {
                        multi_channeling_channels.push((mat.clone(), lmb_to_cb_mat, cut_shifts));
                    }
                }
            }
        }
        multi_channeling_channels
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

        // find the index of the current event
        let event_counter = if let Some(em) = event_manager {
            em.event_buffer.len()
        } else {
            0
        };

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
            let n_bubbles = cutkosky_cuts.n_bubbles;

            if self.settings.general.debug >= 1 {
                println!(
                    "Cut {}:",
                    cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                );
            }

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
                    &mut k_def[..self.n_loops + n_bubbles],
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
                for e in em.event_buffer[event_counter..].iter_mut() {
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

        let mut k_def_lmb = [LorentzVector::<Complex<T>>::default(); MAX_LOOP + 4];
        let mut shifts = [LorentzVector::default(); MAX_LOOP + 4];

        // evaluate the cuts with the proper scaling
        for ((cut_mom, shift), cut) in cut_momenta[..cutkosky_cuts.cuts.len()]
            .iter_mut()
            .zip(shifts.iter_mut())
            .zip(cutkosky_cuts.cuts.iter())
        {
            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            *shift = utils::evaluate_signature(
                &cut.signature.1,
                &external_momenta[..self.external_momenta.len()],
            );
            *cut_mom = k * scaling + *shift;
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
            let (h, h_norm) = match self.settings.cross_section.normalising_function.name {
                NormalisingFunction::RightExponential => {
                    // Only support center at 1 for now
                    assert_eq!(self.settings.cross_section.normalising_function.center, 1.,
                        "For now the right exponential normalising function only support centers at 1.0.");
                    (
                        (-Float::powi(
                            scaling
                                / Into::<T>::into(
                                    self.settings.cross_section.normalising_function.spread,
                                ),
                            2,
                        ))
                        .exp(),
                        (<T as FloatConst>::PI().sqrt() / Into::<T>::into(2.0)
                            * Into::<T>::into(
                                self.settings.cross_section.normalising_function.spread,
                            )),
                    )
                }
                NormalisingFunction::LeftRightExponential => {
                    // Only support center and spread at 1 for now
                    assert_eq!(self.settings.cross_section.normalising_function.center, 1.,
                        "For now the left-right exponential normalising function only support a center at 1.0.");
                    assert_eq!(self.settings.cross_section.normalising_function.spread, 1.,
                        "For now the left-right exponential normalising function only support a spread set to 1.0.");
                    (
                        (-((Float::powi(scaling, 2)
                            + Float::powi(
                                Into::<T>::into(
                                    self.settings.cross_section.normalising_function.spread,
                                ),
                                2,
                            ))
                            / scaling))
                            .exp(),
                        // 2 Sqrt[\[Sigma]] BesselK[1, 2 Sqrt[\[Sigma]]], with \[Sigma]=1
                        Into::<T>::into(0.27973176363304485456919761407082204777),
                    )
                }
                NormalisingFunction::None => (Into::<T>::into(1.0), Into::<T>::into(1.0)),
            };
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

        if self.settings.general.debug >= 2 {
            println!(
                "  | rescaled loop momenta = {:?}",
                &rescaled_loop_momenta[..self.n_loops]
            );
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

        // for the evaluation of the numerator we need complex loop momenta of the supergraph.
        // the order is: cut momenta, momenta left graph, momenta right graph
        // NOTE: this requires that the last level 0 cutkosky cut is at the end of the list
        for (kd, cut_mom) in k_def
            .iter_mut()
            .zip(&cut_momenta[..cutkosky_cuts.cuts.len() - 1])
        {
            *kd = cut_mom.to_complex(true);
        }

        // now apply the same procedure for all uv limits
        let mut diag_and_num_contributions = Complex::zero();
        let mut def_jacobian = Complex::one();
        // regenerate the evaluation of the exponent map of the numerator since the loop momenta have changed
        let mut regenerate_momenta = true;
        for (uv_index, (cut_uv_limit, diag_cache)) in
            cutkosky_cuts.uv_limits.iter_mut().zip(cache).enumerate()
        {
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
                    if uv_index == 0 {
                        subgraph.fixed_deformation =
                            subgraph.determine_ellipsoid_overlap_structure(true);
                    } else {
                        subgraph.fixed_deformation = vec![];
                    }

                    if self.settings.general.debug > 0 {
                        // check if the overlap structure makes sense
                        subgraph.check_fixed_deformation();
                    }
                }
            }

            if !self
                .settings
                .cross_section
                .inherit_deformation_for_uv_counterterm
                || uv_index == 0
            {
                def_jacobian = Complex::one();
            }

            // compute the deformation vectors
            let mut k_def_index = cutkosky_cuts.cuts.len() - 1;
            for ((subgraph, subgraph_cache), &conjugate_deformation) in cut_uv_limit
                .diagrams
                .iter_mut()
                .zip_eq(diag_cache.iter_mut())
                .zip_eq(cut_uv_limit.conjugate_deformation.iter())
            {
                if !self
                    .settings
                    .cross_section
                    .inherit_deformation_for_uv_counterterm
                    || uv_index == 0
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

                        shifts[k_def_index] = utils::evaluate_signature(
                            &lmm.1,
                            &external_momenta[..self.external_momenta.len()],
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
                        def_jacobian *= jac_def.conj();
                    } else {
                        def_jacobian *= jac_def;
                    }
                } else {
                    k_def_index += subgraph.n_loops;
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

                subgraph_cache.cached_topology_integrand.clear();
            }

            // convert from the cut basis to the loop momentum basis
            if self.numerator_in_loop_momentum_basis || cfg!(feature = "mg_numerator") {
                if let Some(m) = &cut_uv_limit.cb_to_lmb {
                    for (kl, r) in k_def_lmb[..self.n_loops]
                        .iter_mut()
                        .zip_eq(m.chunks(self.n_loops))
                    {
                        *kl = LorentzVector::default();
                        for ((&sign, mom), shift) in r
                            .iter()
                            .zip_eq(k_def.iter())
                            .zip_eq(&shifts[..self.n_loops])
                        {
                            if sign != 0 {
                                *kl += (*mom - shift).multiply_sign(sign);
                            }
                        }
                    }
                } else {
                    panic!("The numerator is in the loop momentum basis but no map is defined.");
                }
            }

            // evaluate the MG numerator
            let mut num_computed = false;
            #[cfg(feature = "mg_numerator")]
            {
                let mut mg_numerator = mem::replace(&mut self.mg_numerator.mg_numerator, None);
                let mut energy_polynomial_coefficients =
                    mem::replace(&mut cut_uv_limit.energy_polynomial_coefficients, vec![]);

                if let Some(call_signature) = &self.mg_numerator.call_signature {
                    let mut all_external_momenta: ArrayVec<[f64; (MAX_LOOP + 4) * 8]> =
                        ArrayVec::new();

                    for m in &self.external_momenta[..self.n_incoming_momenta] {
                        all_external_momenta
                            .try_extend_from_slice(&[m.t, 0., m.x, 0., m.y, 0., m.z, 0.])
                            .unwrap();
                    }

                    for m in &k_def_lmb[..self.n_loops] {
                        all_external_momenta
                            .try_extend_from_slice(&[
                                m.t.re.to_f64().unwrap(),
                                m.t.im.to_f64().unwrap(),
                                m.x.re.to_f64().unwrap(),
                                m.x.im.to_f64().unwrap(),
                                m.y.re.to_f64().unwrap(),
                                m.y.im.to_f64().unwrap(),
                                m.z.re.to_f64().unwrap(),
                                m.z.im.to_f64().unwrap(),
                            ])
                            .unwrap();
                    }

                    let n_incoming_momenta = self.n_incoming_momenta;
                    let n_loops = self.n_loops;
                    let ltd_start_index = cutkosky_cuts.cuts.len() - 1;

                    energy_polynomial_coefficients.clear();
                    for s in &cut_uv_limit.energy_polynomial_samples {
                        // set the energies of the LTD momenta
                        if let Some(m) = &cut_uv_limit.cb_to_lmb {
                            for (i, r) in m.chunks(n_loops).enumerate() {
                                all_external_momenta[(n_incoming_momenta + i) * 8] = 0.;
                                all_external_momenta[(n_incoming_momenta + i) * 8 + 1] = 0.;
                                // x is only for the LTD momenta, so we need to
                                for (j, ((&sign, mom), shift)) in r
                                    .iter()
                                    .zip_eq(k_def.iter())
                                    .zip_eq(&shifts[..n_loops])
                                    .enumerate()
                                {
                                    if sign != 0 {
                                        if j < ltd_start_index {
                                            all_external_momenta[(n_incoming_momenta + i) * 8] +=
                                                (mom.t - shift.t)
                                                    .multiply_sign(sign)
                                                    .re
                                                    .to_f64()
                                                    .unwrap();
                                        } else {
                                            all_external_momenta[(n_incoming_momenta + i) * 8] +=
                                                (s[j - ltd_start_index]
                                                    - shift.t.to_f64().unwrap())
                                                .multiply_sign(sign);
                                        }
                                    }
                                }
                            }
                        }

                        let num = MGNumeratorMod::get_numerator(
                            mg_numerator.as_mut().unwrap(),
                            &all_external_momenta,
                            call_signature.proc_id,
                            call_signature.left_diagram_id,
                            call_signature.right_diagram_id,
                        );

                        energy_polynomial_coefficients.push(num);
                    }

                    // multiply with the matrix to get the coefficients of the energy polynomial
                    // and put them in the correct place in the reduced polynomial
                    for c in &mut diag_cache[0].reduced_coefficient_lb[0] {
                        *c = Complex::zero();
                    }

                    for (pos, r) in cut_uv_limit
                        .energy_polynomial_map
                        .iter_mut()
                        .zip_eq(&cut_uv_limit.energy_polynomial_matrix)
                    {
                        let mut c: Complex<f64> = Complex::zero();
                        for (v, eval) in r.iter().zip_eq(&energy_polynomial_coefficients) {
                            c += eval * v;
                        }

                        if *pos >= diag_cache[0].reduced_coefficient_lb[0].len() {
                            diag_cache[0].reduced_coefficient_lb[0]
                                .resize(*pos + 1, Complex::zero());
                        }

                        diag_cache[0].reduced_coefficient_lb[0][*pos] =
                            Complex::new(Into::<T>::into(c.re), Into::<T>::into(c.im));
                    }

                    num_computed = true;
                }

                if self.form_numerator.form_numerator.is_some() {
                    let mut form_numerator =
                        mem::replace(&mut self.form_numerator.form_numerator, None);
                    let mut form_numerator_buffer =
                        mem::replace(&mut self.form_numerator.form_numerator_buffer, vec![]);
                    if let Some(call_signature) = &self.form_numerator.call_signature {
                        let mut scalar_products: ArrayVec<[f64; MAX_LOOP * MAX_LOOP * 6]> =
                            ArrayVec::new();

                        for (i1, e1) in self.external_momenta[..self.n_incoming_momenta]
                            .iter()
                            .enumerate()
                        {
                            scalar_products.push(e1.t);
                            scalar_products.push(0.);
                            for e2 in &self.external_momenta[i1..self.n_incoming_momenta] {
                                scalar_products.push(e1.dot(e2));
                                scalar_products.push(0.);
                            }
                        }

                        for (i1, m1) in k_def_lmb[..self.n_loops].iter().enumerate() {
                            scalar_products.push(m1.t.re.to_f64().unwrap());
                            scalar_products.push(m1.t.im.to_f64().unwrap());
                            for e1 in &self.external_momenta[..self.n_incoming_momenta] {
                                let d = m1.dot(&e1.cast());
                                scalar_products.push(d.re.to_f64().unwrap());
                                scalar_products.push(d.im.to_f64().unwrap());
                                let d = m1.spatial_dot(&e1.cast());
                                scalar_products.push(d.re.to_f64().unwrap());
                                scalar_products.push(d.im.to_f64().unwrap());
                            }

                            for m2 in k_def_lmb[i1..self.n_loops].iter() {
                                let d = m1.dot(m2);
                                scalar_products.push(d.re.to_f64().unwrap());
                                scalar_products.push(d.im.to_f64().unwrap());
                                let d = m1.spatial_dot(m2);
                                scalar_products.push(d.re.to_f64().unwrap());
                                scalar_products.push(d.im.to_f64().unwrap());
                            }
                        }

                        if form_numerator_buffer.len()  < 2 {
                            form_numerator_buffer.resize(2, 0.);
                        }

                        FORMNumeratorMod::get_numerator(
                            form_numerator.as_mut().unwrap(),
                            &scalar_products,
                            call_signature.id,
                            0, // only evaluate the constant part for now
                            &mut form_numerator_buffer,
                        );
                        let result = Complex::new(Into::<T>::into(form_numerator_buffer[0]), Into::<T>::into(form_numerator_buffer[1]));

                        // compare against the MG numerator
                        if self.mg_numerator.call_signature.is_some() {
                            if (diag_cache[0].reduced_coefficient_lb[0][0] - result).norm_sqr()
                                > Into::<T>::into(1e-10 * self.e_cm_squared)
                            {
                                println!(
                                    "Mismatch between MG and RUST numerator: {} vs {} for {:?}",
                                    diag_cache[0].reduced_coefficient_lb[0][0],
                                    result,
                                    loop_momenta
                                );
                            }
                        }

                        diag_cache[0].reduced_coefficient_lb[0].clear();
                        diag_cache[0].reduced_coefficient_lb[0].push(Complex::new(
                            Into::<T>::into(result.re),
                            Into::<T>::into(result.im),
                        ));
                        num_computed = true;
                    }

                    mem::swap(&mut form_numerator, &mut self.form_numerator.form_numerator);
                    mem::swap(&mut form_numerator_buffer, &mut self.form_numerator.form_numerator_buffer);
                }

                mem::swap(&mut mg_numerator, &mut self.mg_numerator.mg_numerator);
                mem::swap(
                    &mut energy_polynomial_coefficients,
                    &mut cut_uv_limit.energy_polynomial_coefficients,
                );
            }

            if !num_computed {
                // evaluate the cut numerator with the spatial parts of all loop momenta and
                // the energy parts of the Cutkosky cuts.
                // Store the reduced numerator in the left graph cache for now
                // TODO: this is impractical
                cut_uv_limit.numerator.evaluate_reduced_in_lb(
                    if self.numerator_in_loop_momentum_basis {
                        &k_def_lmb[..self.n_loops]
                    } else {
                        &k_def
                    },
                    cutkosky_cuts.cuts.len() - 1, // FIXME: this is not correct in the loop momentum basis!
                    &mut diag_cache[0],
                    0,
                    regenerate_momenta,
                    self.settings.general.partial_fractioning_threshold > 0.0,
                );
                regenerate_momenta = false;
            }

            let left_cache = &mut diag_cache[0];
            mem::swap(
                &mut left_cache.reduced_coefficient_lb_supergraph[0],
                &mut left_cache.reduced_coefficient_lb[0],
            );

            let mut supergraph_coeff = mem::replace(
                &mut diag_cache[0].reduced_coefficient_lb_supergraph[0],
                vec![],
            );

            // evaluate the subgraphs for every monomial in the numerator
            let mut result_complete_numerator = Complex::default();
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

                    for n in &mut subgraph.numerator.coefficients {
                        *n = Complex::zero();
                    }

                    let mut monomial_index = 0;

                    if subgraph.n_loops == 0 {
                        subgraph.numerator.coefficients[0] = Complex::one();
                        subgraph.numerator.max_rank = 0;
                    } else {
                        // now set the coefficient to 1 for the current monomial in the subgraph
                        // by mapping the powers in the reduced numerator back to the power pattern
                        // of the complete numerator
                        // there is partial support for numerators in the loop momentum basis: only permutations of the
                        // order from lmb to cb are allowed
                        let mut subgraph_powers = [0; 10]; // TODO: make max rank a constant
                        let mut power_index = 0;

                        if self.numerator_in_loop_momentum_basis {
                            for (lmi, (lmm, ext)) in subgraph.loop_momentum_map.iter().enumerate() {
                                debug_assert!(ext.iter().all(|c| *c == 0));
                                debug_assert!(lmm.iter().filter(|c| **c != 0).count() == 1);

                                let lmb_index = lmm.iter().filter(|c| **c != 0).next().unwrap();
                                for _ in 0..powers[*lmb_index as usize] {
                                    subgraph_powers[power_index] = lmi * 4;
                                    power_index += 1;
                                }
                            }
                        } else {
                            for (lmi, p) in powers[def_mom_index..def_mom_index + subgraph.n_loops]
                                .iter()
                                .enumerate()
                            {
                                for _ in 0..*p {
                                    subgraph_powers[power_index] = lmi * 4;
                                    power_index += 1;
                                }
                            }
                        }

                        subgraph.numerator.max_rank = power_index; // FIXME: is this right?

                        if power_index == 0 {
                            subgraph.numerator.coefficients[0] = Complex::one();
                        } else {
                            monomial_index = subgraph.numerator.sorted_linear[..]
                                .binary_search_by(|x| {
                                    utils::compare_slice(&x[..], &subgraph_powers[..power_index])
                                })
                                .unwrap()
                                + 1;
                            subgraph.numerator.coefficients[monomial_index] += Complex::one();
                        }
                    }

                    let reduced_pos = LTDNumerator::powers_to_position(
                        &subgraph.numerator.coefficient_index_to_powers[monomial_index]
                            [..subgraph.numerator.n_loops],
                        subgraph.numerator.n_loops,
                        &subgraph.numerator.reduced_blocks,
                    );

                    // update the non-empty coeff map to the single entry
                    subgraph.numerator.non_empty_coeff_map_to_reduced_numerator[0].clear();
                    subgraph.numerator.non_empty_coeff_map_to_reduced_numerator[0].push((
                        monomial_index,
                        reduced_pos,
                        Complex::one(),
                    ));

                    // only compute the subdiagram with this numerator once
                    let mut cached_res = None;
                    for (i, r) in &subgraph_cache.cached_topology_integrand {
                        if *i == monomial_index {
                            cached_res = Some(*r);
                            break;
                        }
                    }

                    let mut res = if cached_res.is_none() {
                        let res = if subgraph.loop_lines.len() > 0 {
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

                        subgraph_cache
                            .cached_topology_integrand
                            .push((monomial_index, res));
                        res
                    } else {
                        cached_res.unwrap()
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
                            "  | monomial res {} ({}l) = {:e}",
                            subgraph.name, subgraph.n_loops, res
                        );
                    }

                    def_mom_index += subgraph.n_loops;
                }
                result_complete_numerator += num_result * def_jacobian;
            }

            mem::swap(
                &mut supergraph_coeff,
                &mut diag_cache[0].reduced_coefficient_lb_supergraph[0],
            );

            if self.settings.general.debug >= 1 {
                println!(
                    "  | res uv limit {}: = {:e}",
                    uv_index, result_complete_numerator
                );
            }

            diag_and_num_contributions += result_complete_numerator;
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
}

impl IntegrandImplementation for SquaredTopologySet {
    type Cache = SquaredTopologyCache;

    /// Create a rotated version of this squared topology. The axis needs to be normalized.
    fn rotate(&self, angle: float, axis: (float, float, float)) -> SquaredTopologySet {
        let rotated_topologies: Vec<_> = self
            .topologies
            .iter()
            .map(|t| t.rotate(angle, axis))
            .collect();

        let mut c = self.multi_channeling_channels.clone();

        let rot_matrix = &self.rotation_matrix;
        for (_, _, shifts) in &mut c {
            for shift in shifts.iter_mut() {
                shift.x = (rot_matrix[0][0] * shift.x
                    + rot_matrix[0][1] * shift.y
                    + rot_matrix[0][2] * shift.z)
                    .to_f64()
                    .unwrap();
                shift.y = (rot_matrix[1][0] * shift.x
                    + rot_matrix[1][1] * shift.y
                    + rot_matrix[1][2] * shift.z)
                    .to_f64()
                    .unwrap();
                shift.z = (rot_matrix[2][0] * shift.x
                    + rot_matrix[2][1] * shift.y
                    + rot_matrix[2][2] * shift.z)
                    .to_f64()
                    .unwrap();
            }
        }

        SquaredTopologySet {
            name: self.name.clone(),
            multiplicity: self.multiplicity.clone(),
            e_cm_squared: self.e_cm_squared,
            n_loops: self.n_loops.clone(),
            settings: self.settings.clone(),
            rotation_matrix: rotated_topologies[0].rotation_matrix.clone(),
            topologies: rotated_topologies,
            multi_channeling_channels: c,
        }
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
            float_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
            quad_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
        }
    }
}
