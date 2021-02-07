use crate::{dashboard::{StatusUpdate, StatusUpdateSender}, topologies};
use crate::integrand::{IntegrandImplementation, IntegrandSample};
use crate::observables::EventManager;
use crate::topologies::FixedDeformationLimit;
use crate::topologies::{Cut, LTDCache, LTDNumerator, Topology};
use crate::utils;
use crate::{
    float, DeformationStrategy, FloatLike, IRHandling, IntegrandType, NormalisingFunction,
    NumeratorSource, Settings,
};
use arrayvec::ArrayVec;
use color_eyre::{Help, Report};
use dlopen::wrapper::Container;
use eyre::WrapErr;
use f128::f128;
use havana::{ContinuousGrid, DiscreteGrid, Grid};
use itertools::Itertools;
use lorentz_vector::{LorentzVector, RealNumberLike};
use num::Complex;
use num_traits::{Float, FloatConst, FromPrimitive, Inv, NumCast, One, ToPrimitive, Zero};
use rand::{thread_rng, Rng};
use serde::{Deserialize, __private::de};
use std::collections::HashMap;
use std::fs::File;
use std::mem;
use std::path::Path;
use std::time::Instant;
use utils::Signum;

#[cfg(not(feature = "higher_loops"))]
pub const MAX_SG_LOOP: usize = 4;
#[cfg(feature = "higher_loops")]
pub const MAX_SG_LOOP: usize = 10;

mod form_numerator {
    use dlopen::wrapper::{Container, WrapperApi};
    use libc::{c_double, c_int};
    use std::path::PathBuf;

    pub trait GetNumerator {
        fn get_numerator(
            api_container: &mut Container<FORMNumeratorAPI>,
            p: &[Self],
            params: &[Self],
            diag: usize,
            conf: usize,
            poly: &mut [Self],
        ) -> usize
        where
            Self: std::marker::Sized;
    }

    impl GetNumerator for f64 {
        fn get_numerator(
            api_container: &mut Container<FORMNumeratorAPI>,
            p: &[f64],
            params: &[f64],
            diag: usize,
            conf: usize,
            poly: &mut [f64],
        ) -> usize {
            unsafe {
                api_container.evaluate(
                    &p[0] as *const f64,
                    &params[0] as *const f64,
                    diag as i32,
                    conf as i32,
                    &mut poly[0] as *mut f64,
                ) as usize
            }
        }
    }

    impl GetNumerator for f128::f128 {
        fn get_numerator(
            api_container: &mut Container<FORMNumeratorAPI>,
            p: &[f128::f128],
            params: &[f128::f128],
            diag: usize,
            conf: usize,
            poly: &mut [f128::f128],
        ) -> usize {
            unsafe {
                api_container.evaluate_f128(
                    &p[0] as *const f128::f128,
                    &params[0] as *const f128::f128,
                    diag as i32,
                    conf as i32,
                    &mut poly[0] as *mut f128::f128,
                ) as usize
            }
        }
    }

    #[derive(WrapperApi)]
    pub struct FORMNumeratorAPI {
        evaluate: unsafe extern "C" fn(
            p: *const c_double,
            params: *const c_double,
            diag: c_int,
            conf: c_int,
            out: *mut c_double,
        ) -> c_int,
        evaluate_f128: unsafe extern "C" fn(
            p: *const f128::f128,
            params: *const f128::f128,
            diag: c_int,
            conf: c_int,
            out: *mut f128::f128,
        ) -> c_int,
        get_buffer_size: unsafe extern "C" fn() -> c_int,
        get_rank: unsafe extern "C" fn(diag: c_int, conf: c_int) -> c_int,
    }

    pub fn get_buffer_size(api_container: &mut Container<FORMNumeratorAPI>) -> usize {
        unsafe { api_container.get_buffer_size() as usize }
    }

    pub fn get_rank(
        api_container: &mut Container<FORMNumeratorAPI>,
        diag: usize,
        conf: usize,
    ) -> usize {
        unsafe { api_container.get_rank(diag as i32, conf as i32) as usize }
    }

    pub fn load(base_path: &str) -> Container<FORMNumeratorAPI> {
        let mut lib_path = PathBuf::from(&base_path);
        lib_path.push("lib/libFORM_numerators.so");
        let container: Container<FORMNumeratorAPI> =
            unsafe { Container::load(lib_path) }.expect("Could not open library or load symbols");

        container
    }
}

mod form_integrand {
    use dlopen::wrapper::{Container, WrapperApi};
    use libc::{c_double, c_int};
    use num::Complex;
    use std::path::PathBuf;

    pub trait GetIntegrand {
        fn get_integrand_ltd(
            api_container: &mut Container<FORMIntegrandAPI>,
            p: &[Self],
            params: &[Self],
            diag: usize,
            conf: usize,
        ) -> Complex<Self>
        where
            Self: std::marker::Sized;

        fn get_integrand_pf(
            api_container: &mut Container<FORMIntegrandAPI>,
            p: &[Self],
            params: &[Self],
            diag: usize,
            conf: usize,
        ) -> Complex<Self>
        where
            Self: std::marker::Sized;
    }

    impl GetIntegrand for f64 {
        fn get_integrand_ltd(
            api_container: &mut Container<FORMIntegrandAPI>,
            p: &[f64],
            params: &[f64],
            diag: usize,
            conf: usize,
        ) -> Complex<f64> {
            unsafe {
                let mut c = Complex::default();
                api_container.evaluate_ltd(
                    &p[0] as *const f64,
                    &params[0] as *const f64,
                    diag as i32,
                    conf as i32,
                    &mut c as *mut Complex<f64>,
                );
                c
            }
        }
        fn get_integrand_pf(
            api_container: &mut Container<FORMIntegrandAPI>,
            p: &[f64],
            params: &[f64],
            diag: usize,
            conf: usize,
        ) -> Complex<f64> {
            unsafe {
                let mut c = Complex::default();
                api_container.evaluate_pf(
                    &p[0] as *const f64,
                    &params[0] as *const f64,
                    diag as i32,
                    conf as i32,
                    &mut c as *mut Complex<f64>,
                );
                c
            }
        }
    }

    impl GetIntegrand for f128::f128 {
        fn get_integrand_ltd(
            api_container: &mut Container<FORMIntegrandAPI>,
            p: &[f128::f128],
            params: &[f128::f128],
            diag: usize,
            conf: usize,
        ) -> Complex<f128::f128> {
            unsafe {
                let mut c = Box::new(Complex::default());
                api_container.evaluate_ltd_f128(
                    &p[0] as *const f128::f128,
                    &params[0] as *const f128::f128,
                    diag as i32,
                    conf as i32,
                    c.as_mut() as *mut Complex<f128::f128>,
                );
                *c
            }
        }

        fn get_integrand_pf(
            api_container: &mut Container<FORMIntegrandAPI>,
            p: &[f128::f128],
            params: &[f128::f128],
            diag: usize,
            conf: usize,
        ) -> Complex<f128::f128> {
            unsafe {
                let mut c = Box::new(Complex::default());
                api_container.evaluate_pf_f128(
                    &p[0] as *const f128::f128,
                    &params[0] as *const f128::f128,
                    diag as i32,
                    conf as i32,
                    c.as_mut() as *mut Complex<f128::f128>,
                );
                *c
            }
        }
    }

    #[derive(WrapperApi)]
    pub struct FORMIntegrandAPI {
        #[dlopen_name = "evaluate_LTD"]
        evaluate_ltd: unsafe extern "C" fn(
            p: *const c_double,
            params: *const c_double,
            diag: c_int,
            conf: c_int,
            out: *mut Complex<c_double>,
        ),
        #[dlopen_name = "evaluate_LTD_f128"]
        evaluate_ltd_f128: unsafe extern "C" fn(
            p: *const f128::f128,
            params: *const f128::f128,
            diag: c_int,
            conf: c_int,
            out: *mut Complex<f128::f128>,
        ),
        #[dlopen_name = "evaluate_PF"]
        evaluate_pf: unsafe extern "C" fn(
            p: *const c_double,
            params: *const c_double,
            diag: c_int,
            conf: c_int,
            out: *mut Complex<c_double>,
        ),
        #[dlopen_name = "evaluate_PF_f128"]
        evaluate_pf_f128: unsafe extern "C" fn(
            p: *const f128::f128,
            params: *const f128::f128,
            diag: c_int,
            conf: c_int,
            out: *mut Complex<f128::f128>,
        ),
    }

    pub fn load(base_path: &str) -> Container<FORMIntegrandAPI> {
        let mut lib_path = PathBuf::from(&base_path);
        lib_path.push("lib/libFORM_integrands.so");
        let container: Container<FORMIntegrandAPI> =
            unsafe { Container::load(lib_path) }.expect("Could not open library or load symbols");

        container
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCut {
    pub name: String,
    pub sign: i8,
    pub power: usize,
    pub particle_id: isize,
    pub signature: (Vec<i8>, Vec<i8>),
    pub m_squared: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCutDiagramInfo {
    pub graph: Topology,
    pub conjugate_deformation: bool,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCutDiagramSet {
    pub id: usize,
    pub diagram_info: Vec<CutkoskyCutDiagramInfo>,
    #[serde(default)]
    pub numerator_tensor_coefficients_sparse: Vec<(Vec<usize>, (f64, f64))>,
    #[serde(default)]
    pub numerator_tensor_coefficients: Vec<(f64, f64)>,
    #[serde(skip_deserializing)]
    pub numerator: LTDNumerator,
    pub cb_to_lmb: Option<Vec<i8>>,
    #[serde(skip_deserializing)]
    pub energy_polynomial_matrix: Vec<Vec<f64>>,
    #[serde(skip_deserializing)]
    pub energy_polynomial_map: Vec<usize>,
    #[serde(skip_deserializing)]
    pub energy_polynomial_samples: Vec<Vec<f64>>,
    #[serde(skip_deserializing)]
    pub energy_polynomial_coefficients: Vec<Complex<f64>>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct CutkoskyCuts {
    pub cuts: Vec<CutkoskyCut>,
    pub diagram_sets: Vec<CutkoskyCutDiagramSet>,
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

#[derive(Debug, Clone, Deserialize)]
pub struct FORMIntegrandCallSignature {
    pub id: usize,
}

#[derive(Deserialize)]
pub struct FORMNumerator {
    call_signature: Option<FORMNumeratorCallSignature>,
    #[serde(skip_deserializing)]
    pub form_numerator: Option<Container<form_numerator::FORMNumeratorAPI>>,
    #[serde(skip_deserializing)]
    pub form_numerator_buffer: Vec<f64>,
    #[serde(skip_deserializing)]
    pub form_numerator_buffer_size: usize,
    #[serde(default)]
    pub base_path: String,
}

impl Clone for FORMNumerator {
    fn clone(&self) -> Self {
        FORMNumerator {
            call_signature: self.call_signature.clone(),
            base_path: self.base_path.clone(),
            form_numerator: if self.form_numerator.is_some() {
                Some(form_numerator::load(&self.base_path))
            } else {
                None
            },
            form_numerator_buffer: self.form_numerator_buffer.clone(),
            form_numerator_buffer_size: self.form_numerator_buffer_size,
        }
    }
}

#[derive(Deserialize)]
pub struct FORMIntegrand {
    call_signature: Option<FORMIntegrandCallSignature>,
    #[serde(skip_deserializing)]
    pub form_integrand: Option<Container<form_integrand::FORMIntegrandAPI>>,
    #[serde(default)]
    pub base_path: String,
}

impl Clone for FORMIntegrand {
    fn clone(&self) -> Self {
        FORMIntegrand {
            call_signature: self.call_signature.clone(),
            base_path: self.base_path.clone(),
            form_integrand: if self.form_integrand.is_some() {
                Some(form_integrand::load(&self.base_path))
            } else {
                None
            },
        }
    }
}

#[derive(Clone, Deserialize)]
pub struct MultiChannelingBasis {
    pub defining_propagators: Vec<(usize, usize)>,
}
// MODIFY FOR AMPLITUDES -> copy dict in form of structs: remember clone and desirialize #derive()  : One SG
#[derive(Clone, Deserialize)]
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
    pub topo: Topology,
    pub default_fixed_cut_momenta: (Vec<LorentzVector<f64>>, Vec<LorentzVector<f64>>), // If amplitude: overwrite with info from external data: (incoming, outgoing), nfinal-1 indpendent fixed
    #[serde(default)]
    pub multi_channeling_bases: Vec<MultiChannelingBasis>,
    #[serde(skip_deserializing)]
    pub multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f64>>)>,
    #[serde(rename = "FORM_numerator")]
    pub form_numerator: FORMNumerator,
    #[serde(rename = "FORM_integrand")]
    pub form_integrand: FORMIntegrand,
    #[serde(skip_deserializing)]
    pub is_stability_check_topo: bool,
    #[serde(skip_deserializing)]
    pub pol: Option<Vec<LorentzVector<Complex<f64>>>>,
    #[serde(skip_deserializing)]
    pub cpol: Option<Vec<LorentzVector<Complex<f64>>>>,
    #[serde(skip_deserializing)]
    pub spinor_u: Option<Vec<Vec<Complex<f64>>>>,
    #[serde(skip_deserializing)]
    pub spinor_ubar: Option<Vec<Vec<Complex<f64>>>>,
    #[serde(skip_deserializing)]
    pub spinor_v: Option<Vec<Vec<Complex<f64>>>>,
    #[serde(skip_deserializing)]
    pub spinor_vbar: Option<Vec<Vec<Complex<f64>>>>,
    pub is_amplitude: Option<bool>,
}

// accomodate <- external data from here: read in set -> feed to squared topology (set up struc in squaredtopology as well)
#[derive(Clone)]
pub struct SquaredTopologySet {
    pub name: String,
    pub e_cm_squared: f64,
    topologies: Vec<SquaredTopology>,
    additional_topologies: Vec<Vec<(SquaredTopology, Vec<(Vec<i8>, Vec<i8>)>)>>,
    pub multiplicity: Vec<f64>,
    pub settings: Settings,
    pub rotation_matrix: [[float; 3]; 3],
    pub multi_channeling_channels: Vec<(Vec<i8>, Vec<i8>, Vec<LorentzVector<f64>>)>,
    pub stability_check_topologies: Vec<Vec<SquaredTopology>>,
    pub is_stability_check_topo: bool,
}


#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetAdditionalTopology {
    pub name: String,
    pub this_lmb_to_defining_lmb: Vec<(Vec<i8>, Vec<i8>)>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetTopology {
    pub name: String,
    pub multiplicity: f64,
    #[serde(default, rename = "additional_LMBs")]
    pub additional_lmbs: Vec<SquaredTopologySetAdditionalTopology>,
    #[serde(default)]
    pub stability_check_topologies: Vec<String>,
}


// FOR amplitudes
#[derive(Debug, Deserialize,Clone)]
pub struct ExternalData {
    pub in_momenta: Vec<LorentzVector<f64>>,
    pub out_momenta: Vec<LorentzVector<f64>>,
    pub pol: Vec<Vec<Vec<f64>>> ,
    pub cpol: Vec<Vec<Vec<f64>>>,
    pub spinor_u: Vec<Vec<Vec<f64>>> ,
    pub spinor_ubar: Vec<Vec<Vec<f64>>> ,
    pub spinor_v: Vec<Vec<Vec<f64>>>  ,
    pub spinor_vbar: Vec<Vec<Vec<f64>>>  ,
    pub  n_in : u32,
    pub n_out : u32,
}
impl Default for ExternalData {
    fn default() -> ExternalData {
        ExternalData {
            in_momenta: vec![],
            out_momenta: vec![],
            pol: vec![],
            cpol: vec![],
            spinor_u: vec![],
            spinor_ubar: vec![],
            spinor_v: vec![],
            spinor_vbar: vec![],
            n_in: 0,
            n_out:0
        }
    }
}
// HELPER FUNCTIONS FOR INPUT CONVERSION
pub fn from_vec_vec_to_cmplx_lv(input:Vec<Vec<f64>>) -> LorentzVector<Complex<f64>>  {
    let v = input.clone();
    let (t_re, x_re, y_re, z_re) = (v[0][0], v[1][0], v[2][0], v[3][0]);
    let (t_im, x_im, y_im, z_im) = (v[0][1], v[1][1], v[2][1], v[3][1]);
    
    let l_vec = LorentzVector {
            t: Complex::new(t_re, t_im),
            x: Complex::new(x_re, x_im),
            y: Complex::new(y_re, y_im),
            z: Complex::new(z_re, z_im),
        };
    l_vec
    
    }
pub fn from_vec_vec_to_cmplx_v(input:Vec<Vec<f64>>) -> Vec<Complex<f64>>  {
        let v = input.clone();
        let (t_re, x_re, y_re, z_re) = (v[0][0], v[1][0], v[2][0], v[3][0]);
        let (t_im, x_im, y_im, z_im) = (v[0][1], v[1][1], v[2][1], v[3][1]);
        
        let vec = vec![ 
                Complex::new(t_re, t_im),
                Complex::new(x_re, x_im),
                Complex::new(y_re, y_im),
                Complex::new(z_re, z_im)];
        vec
        
        }
// HELPER FUNCTIONS FOR INPUT CONVERSION: END

#[derive(Debug, Clone, Deserialize)]
pub struct SquaredTopologySetInput {
    pub name: String,
    topologies: Vec<SquaredTopologySetTopology>,
    pub external_data: Option<ExternalData>,
}
impl Default for SquaredTopologySetInput {
    fn default() -> SquaredTopologySetInput {
        SquaredTopologySetInput {
            name : Default::default(),
            external_data: Default::default(),
            topologies: vec![],
        }
    }
}


impl SquaredTopologySet {
    pub fn from_one(mut squared_topology: SquaredTopology) -> SquaredTopologySet {
        squared_topology.generate_multi_channeling_channels();

        SquaredTopologySet {
            name: squared_topology.name.clone(),
            settings: squared_topology.settings.clone(),
            e_cm_squared: squared_topology.e_cm_squared,
            rotation_matrix: squared_topology.rotation_matrix.clone(),
            multi_channeling_channels: squared_topology.multi_channeling_channels.clone(),
            topologies: vec![squared_topology],
            additional_topologies: vec![vec![]],
            multiplicity: vec![1.],
            stability_check_topologies: vec![vec![]],
            is_stability_check_topo: false,
        }
    }

    pub fn from_file(filename: &str, settings: &Settings) -> Result<SquaredTopologySet, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open squared topology set file {}", filename))
            .suggestion("Does the path exist?")?;
        // Q: Is that enough for having ExternalStructure set, even if it does not exist?
    
        
        let squared_topology_set_input:SquaredTopologySetInput = serde_yaml::from_reader(f)
            .wrap_err("Could not parse squared topology set file")
            .suggestion("Is it a correct yaml file")?;


        let mut topologies: Vec<SquaredTopology> = vec![];
        let mut additional_topologies: Vec<Vec<(SquaredTopology, Vec<(Vec<i8>, Vec<i8>)>)>> =
            vec![];
        let mut multiplicity: Vec<f64> = vec![];
        let mut stability_topologies = vec![];
        let amp_input = squared_topology_set_input.clone();
        

        for topo in squared_topology_set_input.topologies {
            let filename = std::path::Path::new(&filename)
                .with_file_name(topo.name)
                .with_extension("yaml");
            
            let mut squared_topology = SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                .wrap_err("Could not load subtopology file")?;
            // we overwrite here for amplitudes
            squared_topology.set_external_data(& amp_input);
            squared_topology.overwrite_momenta(& amp_input);
            


            let mut additional_topologies_for_topo = vec![];
            for t in topo.additional_lmbs {
                let filename = std::path::Path::new(&filename)
                    .with_file_name(t.name)
                    .with_extension("yaml");
                let additional_squared_topology =
                    SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                        .wrap_err("Could not load subtopology file")?;
                additional_topologies_for_topo
                    .push((additional_squared_topology, t.this_lmb_to_defining_lmb));
            }
            additional_topologies.push(additional_topologies_for_topo);

            let mut stability_topologies_for_topo = vec![];
            for t in topo.stability_check_topologies {
                let filename = std::path::Path::new(&filename)
                    .with_file_name(t)
                    .with_extension("yaml");
                let mut additional_squared_topology =
                    SquaredTopology::from_file(filename.to_str().unwrap(), settings)
                        .wrap_err("Could not load subtopology file")?;
                additional_squared_topology.is_stability_check_topo = true;
                stability_topologies_for_topo.push(additional_squared_topology);
            }
            stability_topologies.push(stability_topologies_for_topo);

            topologies.push(squared_topology);
            multiplicity.push(topo.multiplicity);
        }

        let rotation_matrix = [
            [float::one(), float::zero(), float::zero()],
            [float::zero(), float::one(), float::zero()],
            [float::zero(), float::zero(), float::one()],
        ];

        let mut sts = SquaredTopologySet {
            name: squared_topology_set_input.name,
            e_cm_squared: topologies[0].e_cm_squared,
            topologies,
            additional_topologies,
            rotation_matrix,
            settings: settings.clone(),
            multiplicity,
            multi_channeling_channels: vec![],
            stability_check_topologies: stability_topologies,
            is_stability_check_topo: false,
        };
        sts.create_multi_channeling_channels();
        Ok(sts)
    }

    fn create_multi_channeling_channels(&mut self) {
        let mut multi_channeling_channels = vec![];
        let mut max_loops = 0;
        for t in &mut self.topologies {
            t.generate_multi_channeling_channels();

            if t.n_loops > max_loops {
                multi_channeling_channels.clear();
                max_loops = t.n_loops;
            }
            if t.n_loops == max_loops {
                multi_channeling_channels.extend(t.multi_channeling_channels.iter().cloned());
            }
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

        self.multi_channeling_channels = unique_multi_channeling_channels;
    }

    /// Create a rotated version of this squared topology. The axis needs to be normalized.
    fn rotate(&self, angle: float, axis: (float, float, float)) -> SquaredTopologySet {
        let rotated_topologies: Vec<_> = self
            .topologies
            .iter()
            .map(|t| {
                let mut stab_t = t.rotate(angle, axis);
                stab_t.is_stability_check_topo = true;
                stab_t
            })
            .collect();
        let mut c = self.multi_channeling_channels.clone();

        let rot_matrix = &rotated_topologies[0].rotation_matrix;
        for (_, _, shifts) in &mut c {
            for shift in shifts.iter_mut() {
                let old_shift = shift.clone();
                shift.x = (rot_matrix[0][0] * old_shift.x
                    + rot_matrix[0][1] * old_shift.y
                    + rot_matrix[0][2] * old_shift.z)
                    .to_f64()
                    .unwrap();
                shift.y = (rot_matrix[1][0] * old_shift.x
                    + rot_matrix[1][1] * old_shift.y
                    + rot_matrix[1][2] * old_shift.z)
                    .to_f64()
                    .unwrap();
                shift.z = (rot_matrix[2][0] * old_shift.x
                    + rot_matrix[2][1] * old_shift.y
                    + rot_matrix[2][2] * old_shift.z)
                    .to_f64()
                    .unwrap();
            }
        }

        SquaredTopologySet {
            name: self.name.clone(),
            multiplicity: self.multiplicity.clone(),
            e_cm_squared: self.e_cm_squared,
            settings: self.settings.clone(),
            rotation_matrix: rotated_topologies[0].rotation_matrix.clone(),
            additional_topologies: vec![vec![]; rotated_topologies.len()], // rotated versions of additional topologies are not supported
            topologies: rotated_topologies,
            multi_channeling_channels: c,
            stability_check_topologies: vec![vec![]; self.topologies.len()],
            is_stability_check_topo: true,
            // // For amplitudes
            // cpol: self.cpol.clone(),
            // pol: self.pol.clone(),
            // spinor_u: self.spinor_u.clone(),
            // spinor_ubar: self.spinor_ubar.clone(),
            // spinor_v: self.spinor_v.clone(),
            // spinor_vbar: self.spinor_vbar.clone(),
            // in_momenta:self.in_momenta.clone(),
            // out_momenta:self.out_momenta.clone(),
            

        }
    }

    pub fn get_maximum_loop_count(&self) -> usize {
        self.topologies
            .iter()
            .map(|t| t.n_loops - self.settings.cross_section.fixed_cut_momenta.len())
            .max()
            .unwrap()
    }

    pub fn create_caches<T: FloatLike>(&self) -> SquaredTopologyCache<T> {
        SquaredTopologyCache {
            topology_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
            deformation_vector_cache: vec![],
            scalar_products: vec![],
            current_supergraph: 0,
            current_deformation_index: 0,
        }
    }

    pub fn print_info(&self, status_update_sender: &mut StatusUpdateSender) {
        /*
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Overall number of loops in all supergraphs: {}",
                self.n_loops
            )))
            .unwrap();
        */
        let n_topologies = self.topologies.len();

        let mut n_cutkosky_cuts_per_cut_cardinality = HashMap::new();
        let mut n_diagrams_per_loop = HashMap::new();
        let mut n_ltd_cuts_per_loop = HashMap::new();

        for topology in &self.topologies {
            for cutkosky_cuts in &topology.cutkosky_cuts {
                *n_cutkosky_cuts_per_cut_cardinality
                    .entry(cutkosky_cuts.cuts.len())
                    .or_insert(0) += 1;
                for diagram_set in &cutkosky_cuts.diagram_sets {
                    for d_info in &diagram_set.diagram_info {
                        *n_diagrams_per_loop.entry(d_info.graph.n_loops).or_insert(0) += 1;
                        *n_ltd_cuts_per_loop.entry(d_info.graph.n_loops).or_insert(0) += d_info
                            .graph
                            .ltd_cut_options
                            .iter()
                            .map(|v| v.len())
                            .sum::<usize>()
                            .max(1);
                    }
                }
            }
        }

        let mut tmp_vec: Vec<_> = n_ltd_cuts_per_loop.into_iter().collect();
        let tmp_sum: usize = tmp_vec.iter().map(|(_, v)| *v).sum();
        tmp_vec.sort_by(|a, b| a.0.cmp(&b.0));
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Loop count -> number of LTD cuts: {} total: {}",
                tmp_vec
                    .iter()
                    .map(|(k, v)| format!("{}->{}", k, v))
                    .collect::<Vec<String>>()
                    .join(", "),
                tmp_sum
            )))
            .unwrap();

        let mut tmp_vec: Vec<_> = n_diagrams_per_loop.into_iter().collect();
        let tmp_sum: usize = tmp_vec.iter().map(|(_, v)| *v).sum();
        tmp_vec.sort_by(|a, b| a.0.cmp(&b.0));
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Loop count -> number of subdiagrams: {} total: {}",
                tmp_vec
                    .iter()
                    .map(|(k, v)| format!("{}->{}", k, v))
                    .collect::<Vec<String>>()
                    .join(", "),
                tmp_sum
            )))
            .unwrap();

        let mut tmp_vec: Vec<(usize, usize)> =
            n_cutkosky_cuts_per_cut_cardinality.into_iter().collect();
        let tmp_sum: usize = tmp_vec.iter().map(|(_, v)| *v).sum();
        tmp_vec.sort_by(|a, b| a.0.cmp(&b.0));
        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Cut cardinality -> number of Cutkosky cuts: {} total: {}",
                tmp_vec
                    .iter()
                    .map(|(k, v)| format!("{}->{}", k, v))
                    .collect::<Vec<String>>()
                    .join(", "),
                tmp_sum
            )))
            .unwrap();

        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Number of integration channels: {}",
                self.multi_channeling_channels.len()
            )))
            .unwrap();

        status_update_sender
            .send(StatusUpdate::Message(format!(
                "Number of supergraphs: {}",
                n_topologies
            )))
            .unwrap();

        if n_topologies == 1 {
            if let Some(imag_res) = self.topologies[0].topo.analytical_result_imag {
                if imag_res != 0.0 {
                    status_update_sender
                        .send(StatusUpdate::Message(format!(
                            "Target benchmark result (imag) : {:+e}",
                            imag_res
                        )))
                        .unwrap();
                }
            }
            if let Some(real_res) = self.topologies[0].topo.analytical_result_real {
                if real_res != 0.0 {
                    status_update_sender
                        .send(StatusUpdate::Message(format!(
                            "Target benchmark result (real) : {:+e}",
                            real_res
                        )))
                        .unwrap();
                }
            }
        }
    }

    pub fn multi_channeling<
        'a,
        T: form_integrand::GetIntegrand + form_numerator::GetNumerator + FloatLike,
    >(
        &mut self,
        selected_topology: Option<usize>,
        sample: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCache<T>,
        mut event_manager: Option<&mut EventManager>,
    ) -> Complex<T> {
        let mut rng = thread_rng();
        let mut xrot = [0.; MAX_SG_LOOP * 3];

        // obtain the sampled channel if one is provided
        let (selected_channel, x) = match sample {
            IntegrandSample::Flat(_, x) => (None, x),
            IntegrandSample::Nested(x) => match x {
                havana::Sample::ContinuousGrid(_, x) => (None, x.as_slice()),
                havana::Sample::DiscreteGrid(_, selected_channel, sub_sample) => (
                    Some(selected_channel[0]),
                    match sub_sample.as_ref().unwrap().as_ref() {
                        havana::Sample::ContinuousGrid(_, x) => x.as_slice(),
                        _ => unreachable!(),
                    },
                ),
                havana::Sample::MultiChannel(_, _, _) => unimplemented!(),
            },
        };

        let mut mc_channels = if let Some(t) = selected_topology {
            std::mem::replace(&mut self.topologies[t].multi_channeling_channels, vec![])
        } else {
            std::mem::replace(&mut self.multi_channeling_channels, vec![])
        };

        // paramaterize and consider the result in a channel basis
        let n_fixed = self.settings.cross_section.fixed_cut_momenta.len(); // <- means nfinal-1 is assumed
        let n_loops = x.len() / 3 + n_fixed;
        let mut k_channel = [LorentzVector::default(); MAX_SG_LOOP];
        for i in 0..n_loops {
            xrot[..x.len()].copy_from_slice(x);

            if self.is_stability_check_topo && self.settings.general.stability_nudge_size > 0. {
                for xi in xrot.iter_mut() {
                    if rng.gen_bool(0.5) {
                        *xi += self.settings.general.stability_nudge_size;
                    } else {
                        *xi -= self.settings.general.stability_nudge_size;
                    }

                    *xi = xi.max(0.).min(1.0);
                }
            }
            // parametrize
            let (l_energy, l_space) = if i < n_loops - n_fixed {
                // set the loop index to i + 1 so that we can also shift k
                (
                    T::zero(),
                    Topology::parameterize(
                        &xrot[i * 3..(i + 1) * 3],
                        self.e_cm_squared,
                        i,
                        &self.settings,
                    )
                    .0,
                )
            } else { // FIXED CUT MOMENTA ARE USED HERE nfinal -1
                let m = self.settings.cross_section.fixed_cut_momenta[i + n_fixed - n_loops].cast();
                (m.t, [m.x, m.y, m.z])
            };

            let rot = &self.rotation_matrix;
            k_channel[i] = LorentzVector::from_args(
                l_energy,
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

        let mut k_lmb = [LorentzVector::default(); MAX_SG_LOOP];
        let mut k_other_channel = [LorentzVector::default(); MAX_SG_LOOP];
        let mut event_counter = 0;
        let mut result = Complex::zero();
        for (channel_id, (channel, _, channel_shift)) in mc_channels.iter().enumerate() {
            if let Some(selected_channel) = selected_channel {
                if selected_channel != channel_id {
                    continue;
                }
            }

            if let Some(selected_channel) = self.settings.general.multi_channeling_channel {
                if selected_channel != channel_id as isize {
                    continue;
                }
            }

            // transform to the loop momentum basis
            for (kk, r) in k_lmb[..n_loops].iter_mut().zip_eq(channel.chunks(n_loops)) {
                *kk = LorentzVector::default();

                for ((ko, s), shift) in k_channel[..n_loops]
                    .iter()
                    .zip_eq(r.iter())
                    .zip_eq(channel_shift)
                {
                    *kk += (ko - shift.cast()).multiply_sign(*s);
                }
            }

            // determine the normalization constant
            let mut normalization = T::zero();
            for (_, other_channel_inv, other_channel_shift) in &mc_channels {
                // transform from the loop momentum basis to the other channel basis
                for ((kk, r), shift) in k_other_channel[..n_loops]
                    .iter_mut()
                    .zip_eq(other_channel_inv.chunks(n_loops))
                    .zip_eq(other_channel_shift)
                {
                    *kk = shift.cast();
                    for (ko, s) in k_lmb[..n_loops].iter().zip_eq(r.iter()) {
                        *kk += ko.multiply_sign(*s);
                    }
                }

                let mut inv_jac_para = T::one();
                for i in 0..n_loops {
                    let jac = if i < n_loops - n_fixed {
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

                        Topology::inv_parametrize(&rotated, self.e_cm_squared, i, &self.settings).1
                    } else {
                        (Into::<T>::into(2.) * <T as FloatConst>::PI())
                            .powi(4)
                            .inv()
                    };

                    inv_jac_para *= jac;
                }

                normalization += inv_jac_para;
            }

            // evaluate the integrand in this channel
            let mut channel_result = Complex::zero();
            for (current_supergraph, (t, &m)) in self
                .topologies
                .iter_mut()
                .zip_eq(&self.multiplicity)
                .enumerate()
            {
                if let Some(tt) = selected_topology {
                    if current_supergraph != tt {
                        continue;
                    }
                }

                cache.current_supergraph = current_supergraph;

                if self.settings.general.debug > 0 {
                    println!(
                        "Evaluating supergraph {} for channel {}",
                        t.name, channel_id
                    );

                    print!("Point in channel: ");
                    for i in 0..n_loops {
                        let (x, _) = Topology::inv_parametrize(
                            &k_lmb[i],
                            self.e_cm_squared,
                            i,
                            &self.settings,
                        );
                        if i > 0 {
                            print!(", {}", x.iter().join(","));
                        } else {
                            print!("{}", x.iter().join(","));
                        }
                    }
                    println!("");
                }

                // undo the jacobian for unused dimensions
                let mut jac_correction = T::one();
                for i in t.n_loops..n_loops {
                    let (_, jac) =
                        Topology::inv_parametrize(&k_lmb[i], self.e_cm_squared, i, &self.settings);
                    jac_correction *= jac;
                }

                channel_result += t.evaluate_mom(&k_lmb[..t.n_loops], cache, &mut event_manager)
                    / normalization
                    * jac_correction
                    * Into::<T>::into(m as f64);
                if let Some(em) = &mut event_manager {
                    if em.track_events {
                        for e in em.event_buffer[event_counter..].iter_mut() {
                            e.integrand *= m as f64 / normalization.to_f64().unwrap()
                                * jac_correction.to_f64().unwrap();
                        }
                        event_counter = em.event_buffer.len();
                    }
                }
            }

            result += channel_result;
        }

        if let Some(t) = selected_topology {
            std::mem::swap(
                &mut self.topologies[t].multi_channeling_channels,
                &mut mc_channels,
            );
        } else {
            std::mem::swap(&mut self.multi_channeling_channels, &mut mc_channels);
        };

        result
    }

    pub fn evaluate<
        'a,
        T: form_integrand::GetIntegrand + form_numerator::GetNumerator + FloatLike,
    >(
        &mut self,
        sample: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCache<T>,
        mut event_manager: Option<&mut EventManager>,
    ) -> Complex<T> {
        if self.settings.general.debug > 0 {
            println!("x-space point: {:?}", sample);
        }

        // obtain the sampled topology if one is provided
        let (selected_topology, sub_sample) = if self.topologies.len() > 1 {
            match sample {
                IntegrandSample::Flat(w, x) => (None, IntegrandSample::Flat(w, x)),
                IntegrandSample::Nested(x) => match x {
                    havana::Sample::ContinuousGrid(_, _) => (None, IntegrandSample::Nested(x)),
                    havana::Sample::DiscreteGrid(_, selected_topology, sub_sample) => (
                        Some(selected_topology[0]),
                        IntegrandSample::Nested(sub_sample.as_ref().unwrap().as_ref()),
                    ),
                    havana::Sample::MultiChannel(_, _, _) => unimplemented!(),
                },
            }
        } else {
            (None, sample)
        };

        // clear the deformation cache for non-stability topologies, unless
        // we are integrating an amplitude, since the deformation is fixed
        cache.current_deformation_index = 0;
        if !self.is_stability_check_topo && self.settings.cross_section.fixed_cut_momenta.is_empty()
        {
            cache.deformation_vector_cache.clear();
        }

        if self.settings.general.multi_channeling && !self.multi_channeling_channels.is_empty() {
            return self.multi_channeling(selected_topology, sub_sample, cache, event_manager);
        }

        let x = sub_sample.to_flat();

        let mut rng = thread_rng();
        let mut xrot = [0.; MAX_SG_LOOP * 3];

        // jointly parameterize all squared topologies
        let n_fixed = self.settings.cross_section.fixed_cut_momenta.len();
        let n_loops = x.len() / 3 + n_fixed;
        let mut k = [LorentzVector::default(); MAX_SG_LOOP];
        let mut para_jacs = [T::one(); MAX_SG_LOOP];
        let mut jac_para = T::one();
        for i in 0..n_loops {
            xrot[..x.len()].copy_from_slice(x);

            if self.is_stability_check_topo && self.settings.general.stability_nudge_size > 0. {
                for xi in xrot.iter_mut() {
                    if rng.gen_bool(0.5) {
                        *xi += self.settings.general.stability_nudge_size;
                    } else {
                        *xi -= self.settings.general.stability_nudge_size;
                    }

                    *xi = xi.max(0.).min(1.0);
                }
            }

            let (l_energy, (l_space, jac)) = if i < n_loops - n_fixed {
                // set the loop index to i + 1 so that we can also shift k
                (
                    T::zero(),
                    Topology::parameterize(
                        &xrot[i * 3..(i + 1) * 3],
                        self.e_cm_squared, // NOTE: taking e_cm from the first graph
                        i,
                        &self.settings,
                    ),
                )
            } else {
                let m = self.settings.cross_section.fixed_cut_momenta[i + n_fixed - n_loops].cast();
                (
                    m.t,
                    (
                        [m.x, m.y, m.z],
                        (Into::<T>::into(2.) * <T as FloatConst>::PI()).powi(4),
                    ),
                )
            };

            // there could be some rounding here
            let rot = &self.rotation_matrix;
            k[i] = LorentzVector::from_args(
                l_energy,
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
            para_jacs[i] = jac_para;
        }

        let mut result = Complex::zero();
        let mut event_counter = 0;
        for (current_supergraph, ((t, &m), extra_topos)) in self
            .topologies
            .iter_mut()
            .zip_eq(&self.multiplicity)
            .zip_eq(self.additional_topologies.iter_mut())
            .enumerate()
        {
            if let Some(tt) = selected_topology {
                if current_supergraph != tt {
                    continue;
                }
            }

            cache.current_supergraph = current_supergraph;

            if self.settings.general.debug > 0 {
                println!("Evaluating supergraph {}", t.name);
            }

            // a squared topology may not use all loop variables, so we set the unused jacobians to 1
            let r = t.evaluate_mom(&k[..t.n_loops], cache, &mut event_manager)
                * para_jacs[t.n_loops - 1]
                * Into::<T>::into(m as f64);
            result += r;

            if self
                .settings
                .cross_section
                .compare_with_additional_topologies
            {
                let mut track_setting = false;
                if let Some(em) = event_manager.as_mut() {
                    track_setting = em.track_events;
                    em.track_events = false;
                }
                for (add_id, (topo, mat)) in extra_topos.iter_mut().enumerate() {
                    if self.settings.general.debug > 0 {
                        println!("Evaluating additional topology {}", add_id);
                    }

                    let mut add_ks = [LorentzVector::default(); MAX_SG_LOOP];
                    let external_momenta: ArrayVec<[LorentzVector<T>; MAX_SG_LOOP]> =
                        t.external_momenta.iter().map(|e| e.cast()).collect();

                    for (add_k, row) in add_ks.iter_mut().zip(mat) {
                        *add_k = SquaredTopology::evaluate_signature(
                            row,
                            &external_momenta,
                            &k[..t.n_loops],
                        );
                    }

                    let add_r = topo.evaluate_mom(&add_ks[..t.n_loops], cache, &mut event_manager)
                        * para_jacs[t.n_loops - 1]
                        * Into::<T>::into(m as f64);

                    if r.is_zero() || add_r.is_zero() {
                        if (r - add_r).norm_sqr() > Into::<T>::into(1e-10 * self.e_cm_squared) {
                            println!(
                            "Mismatch for {} between standard and additional topology {}: r={} vs add_r={}", t.name, add_id, r, add_r
                        );
                        }
                    } else if (r - add_r).norm_sqr() / (r.norm_sqr() + add_r.norm_sqr())
                        > Into::<T>::into(1e-10)
                    {
                        println!(
                        "Mismatch for {} between standard and additional topology {}: r={} vs add_r={}",
                        t.name, add_id, r, add_r
                    );
                    }
                }
                if let Some(em) = event_manager.as_mut() {
                    em.track_events = track_setting;
                }
            }

            if let Some(em) = &mut event_manager {
                if em.track_events {
                    for e in em.event_buffer[event_counter..].iter_mut() {
                        e.integrand *= para_jacs[t.n_loops - 1].to_f64().unwrap() * m as f64;
                    }
                    event_counter = em.event_buffer.len();
                }
            }
        }

        result
    }
}

#[derive(Default)]
pub struct SquaredTopologyCache<T: FloatLike> {
    topology_cache: Vec<Vec<Vec<Vec<LTDCache<T>>>>>,
    scalar_products: Vec<T>,
    deformation_vector_cache: Vec<Vec<FixedDeformationLimit>>,
    current_supergraph: usize,
    current_deformation_index: usize,
}

impl<T: FloatLike> SquaredTopologyCache<T> {
    #[inline]
    pub fn get_topology_cache(
        &mut self,
        cut_index: usize,
        diagram_set_index: usize,
        diagram_index: usize,
    ) -> &mut LTDCache<T> {
        &mut self.topology_cache[self.current_supergraph][cut_index][diagram_set_index]
            [diagram_index]
    }
}

/// Cache for squared topology sets.
#[derive(Default)]
pub struct SquaredTopologyCacheCollection {
    pub float_cache: SquaredTopologyCache<float>,
    pub quad_cache: SquaredTopologyCache<f128>,
}

pub trait CachePrecisionSelector<T: FloatLike> {
    fn get(&mut self) -> &mut SquaredTopologyCache<T>;
}

impl CachePrecisionSelector<float> for SquaredTopologyCacheCollection {
    #[inline]
    fn get(&mut self) -> &mut SquaredTopologyCache<float> {
        &mut self.float_cache
    }
}

impl CachePrecisionSelector<f128> for SquaredTopologyCacheCollection {
    #[inline]
    fn get(&mut self) -> &mut SquaredTopologyCache<f128> {
        &mut self.quad_cache
    }
}

impl SquaredTopology {
    // 
    // FOR AMPLITUDES: extensions start
    // 
    pub fn set_external_data(&mut self,input:& SquaredTopologySetInput) {
        

        if input.external_data.is_some() {
            let e_data= input.external_data.clone().unwrap();
            let mut vv = Vec::new();
            for elem in e_data.pol {
                vv.push(from_vec_vec_to_cmplx_lv(elem))
            }
            self.pol.get_or_insert_with(||vv);

            let mut vv = Vec::new();
            for elem in e_data.cpol {
                vv.push(from_vec_vec_to_cmplx_lv(elem))
            }
            self.cpol.get_or_insert_with(||vv);

            let mut vv = Vec::new();
            for elem in e_data.spinor_u {
                vv.push(from_vec_vec_to_cmplx_v(elem))
            }
            self.spinor_u.get_or_insert_with(||vv);

            let mut vv = Vec::new();
            for elem in e_data.spinor_ubar {
                vv.push(from_vec_vec_to_cmplx_v(elem))
            }            
            self.spinor_ubar.get_or_insert_with(||vv);

            let mut vv = Vec::new();
            for elem in e_data.spinor_v {
                vv.push(from_vec_vec_to_cmplx_v(elem))
            }            
            self.spinor_v.get_or_insert_with(||vv);
            

            let mut vv = Vec::new();
            for elem in e_data.spinor_vbar {
                vv.push(from_vec_vec_to_cmplx_v(elem))
            }            
            self.spinor_vbar.get_or_insert_with(||vv);

            self.is_amplitude.get_or_insert_with(||true);
        }

    }

    pub fn overwrite_momenta(&mut self,input:& SquaredTopologySetInput) {
        
        if input.external_data.is_some() {
            let e_data= input.external_data.clone().unwrap();
            let in_momenta = &(e_data.in_momenta.clone());
            let out_momenta = &(e_data.out_momenta.clone());
            // overwrite incoming and e_cm
            self.external_momenta = in_momenta.clone();
            self.external_momenta.extend(in_momenta);
            let mut sum_incoming: lorentz_vector::LorentzVector<f64> = LorentzVector::default();
            for m in in_momenta {
                sum_incoming += *m;
            }
            self.e_cm_squared = sum_incoming.square().abs();
            
            // I get some problem if I try it otherwise: regarding slicing on vectors -> need to ask ben.
            let c_mom = (*out_momenta).clone();
            let mut c_mom_fix = vec![];
            let mut count = 0;
            for cc in c_mom {
                if count < e_data.n_out - 1 {
                    c_mom_fix.push(cc);
                    count +=1;
                }
                
            }
            self.settings.cross_section.incoming_momenta = (*in_momenta).clone();
            self.settings.cross_section.fixed_cut_momenta= c_mom_fix.clone();
            self.default_fixed_cut_momenta = ((*in_momenta).clone(),c_mom_fix);
            self.is_amplitude.get_or_insert_with(||true);

            

        }
        

    }
    // 
    // FOR AMPLITUDES: extensions end
    // 
    pub fn from_file(filename: &str, settings: &Settings) -> Result<SquaredTopology, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open squared topology file {}", filename))
            .suggestion("Does the path exist?")?;

        let mut squared_topo: SquaredTopology = serde_yaml::from_reader(f)
            .wrap_err("Could not parse squared topology file")
            .suggestion("Is it a correct yaml file")?;

        assert!(
            squared_topo.n_loops <= MAX_SG_LOOP,
            "MAX_SG_LOOP is too small: it should be at least {}",
            squared_topo.n_loops
        );

        // update the UV mass and small mass in the squared topology
        // this affects the multi-channeling
        for ll in &mut squared_topo.topo.loop_lines {
            for p in &mut ll.propagators {
                if p.uv {
                    p.m_squared = settings.cross_section.m_uv_sq;
                }
                if p.m_squared == 0. {
                    p.m_squared = settings.cross_section.small_mass_sq;
                }
            }
        }

        squared_topo.settings = settings.clone();
        for cutkosky_cuts in &mut squared_topo.cutkosky_cuts {
            for cut in &mut cutkosky_cuts.cuts {
                if cut.m_squared == 0. {
                    cut.m_squared = settings.cross_section.small_mass_sq;
                }
            }

            for diagram_set in &mut cutkosky_cuts.diagram_sets {
                diagram_set.numerator = if diagram_set.numerator_tensor_coefficients.len() == 0
                    && diagram_set.numerator_tensor_coefficients_sparse.len() == 0
                {
                    LTDNumerator::one(squared_topo.n_loops)
                } else {
                    if !diagram_set.numerator_tensor_coefficients_sparse.is_empty() {
                        LTDNumerator::from_sparse(
                            squared_topo.n_loops,
                            &diagram_set.numerator_tensor_coefficients_sparse,
                        )
                    } else {
                        LTDNumerator::new(
                            squared_topo.n_loops,
                            &diagram_set
                                .numerator_tensor_coefficients
                                .iter()
                                .map(|x| Complex::new(x.0, x.1))
                                .collect::<Vec<_>>(),
                        )
                    }
                };

                for d in &mut diagram_set.diagram_info {
                    d.graph.settings = settings.clone();
                    d.graph.process(false);

                    // update the UV mass
                    for ll in &mut d.graph.loop_lines {
                        for p in &mut ll.propagators {
                            if p.uv {
                                p.m_squared = settings.cross_section.m_uv_sq;
                            }
                            if p.m_squared == 0. {
                                p.m_squared = settings.cross_section.small_mass_sq;
                            }
                        }
                    }
                }
            }
        }

        // overwrite an empty list of fixed cut momenta in the topology settings if there is a default
        if settings.cross_section.fixed_cut_momenta.is_empty() {
            squared_topo.settings.cross_section.fixed_cut_momenta =
                squared_topo.default_fixed_cut_momenta.1.clone();
        }

        // set the external momenta and e_cm
        let incoming_momenta = if settings.cross_section.incoming_momenta.is_empty() {
            &squared_topo.default_fixed_cut_momenta.0
        } else {
            &settings.cross_section.incoming_momenta
        };

        squared_topo.external_momenta = incoming_momenta.clone();
        squared_topo.external_momenta.extend(incoming_momenta);
        let mut sum_incoming: lorentz_vector::LorentzVector<f64> = LorentzVector::default();
        for m in incoming_momenta {
            sum_incoming += *m;
        }
        
        squared_topo.e_cm_squared = sum_incoming.square().abs();
        
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

        let base_path = if settings.cross_section.numerator_source == NumeratorSource::Yaml {
            String::new()
        } else {
            std::env::var("MG_NUMERATOR_PATH")
                .or_else(|_| {
                    let mut pb = Path::new(filename)
                        .parent()
                        .ok_or("no parent")?
                        .parent()
                        .ok_or("no parent")?
                        .to_path_buf();
                    let mut pb_alt = Path::new(filename)
                        .parent()
                        .ok_or("no parent")?
                        .parent()
                        .ok_or("no parent")?
                        .parent()
                        .ok_or("no parent")?
                        .to_path_buf();
                    pb.push("lib");
                    pb_alt.push("lib");
                    if pb.exists() {
                        pb.pop();
                        Ok(pb.into_os_string().into_string().map_err(|_| "bad path")?)
                    } else if pb_alt.exists() {
                        pb_alt.pop();
                        Ok(pb_alt
                            .into_os_string()
                            .into_string()
                            .map_err(|_| "bad path")?)
                    } else {
                        Err("Cannot determine root folder")
                    }
                })
                .expect("Cannot determine base folder from filename. Use MG_NUMERATOR_PATH")
        };

        if settings.cross_section.numerator_source == NumeratorSource::FormIntegrand {
            squared_topo.form_integrand.base_path = base_path.clone();
            squared_topo.form_integrand.form_integrand = Some(form_integrand::load(&base_path));
        }

        if settings.cross_section.numerator_source == NumeratorSource::Form {
            squared_topo.form_numerator.base_path = base_path.clone();
            let diagram_id = squared_topo
                .form_numerator
                .call_signature
                .as_ref()
                .unwrap()
                .id;
            squared_topo.form_numerator.form_numerator = Some(form_numerator::load(&base_path));

            let max_buffer_size = form_numerator::get_buffer_size(
                squared_topo.form_numerator.form_numerator.as_mut().unwrap(),
            );

            for cut in &mut squared_topo.cutkosky_cuts {
                for diagram_set in cut.diagram_sets.iter_mut() {
                    let rank = form_numerator::get_rank(
                        squared_topo.form_numerator.form_numerator.as_mut().unwrap(),
                        diagram_id,
                        diagram_set.id,
                    );

                    // for the FORM numerator, the output is provided in the
                    // basis of the ltd momenta
                    let n_ltd = diagram_set
                        .diagram_info
                        .iter()
                        .map(|diagram_info| diagram_info.graph.n_loops)
                        .sum();

                    if n_ltd > 0 {
                        diagram_set.numerator =
                            LTDNumerator::from_sparse(n_ltd, &[(vec![0; rank], (0., 0.))]);
                    }

                    squared_topo.form_numerator.form_numerator_buffer_size = max_buffer_size;
                    squared_topo.form_numerator.form_numerator_buffer = vec![0.; max_buffer_size];

                    // also give the subgraph numerators the proper size
                    // TODO: set the proper rank of only the variables of the subgraph
                    for diagram_info in &mut diagram_set.diagram_info {
                        if diagram_info.graph.n_loops > 0 {
                            diagram_info.graph.numerator = LTDNumerator::from_sparse(
                                diagram_info.graph.n_loops,
                                &[(vec![0; rank], (0., 0.))],
                            );
                        }
                    }
                }
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

        debug_assert!(signature.1.len() <= external_momenta.len());
        for (&sign, mom) in signature.1.iter().zip(external_momenta) {
            if sign != 0 {
                cut_momentum += mom.multiply_sign(sign);
            }
        }
        cut_momentum
    }

    pub fn generate_multi_channeling_channels(&mut self) {
        // construct the channels from all loop momentum basis (i.e. the LTD cuts)
        self.topo.process(true);

        for mcb in &mut self.multi_channeling_bases {
            mcb.defining_propagators.sort();
        }

        // make sure all except the last of the fixed cut momenta are in the basis
        // to avoid duplicating channels
        let mut fixed_cut_propagators = vec![];
        if self.settings.cross_section.fixed_cut_momenta.len() > 0 {
            // we only have 1 cutkosky cut set and we skip the dependent cut
            'nextcut: for c in
                self.cutkosky_cuts[0].cuts[..self.cutkosky_cuts[0].cuts.len() - 1].iter()
            {
                for (lli, ll) in self.topo.loop_lines.iter().enumerate() {
                    for (pi, p) in ll.propagators.iter().enumerate() {
                        if p.name == c.name {
                            fixed_cut_propagators.push((lli, pi));
                            continue 'nextcut;
                        }
                    }
                }
                println!("Could not find cut in topology: {}", c.name);
            }
            fixed_cut_propagators.sort();
        }

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

                if self.multi_channeling_bases.len() > 0
                    || self.settings.cross_section.fixed_cut_momenta.len() > 0
                {
                    let cut_defining_propagators: Vec<(usize, usize)> = cut
                        .iter()
                        .enumerate()
                        .filter_map(|(i, x)| {
                            if let Cut::PositiveCut(j) | Cut::NegativeCut(j) = x {
                                Some((i, *j))
                            } else {
                                None
                            }
                        })
                        .sorted()
                        .collect();

                    if !fixed_cut_propagators
                        .iter()
                        .all(|c| cut_defining_propagators.contains(c))
                    {
                        continue;
                    }

                    if !self.multi_channeling_bases.iter().any(|mcb| {
                        mcb.defining_propagators
                            .iter()
                            .all(|p| cut_defining_propagators.contains(p))
                    }) {
                        continue;
                    }
                }

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

        self.multi_channeling_channels = multi_channeling_channels;
    }

    pub fn create_caches<T: FloatLike>(&self) -> Vec<Vec<Vec<LTDCache<T>>>> {
        let mut caches = vec![];
        for cutkosky_cuts in &self.cutkosky_cuts {
            let mut lim_cache = vec![];
            for diagram_set in &cutkosky_cuts.diagram_sets {
                let mut diag_cache = vec![];
                for d in &diagram_set.diagram_info {
                    diag_cache.push(LTDCache::<T>::new(&d.graph));
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
        let mut cut_info_cache = [(
            LorentzVector::default(),
            LorentzVector::default(),
            T::zero(),
            T::zero(),
        ); MAX_SG_LOOP + 1];

        // determine an overestimate of the t that solves the energy constraint
        // then -t and t should give two different solutions or do not converge
        let mut t_start = T::zero();
        let mut sum_k = T::zero();
        for (cut_info, cut) in cut_info_cache[..cutkosky_cuts.cuts.len()]
            .iter_mut()
            .zip(&cutkosky_cuts.cuts)
        {
            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            let shift = utils::evaluate_signature(&cut.signature.1, external_momenta);
            let k_norm_sq = k.spatial_squared();
            let k_dot_shift = k.spatial_dot(&shift);
            t_start += Float::abs(k_dot_shift) / k_norm_sq;
            sum_k += k_norm_sq.sqrt();

            *cut_info = (k, shift, k_norm_sq, k_dot_shift);
        }

        t_start += incoming_energy / sum_k;

        // find the two solutions
        let mut solutions = [(T::zero(), T::zero()); 2];
        for (i, &mut mut t) in [-t_start, t_start].iter_mut().enumerate() {
            for it in 0..20 {
                let mut f = incoming_energy;
                let mut df = T::zero();

                for (cut, cut_info) in cutkosky_cuts
                    .cuts
                    .iter()
                    .zip(&cut_info_cache[..cutkosky_cuts.cuts.len()])
                {
                    let energy = ((cut_info.0 * t + cut_info.1).spatial_squared()
                        + Into::<T>::into(cut.m_squared))
                    .sqrt();
                    f -= energy;
                    df -= (t * cut_info.2 + cut_info.3) / energy;
                }

                if debug_level > 4 {
                    println!("  | t{} finder: f={}, df={}, t={}", i, f, df, t);
                }

                if Float::abs(f) < T::epsilon() * Into::<T>::into(10.) * incoming_energy {
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
        if Float::abs(solutions[0].0) + Float::abs(solutions[1].0) == Into::<T>::into(0.0) {
            panic!(
                "Found exact zero solutions: {} for t={} and t={} for k={:?}, ext={:?}",
                solutions[0].0, -t_start, t_start, loop_momenta, external_momenta
            );
        }
        if Float::abs(solutions[0].0 - solutions[1].0)
            / (Float::abs(solutions[0].0) + Float::abs(solutions[1].0))
            < Into::<T>::into(1e-12)
        {
            panic!(
                "Found the same scaling solution twice: {} for t={} and t={} for k={:?}, ext={:?}",
                solutions[0].0, solutions[0].0, solutions[1].0, loop_momenta, external_momenta
            );
        }

        Some(solutions)
    }

    pub fn evaluate_mom<
        T: form_integrand::GetIntegrand + form_numerator::GetNumerator + FloatLike,
    >(
        &mut self,
        loop_momenta: &[LorentzVector<T>],
        cache: &mut SquaredTopologyCache<T>,
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

        let mut external_momenta: ArrayVec<[LorentzVector<T>; MAX_SG_LOOP]> = self
            .external_momenta
            .iter()
            .map(|m| m.map(|c| c.into()))
            .collect();

        let mut cut_momenta = [LorentzVector::default(); MAX_SG_LOOP + 1];
        let mut rescaled_loop_momenta = [LorentzVector::default(); MAX_SG_LOOP];

        let mut subgraph_loop_momenta = [LorentzVector::default(); MAX_SG_LOOP];
        let mut k_def = [LorentzVector::default(); MAX_SG_LOOP];
        let mut result = Complex::zero();
        for cut_index in 0..self.cutkosky_cuts.len() {
            let cutkosky_cuts = &mut self.cutkosky_cuts[cut_index];

            if self.settings.general.debug >= 1 {
                println!(
                    "Cut {}:",
                    cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                );
            }

            let scaling_solutions = if self.settings.cross_section.do_rescaling
                && self.settings.cross_section.fixed_cut_momenta.is_empty()
            {
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

            if scaling_solutions.is_none() || scaling_solutions.unwrap()[1].0 < T::zero() {
                if self.settings.general.debug >= 1 {
                    println!(
                        "Phase space point has no solutions for cut {:?}:",
                        cutkosky_cuts.cuts.iter().map(|c| &c.name).format(", ")
                    );
                }

                if let Some(em) = event_manager {
                    em.no_phase_space_counter += 1;
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
                    &mut k_def[..self.n_loops],
                    cache,
                    event_manager,
                    cut_index,
                    scaling,
                    scaling_jac,
                    None,
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

    // MODIFICATION: Topo level
    pub fn evaluate_cut<
        T: form_integrand::GetIntegrand + form_numerator::GetNumerator + FloatLike,
    >(
        &mut self,
        loop_momenta: &[LorentzVector<T>],
        cut_momenta: &mut [LorentzVector<T>],
        external_momenta: &mut [LorentzVector<T>],
        rescaled_loop_momenta: &mut [LorentzVector<T>],
        subgraph_loop_momenta: &mut [LorentzVector<T>],
        k_def: &mut [LorentzVector<Complex<T>>],
        cache: &mut SquaredTopologyCache<T>,
        event_manager: &mut Option<&mut EventManager>,
        cut_index: usize,
        scaling: T,
        scaling_jac: T,
        selected_diagram_set: Option<usize>,
    ) -> Complex<T> {
        let cutkosky_cuts = &mut self.cutkosky_cuts[cut_index];

        let mut cut_energies_summed = T::zero();
        let mut scaling_result = Complex::one();

        // evaluate the cuts with the proper scaling
        for (cut_mom, cut) in cut_momenta[..cutkosky_cuts.cuts.len()]
            .iter_mut()
            .zip_eq(cutkosky_cuts.cuts.iter())
        {
            let k = utils::evaluate_signature(&cut.signature.0, loop_momenta);
            let shift = utils::evaluate_signature(
                &cut.signature.1,
                &external_momenta[..self.external_momenta.len()],
            );
            *cut_mom = k * scaling + shift;

            if self.settings.cross_section.fixed_cut_momenta.is_empty() {
                let energy = (cut_mom.spatial_squared() + Into::<T>::into(cut.m_squared)).sqrt();
                cut_mom.t = energy.multiply_sign(cut.sign);
                // add (-2 pi i)/(2E)^power for every cut
                scaling_result *=
                    num::Complex::new(T::zero(), -Into::<T>::into(2.0) * <T as FloatConst>::PI())
                        / (Into::<T>::into(2.0) * energy).powi(cut.power as i32);
                cut_energies_summed += energy;
            }
        }

        if self.settings.general.debug >= 1 {
            println!("  | 1/Es = {}", scaling_result);
            println!("  | q0 = {}", cut_energies_summed);
            println!("  | scaling = {}", scaling);
            println!("  | scaling_jac = {}", scaling_jac);
        }

        if self.settings.cross_section.do_rescaling
            && self.settings.cross_section.fixed_cut_momenta.is_empty()
        {
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

                    (
                        (-(Into::<T>::into(
                            self.settings.cross_section.normalising_function.spread,
                        ) * (Float::powi(scaling, 2) + T::one())
                            / scaling))
                            .exp(),
                        if self.settings.cross_section.normalising_function.spread == 1. {
                            Into::<T>::into(0.27973176363304485456919761407082)
                        } else {
                            Into::<T>::into(
                                2. * rgsl::bessel::K1(
                                    2. * self.settings.cross_section.normalising_function.spread,
                                ),
                            )
                        },
                    )
                }
                NormalisingFunction::LeftRightPolynomial => {
                    // Only support center and spread at 1 for now
                    assert_eq!(self.settings.cross_section.normalising_function.center, 1.,
                        "For now the left-right polynomial normalising function only support a center at 1.0.");
                    assert!(self.settings.cross_section.normalising_function.spread>1.,
                        "The left-right polynomial normalising function only support a spread larger than 1.0.");
                    let sigma =
                        Into::<T>::into(self.settings.cross_section.normalising_function.spread);
                    (
                        Float::powf(scaling, sigma)
                            / (Into::<T>::into(1.0)
                                + Float::powf(scaling, Into::<T>::into(2.0) * sigma)),
                        <T as FloatConst>::PI()
                            / (Into::<T>::into(2.0)
                                * sigma
                                * (<T as FloatConst>::PI() / (Into::<T>::into(2.0) * sigma)).cos()),
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
            if self.settings.cross_section.fixed_cut_momenta.is_empty() {
                // set 0th component of external momenta to the sum of the cuts
                // NOTE: we assume 1->1 here and that it is outgoing
                external_momenta[0].t = cut_energies_summed;
                external_momenta[1].t = cut_energies_summed;
            }

            rescaled_loop_momenta[..self.n_loops].copy_from_slice(loop_momenta);
        }

        let mut constants = Complex::one();

        constants *= utils::powi(
            num::Complex::new(
                Into::<T>::into(1.)
                    / <T as Float>::powi(Into::<T>::into(2.) * <T as FloatConst>::PI(), 4),
                T::zero(),
            ),
            self.n_loops,
        );

        if self.settings.cross_section.fixed_cut_momenta.is_empty() {
            // multiply the flux factor
            constants /= if self.n_incoming_momenta == 2 {
                Into::<T>::into(2.)
                    * (SquaredTopology::lambda(
                        (external_momenta[0] + external_momenta[1]).square().abs(),
                        external_momenta[0].square(),
                        external_momenta[1].square(),
                    ))
                    .sqrt()
            } else {
                let e_cm = external_momenta[0].square().sqrt().abs();
                e_cm * Into::<T>::into(2.0)
            };
        }

        if self.settings.cross_section.picobarns {
            // return a weight in picobarns (from GeV^-2)
            constants *= Into::<T>::into(0.389379304e9);
        }

        scaling_result *= constants;

        if self.settings.general.debug >= 2 {
            println!("  | constants={}", constants);
            println!("  | scaling part = {}", scaling_result);
            println!(
                "  | rescaled loop momenta = {:?}",
                &rescaled_loop_momenta[..self.n_loops]
            );
        }

        let e_cm_sq = if self.n_incoming_momenta == 2 {
            (external_momenta[0] + external_momenta[1]).square().abs()
        } else {
            external_momenta[0].square().abs()
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
        // the order is: cut momenta, momenta graph 1, ... graph n
        for (kd, cut_mom) in k_def[..cutkosky_cuts.cuts.len() - 1]
            .iter_mut()
            .zip(&cut_momenta[..cutkosky_cuts.cuts.len() - 1])
        {
            *kd = cut_mom.to_complex(true);
        }
        // PARAMS are read in here. (hyperparams): adjacent reals in memory build complex numbers
        let params = [
            Into::<T>::into(self.settings.cross_section.m_uv_sq.sqrt()),
            T::zero(),
            Into::<T>::into(self.settings.cross_section.mu_r_sq.sqrt()),
            T::zero(),
            Into::<T>::into(self.settings.cross_section.gs),
            T::zero(),
            Into::<T>::into(self.settings.cross_section.small_mass_sq),
            T::zero(),
        ];

        // now apply the same procedure for all uv limits
        let mut diag_and_num_contributions = Complex::zero();
        let mut def_jacobian = Complex::one();

        // regenerate the evaluation of the exponent map of the numerator since the loop momenta have changed
        let mut regenerate_momenta = true;
        for (diag_set_index, diagram_set) in cutkosky_cuts.diagram_sets.iter_mut().enumerate() {
            if self.settings.cross_section.integrand_type == IntegrandType::PF
                && self.settings.cross_section.sum_diagram_sets
                && diag_set_index > 0
            {
                continue;
            }

            if let Some(sid) = selected_diagram_set {
                if sid != diagram_set.id {
                    continue;
                }
            }

            let do_deformation = self.settings.general.deformation_strategy
                == DeformationStrategy::Fixed
                && (diag_set_index == 0
                    || !self
                        .settings
                        .cross_section
                        .inherit_deformation_for_uv_counterterm);

            for diagram_info in &mut diagram_set.diagram_info {
                let subgraph = &mut diagram_info.graph;

                if do_deformation
                    || self.settings.cross_section.numerator_source
                        != NumeratorSource::FormIntegrand
                {
                    // set the shifts, which are expressed in the cut basis
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
                }

                if !do_deformation {
                    subgraph.fixed_deformation.clear();
                    continue;
                }

                subgraph.update_ellipsoids();

                // the stability topologies will inherit the sources
                // amplitudes will also inherit the sources when they are computed once
                if !self.is_stability_check_topo
                    && (cache.current_deformation_index >= cache.deformation_vector_cache.len()
                        || self.settings.cross_section.fixed_cut_momenta.is_empty())
                {
                    subgraph.fixed_deformation =
                        subgraph.determine_ellipsoid_overlap_structure(true);

                    cache
                        .deformation_vector_cache
                        .push(subgraph.fixed_deformation.clone());
                } else {
                    subgraph.fixed_deformation =
                        cache.deformation_vector_cache[cache.current_deformation_index].clone();

                    // rotate the fixed deformation vectors
                    let rot_matrix = &self.rotation_matrix;
                    for d_lim in &mut subgraph.fixed_deformation {
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

                    // update the excluded surfaces
                    for ex in &mut subgraph.all_excluded_surfaces {
                        *ex = false;
                    }

                    for d_lim in &mut subgraph.fixed_deformation {
                        for d in &d_lim.deformation_per_overlap {
                            for surf_id in &d.excluded_surface_indices {
                                subgraph.all_excluded_surfaces[*surf_id] = true;
                            }
                        }
                    }
                }
                cache.current_deformation_index += 1;

                if self.settings.general.debug > 0 {
                    // check if the overlap structure makes sense
                    subgraph.check_fixed_deformation();
                }
            }

            if !self
                .settings
                .cross_section
                .inherit_deformation_for_uv_counterterm
                || diag_set_index == 0
            {
                def_jacobian = Complex::one();
            }

            // compute the deformation vectors
            let mut k_def_index = cutkosky_cuts.cuts.len() - 1;
            for (diagram_index, diagram_info) in diagram_set.diagram_info.iter_mut().enumerate() {
                let subgraph = &diagram_info.graph;
                let subgraph_cache =
                    cache.get_topology_cache(cut_index, diag_set_index, diagram_index);
                if !self
                    .settings
                    .cross_section
                    .inherit_deformation_for_uv_counterterm
                    || diag_set_index == 0
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

                    let (kappas, jac_def) = if self.settings.general.deformation_strategy
                        == DeformationStrategy::Fixed
                    {
                        subgraph.deform(&subgraph_loop_momenta[..subgraph.n_loops], subgraph_cache)
                    } else {
                        (
                            [LorentzVector::default(); crate::MAX_LOOP],
                            Complex::new(T::one(), T::zero()),
                        )
                    };

                    for (lm, kappa) in subgraph_loop_momenta[..subgraph.n_loops]
                        .iter()
                        .zip(&kappas[..subgraph.n_loops])
                    {
                        k_def[k_def_index] = if diagram_info.conjugate_deformation {
                            // take the complex conjugate of the deformation
                            lm.map(|x| Complex::new(x, T::zero()))
                                - kappa.map(|x| Complex::new(T::zero(), x))
                        } else {
                            lm.map(|x| Complex::new(x, T::zero()))
                                + kappa.map(|x| Complex::new(T::zero(), x))
                        };
                        k_def_index += 1;
                    }

                    if diagram_info.conjugate_deformation {
                        def_jacobian *= jac_def.conj();
                    } else {
                        def_jacobian *= jac_def;
                    }
                } else {
                    k_def_index += subgraph.n_loops;
                }

                if self.settings.cross_section.numerator_source != NumeratorSource::FormIntegrand {
                    if subgraph
                        .populate_ltd_cache(
                            &k_def[k_def_index - subgraph.n_loops..k_def_index],
                            subgraph_cache,
                        )
                        .is_err()
                    {
                        if self.settings.deformation.fixed.ir_handling_strategy == IRHandling::None
                        {
                            panic!("NaN on cut energy");
                        }
                    }
                }

                subgraph_cache.cached_topology_integrand.clear();
            }
            // Setup for evaluation of FORM code
            if self.settings.cross_section.numerator_source == NumeratorSource::Form
                || self.settings.cross_section.numerator_source == NumeratorSource::FormIntegrand
            {
                if !self
                    .settings
                    .cross_section
                    .inherit_deformation_for_uv_counterterm
                    || diag_set_index == 0
                {
                    cache.scalar_products.clear();

                    for (i1, e1) in external_momenta[..self.n_incoming_momenta]
                        .iter()
                        .enumerate()
                    {
                        cache.scalar_products.extend_from_slice(&[e1.t, T::zero()]);
                        for e2 in &external_momenta[i1..self.n_incoming_momenta] {
                            let (d, ds) = e1.dot_spatial_dot(e2);
                            cache
                                .scalar_products
                                .extend_from_slice(&[d, T::zero(), ds, T::zero()]);
                        }
                    }

                    for (i1, m1) in k_def[..self.n_loops].iter().enumerate() {
                        cache.scalar_products.extend_from_slice(&[m1.t.re, m1.t.im]);
                        for e1 in &external_momenta[..self.n_incoming_momenta] {
                            let (d, ds) = m1.dot_spatial_dot(&e1.cast());
                            cache
                                .scalar_products
                                .extend_from_slice(&[d.re, d.im, ds.re, ds.im]);
                        }

                        for m2 in k_def[i1..self.n_loops].iter() {
                            let (d, ds) = m1.dot_spatial_dot(m2);
                            cache
                                .scalar_products
                                .extend_from_slice(&[d.re, d.im, ds.re, ds.im]);
                        }
                    }
                    // ADDITIONAL ENTRIES FROM AMPLITUDES                                            
                    if self.is_amplitude.is_some() {
                        println!("I arrived but dont know what to do :(")

                    }

                }

                if self.settings.cross_section.numerator_source == NumeratorSource::Form {
                    let mut form_numerator =
                        mem::replace(&mut self.form_numerator.form_numerator, None);
                    let mut form_numerator_buffer =
                        mem::replace(&mut self.form_numerator.form_numerator_buffer, vec![]);
                    if let Some(call_signature) = &self.form_numerator.call_signature {
                        let mut form_numerator_buffer: ArrayVec<[T; 100]> =
                            (0..2 * self.form_numerator.form_numerator_buffer_size)
                                .map(|_| T::zero())
                                .collect();

                        let time_start = Instant::now();
                        let len = T::get_numerator(
                            form_numerator.as_mut().unwrap(),
                            &cache.scalar_products,
                            &params,
                            call_signature.id,
                            diagram_set.id,
                            &mut form_numerator_buffer,
                        );

                        if let Some(em) = event_manager {
                            em.integrand_evaluation_timing +=
                                Instant::now().duration_since(time_start).as_nanos();
                        }

                        let first_subgraph_cache =
                            cache.get_topology_cache(cut_index, diag_set_index, 0);

                        // the output of the FORM numerator is already in the reduced lb format
                        first_subgraph_cache.reduced_coefficient_lb[0].resize(len, Complex::zero());
                        for (rlb, r) in first_subgraph_cache.reduced_coefficient_lb[0]
                            .iter_mut()
                            .zip_eq(form_numerator_buffer[..2 * len].chunks(2))
                        {
                            *rlb = Complex::new(Into::<T>::into(r[0]), Into::<T>::into(r[1]));
                        }
                    } else {
                        panic!(
                        "No call signature for FORM numerator, but FORM numerator mode is enabled"
                    );
                    }

                    mem::swap(&mut form_numerator, &mut self.form_numerator.form_numerator);
                    mem::swap(
                        &mut form_numerator_buffer,
                        &mut self.form_numerator.form_numerator_buffer,
                    );
                } else {
                    let mut form_integrand =
                        mem::replace(&mut self.form_integrand.form_integrand, None);

                    if let Some(em) = event_manager {
                        if em.time_integrand_evaluation {
                            em.integrand_evaluation_timing_start = Some(Instant::now());
                        }
                    }

                    let mut res = if let Some(call_signature) = &self.form_integrand.call_signature
                    {
                        let res = match self.settings.cross_section.integrand_type {
                            IntegrandType::LTD => T::get_integrand_ltd(
                                form_integrand.as_mut().unwrap(),
                                &cache.scalar_products,
                                &params,
                                call_signature.id,
                                diagram_set.id,
                            ),
                            IntegrandType::PF => T::get_integrand_pf(
                                form_integrand.as_mut().unwrap(),
                                &cache.scalar_products,
                                &params,
                                call_signature.id,
                                if self.settings.cross_section.sum_diagram_sets {
                                    1000 + cut_index
                                } else {
                                    diagram_set.id
                                },
                            ),
                        };

                        Complex::new(Into::<T>::into(res.re), Into::<T>::into(res.im))
                    } else {
                        panic!(
                        "No call signature for FORM integrand, but FORM integrand mode is enabled"
                    );
                    };

                    if let Some(em) = event_manager {
                        if let Some(s) = em.integrand_evaluation_timing_start.take() {
                            em.integrand_evaluation_timing +=
                                Instant::now().duration_since(s).as_nanos();
                        }
                    }

                    mem::swap(&mut form_integrand, &mut self.form_integrand.form_integrand);

                    res *= def_jacobian;

                    if self.settings.general.debug >= 1 {
                        println!("  | res diagram set {}: = {:e}", diagram_set.id, res);
                    }

                    diag_and_num_contributions += res;
                    continue;
                }
            }

            let first_subgraph_cache = cache.get_topology_cache(cut_index, diag_set_index, 0);

            if self.settings.cross_section.numerator_source == NumeratorSource::Yaml {
                // evaluate the cut numerator with the spatial parts of all loop momenta and
                // the energy parts of the Cutkosky cuts.
                // Store the reduced numerator in the left graph cache for now
                // TODO: this is impractical
                diagram_set.numerator.evaluate_reduced_in_lb(
                    &k_def,
                    cutkosky_cuts.cuts.len() - 1,
                    first_subgraph_cache,
                    0,
                    regenerate_momenta,
                    false,
                );
                regenerate_momenta = true; // false;
            }

            mem::swap(
                &mut first_subgraph_cache.reduced_coefficient_lb_supergraph[0],
                &mut first_subgraph_cache.reduced_coefficient_lb[0],
            );

            let mut supergraph_coeff = mem::replace(
                &mut first_subgraph_cache.reduced_coefficient_lb_supergraph[0],
                vec![],
            );

            // evaluate the subgraphs for every monomial in the numerator
            let mut result_complete_numerator = Complex::default();
            for (coeff, powers) in supergraph_coeff.iter().zip(
                &diagram_set.numerator.reduced_coefficient_index_to_powers
                    [..diagram_set.numerator.reduced_size],
            ) {
                if coeff.is_zero() {
                    continue;
                }

                // the yaml gives the numerator in the cb basis, the other methods in the LTD basis, consisting
                // of only the LTD momenta
                let mut ltd_index =
                    if self.settings.cross_section.numerator_source == NumeratorSource::Yaml {
                        cutkosky_cuts.cuts.len() - 1
                    } else {
                        0
                    };

                // only consider the coefficients that have no powers in the cutkosky cuts
                // TODO: make a more efficient way of skipping the other contributions
                assert!(powers[..ltd_index].iter().all(|p| *p == 0));

                if self.settings.general.debug >= 1 {
                    println!("  | monomial {:?} = {}", powers, coeff);
                }

                let mut num_result = *coeff;
                for (diagram_index, diagram_set) in diagram_set.diagram_info.iter_mut().enumerate()
                {
                    let subgraph = &mut diagram_set.graph;
                    let subgraph_cache =
                        cache.get_topology_cache(cut_index, diag_set_index, diagram_index);
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
                        let mut subgraph_powers = [0; 10]; // TODO: make max rank a constant
                        let mut power_index = 0;

                        for (lmi, p) in powers[ltd_index..ltd_index + subgraph.n_loops]
                            .iter()
                            .enumerate()
                        {
                            for _ in 0..*p {
                                subgraph_powers[power_index] = lmi * 4;
                                power_index += 1;
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

                    // Store coefficients in a multi variable polynomial
                    if self.settings.general.partial_fractioning_threshold > 0.0 {
                        subgraph_cache.reduced_coefficient_lb_mpoly.clear();
                        for (c, pows) in subgraph_cache.reduced_coefficient_lb[0].iter().zip(
                            subgraph.numerator.reduced_coefficient_index_to_powers
                                [..subgraph.numerator.reduced_size]
                                .iter(),
                        ) {
                            //println!("#{:?} -> {}", &pows[..subgraph.n_loops], c);
                            subgraph_cache
                                .reduced_coefficient_lb_mpoly
                                .add(&pows[..subgraph.n_loops], *c);
                            //println!(
                            //    "\t->pows: {:?},\n\t->coeffs: {:?},\n\t->n_var: {}",
                            //    subgraph_cache.reduced_coefficient_lb_mpoly.powers,
                            //    subgraph_cache.reduced_coefficient_lb_mpoly.coeffs,
                            //    subgraph_cache.reduced_coefficient_lb_mpoly.n_var
                            //);
                        }
                    }
                    let mut res = if cached_res.is_none() {
                        let res = if subgraph.loop_lines.len() > 0 {
                            subgraph
                                .evaluate_all_dual_integrands::<T>(
                                    if subgraph.n_loops != 0 {
                                        if self.settings.cross_section.numerator_source
                                            == NumeratorSource::Yaml
                                        {
                                            &mut k_def[ltd_index..ltd_index + subgraph.n_loops]
                                        } else {
                                            let start = cutkosky_cuts.cuts.len() - 1 + ltd_index;
                                            &mut k_def[start..start + subgraph.n_loops]
                                        }
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
                            "  | monomial {} res {} ({}l) = {:e}",
                            monomial_index, subgraph.name, subgraph.n_loops, res
                        );
                    }

                    ltd_index += subgraph.n_loops;
                }

                if self.settings.general.debug >= 1 {
                    println!("  | monomial result = {:e}", num_result * def_jacobian);
                    println!(
                        "  | monomial result full weight = {:e}",
                        num_result * def_jacobian * scaling_result
                    )
                }

                result_complete_numerator += num_result * def_jacobian;
            }

            let first_subgraph_cache = cache.get_topology_cache(cut_index, diag_set_index, 0);
            mem::swap(
                &mut supergraph_coeff,
                &mut first_subgraph_cache.reduced_coefficient_lb_supergraph[0],
            );

            if self.settings.general.debug >= 1 {
                println!(
                    "  | res diagram set {}: = {:e}",
                    diagram_set.id, result_complete_numerator
                );
            }

            diag_and_num_contributions += result_complete_numerator;
        }

        scaling_result *= diag_and_num_contributions;

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
        rotated_topology.name += "_rot";
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
            for uv_limit in cut.diagram_sets.iter_mut() {
                uv_limit.numerator = uv_limit.numerator.rotate(rot_matrix);
            }
        }

        rotated_topology
    }
}

impl IntegrandImplementation for SquaredTopologySet {
    type Cache = SquaredTopologyCacheCollection;

    fn create_stability_check(&self, num_checks: usize) -> Vec<SquaredTopologySet> {
        let mut stability_topologies = vec![];

        let mut rng = rand::thread_rng();
        for i in 0..num_checks {
            if self.stability_check_topologies.iter().all(|x| i < x.len()) {
                let mut sts = SquaredTopologySet {
                    name: self.name.clone(),
                    e_cm_squared: self.e_cm_squared,
                    topologies: self
                        .stability_check_topologies
                        .iter()
                        .map(|sct| sct[i].clone())
                        .collect(),
                    additional_topologies: vec![vec![]; self.topologies.len()],
                    multiplicity: self.multiplicity.clone(),
                    settings: self.settings.clone(),
                    rotation_matrix: self.rotation_matrix.clone(),
                    multi_channeling_channels: vec![],
                    stability_check_topologies: vec![vec![]; self.topologies.len()],
                    is_stability_check_topo: true,    
                };
                sts.create_multi_channeling_channels();
                stability_topologies.push(sts);
            } else {
                // we don't have enough stability topologies, so we pad with rotations
                let angle =
                    float::from_f64(rng.gen::<f64>() * 2.).unwrap() * <float as FloatConst>::PI();
                let mut rv = (
                    float::from_f64(rng.gen()).unwrap(),
                    float::from_f64(rng.gen()).unwrap(),
                    float::from_f64(rng.gen()).unwrap(),
                ); // rotation axis
                let inv_norm = (rv.0 * rv.0 + rv.1 * rv.1 + rv.2 * rv.2).sqrt().inv();
                rv = (rv.0 * inv_norm, rv.1 * inv_norm, rv.2 * inv_norm);

                stability_topologies.push(self.rotate(angle, rv))
            }
        }
        stability_topologies
    }

    fn get_target(&self) -> Option<Complex<f64>> {
        let mut target = Complex::zero();

        for (t, m) in self.topologies.iter().zip(&self.multiplicity) {
            if t.topo.analytical_result_real.is_some() || t.topo.analytical_result_imag.is_some() {
                target += Complex::new(
                    t.topo.analytical_result_real.unwrap_or(0.),
                    t.topo.analytical_result_imag.unwrap_or(0.),
                ) * *m;
            }
        }

        if !target.is_zero() {
            Some(target)
        } else {
            None
        }
    }

    fn set_partial_fractioning(&mut self, enable: bool) {
        for t in self.topologies.iter_mut() {
            if enable {
                t.settings.cross_section.integrand_type = IntegrandType::PF;
            } else {
                t.settings.cross_section.integrand_type = IntegrandType::LTD;
            }

            if self.settings.cross_section.numerator_source == NumeratorSource::Form {
                for cs in &mut t.cutkosky_cuts {
                    for ds in &mut cs.diagram_sets {
                        for di in &mut ds.diagram_info {
                            if enable {
                                di.graph.settings.general.partial_fractioning_threshold = 1e-99;
                            } else {
                                di.graph.settings.general.partial_fractioning_threshold = -1.;
                            }
                        }
                    }
                }
            }
        }
    }

    fn create_grid(&self) -> Grid {
        if self.topologies.len() > 1 {
            Grid::DiscreteGrid(DiscreteGrid::new(
                &[self.topologies.len()],
                self.topologies
                    .iter()
                    .map(|t| {
                        if self.settings.general.multi_channeling
                            && !t.multi_channeling_channels.is_empty()
                            && self.settings.general.multi_channeling_channel.is_none()
                        {
                            // construct the discrete grid for multi-channeling
                            // TODO: now we can replace the loop counts by the actual loop counts needed per topology
                            Grid::DiscreteGrid(DiscreteGrid::new(
                                &[t.multi_channeling_channels.len()],
                                (0..t.multi_channeling_channels.len())
                                    .map(|_| {
                                        Grid::ContinuousGrid(ContinuousGrid::new(
                                            3 * t.n_loops,
                                            self.settings.integrator.n_bins,
                                            self.settings.integrator.min_samples_for_update,
                                        ))
                                    })
                                    .collect(),
                                self.settings.integrator.min_probability_per_bin,
                            ))
                        } else {
                            Grid::ContinuousGrid(ContinuousGrid::new(
                                3 * self.get_maximum_loop_count(), // t.n_loops,
                                self.settings.integrator.n_bins,
                                self.settings.integrator.min_samples_for_update,
                            ))
                        }
                    })
                    .collect(),
                self.settings.integrator.min_probability_per_bin,
            ))
        } else {
            if self.settings.general.multi_channeling
                && !self.multi_channeling_channels.is_empty()
                && self.settings.general.multi_channeling_channel.is_none()
            {
                // construct the discrete grid for multi-channeling
                // TODO: now we can replace the loop counts by the actual loop counts needed per topology
                Grid::DiscreteGrid(DiscreteGrid::new(
                    &[self.multi_channeling_channels.len()],
                    (0..self.multi_channeling_channels.len())
                        .map(|_| {
                            Grid::ContinuousGrid(ContinuousGrid::new(
                                3 * self.get_maximum_loop_count(),
                                self.settings.integrator.n_bins,
                                self.settings.integrator.min_samples_for_update,
                            ))
                        })
                        .collect(),
                    self.settings.integrator.min_probability_per_bin,
                ))
            } else {
                Grid::ContinuousGrid(ContinuousGrid::new(
                    3 * self.get_maximum_loop_count(),
                    self.settings.integrator.n_bins,
                    self.settings.integrator.min_samples_for_update,
                ))
            }
        }
    }

    #[inline]
    fn evaluate_float<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCacheCollection,
        event_manager: Option<&mut EventManager>,
    ) -> Complex<float> {
        self.evaluate(x, cache.get(), event_manager)
    }

    #[inline]
    fn evaluate_f128<'a>(
        &mut self,
        x: IntegrandSample<'a>,
        cache: &mut SquaredTopologyCacheCollection,
        event_manager: Option<&mut EventManager>,
    ) -> Complex<f128> {
        self.evaluate(x, cache.get(), event_manager)
    }

    #[inline]
    fn create_cache(&self) -> SquaredTopologyCacheCollection {
        SquaredTopologyCacheCollection {
            float_cache: SquaredTopologyCache {
                topology_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
                deformation_vector_cache: vec![],
                scalar_products: vec![],
                current_supergraph: 0,
                current_deformation_index: 0,
            },
            quad_cache: SquaredTopologyCache {
                topology_cache: self.topologies.iter().map(|t| t.create_caches()).collect(),
                deformation_vector_cache: vec![],
                scalar_products: vec![],
                current_supergraph: 0,
                current_deformation_index: 0,
            },
        }
    }
}
