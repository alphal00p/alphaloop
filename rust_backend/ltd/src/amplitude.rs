use arrayvec::ArrayVec;
use gamma_chain::GammaChain;
use itertools::Itertools;
use num::Complex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use topologies::{Cut, LTDCache, LTDNumerator, Topology};
use utils;
use vector::{Field, LorentzVector};
use FloatLike;
use GeneralSettings;

#[allow(
    non_snake_case,
    non_upper_case_globals,
    non_camel_case_types,
    dead_code
)]
pub mod Parameters {
    pub const g_f: f64 = 1.166390e-05;
    pub const alpha_s: f64 = 1.180000e-01;
    pub const alpha_ew: f64 = 1. / 1.325070e+02;
    pub const C_F: f64 = 4. / 3.;
    pub const q_u: f64 = 2. / 3.;
    pub const q_d: f64 = -1. / 3.;
}

#[allow(non_camel_case_types, dead_code)]
#[derive(Debug, Copy, Deserialize, Clone)]
pub enum Polarizations {
    #[serde(rename = "u+")]
    UPlus,
    #[serde(rename = "u-")]
    UMinus,
    #[serde(rename = "ubar+")]
    UBarPlus,
    #[serde(rename = "ubar-")]
    UBarMinus,
    #[serde(rename = "v+")]
    VPlus,
    #[serde(rename = "v-")]
    VMinus,
    #[serde(rename = "vbar+")]
    VBarPlus,
    #[serde(rename = "vbar-")]
    VBarMinus,
    #[serde(rename = "a+")]
    APlus,
    #[serde(rename = "a-")]
    AMinus,
}

impl Polarizations {
    pub fn complement(&self) -> Polarizations {
        match self {
            Polarizations::UPlus => Polarizations::VMinus,
            Polarizations::UMinus => Polarizations::VPlus,
            Polarizations::UBarPlus => Polarizations::VBarMinus,
            Polarizations::UBarMinus => Polarizations::VBarPlus,
            Polarizations::VPlus => Polarizations::UMinus,
            Polarizations::VMinus => Polarizations::UPlus,
            Polarizations::VBarPlus => Polarizations::UBarMinus,
            Polarizations::VBarMinus => Polarizations::UBarPlus,
            Polarizations::APlus => Polarizations::AMinus,
            Polarizations::AMinus => Polarizations::APlus,
        }
    }
}

// Giving Serde some love
#[derive(Serialize, Deserialize)]
#[serde(remote = "Complex")]
pub struct ComplexDef<T> {
    pub re: T,
    pub im: T,
}
// This should be enough love

#[derive(Default, Debug, Clone, Deserialize)]
pub struct DiagramFullRust {
    pub name: String,
    denominators: Vec<(usize, usize)>,
    pows: Vec<usize>,
    chain: Vec<i8>,
    positions: Vec<i8>,
    loop_signature: i8,
    #[serde(with = "ComplexDef", rename = "factor")]
    factor_f64: Complex<f64>,
    ct: bool,
    tensor_coefficients_split: Vec<[f64; 2]>,
    #[serde(skip_deserializing)]
    tensor_coefficients: Vec<Complex<f64>>,
    #[serde(skip_deserializing)]
    pub numerator: LTDNumerator,
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct SlashedMomentum<T: Field> {
    pub loop_mom: Vec<i8>,
    //#[serde(with = "LorentzVectorDef")]
    pub const_mom: LorentzVector<T>,
}

impl SlashedMomentum<Complex<f64>> {
    fn evaluate<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<Complex<T>>],
    ) -> LorentzVector<Complex<T>> {
        let mut res = LorentzVector::default();
        if self.loop_mom.len() > 0 {
            for (l, c) in loop_momenta.iter().zip(self.loop_mom.iter()) {
                res += l * Complex::new(T::from_i8(*c).unwrap(), T::zero());
            }
            res += self
                .const_mom
                .map(|x| Complex::new(Into::<T>::into(x.re), Into::<T>::into(x.im)));
        } else {
            res += self
                .const_mom
                .map(|x| Complex::new(Into::<T>::into(x.re), Into::<T>::into(x.im)));
        }
        res
    }
}

impl<T: Field> std::fmt::Display for SlashedMomentum<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "loop_mom: {:?},\nconst_mom: {}",
            self.loop_mom, self.const_mom
        )
    }
}

#[derive(Default, Debug, Clone, Deserialize)]
pub struct Amplitude {
    pub name: String,
    pub amp_type: String,
    pub topology: String,
    pub n_loops: usize,
    pub sets: Vec<Vec<String>>,
    #[serde(rename = "vectors")]
    pub vectors_f64: Vec<SlashedMomentum<f64>>,
    #[serde(skip_deserializing)]
    pub vectors: Vec<SlashedMomentum<Complex<f64>>>,
    pub pols_type: Vec<Polarizations>,
    pub ps: Vec<LorentzVector<f64>>,
    pub mu_r_sq: f64,
    #[serde(skip_deserializing)]
    pub born: Complex<f64>,
    #[serde(skip_deserializing)]
    pub int_ct: Vec<Complex<f64>>,
    #[serde(skip_deserializing)]
    pub polarizations: Vec<ArrayVec<[Complex<f64>; 4]>>,
    pub diagrams: Vec<DiagramFullRust>,
    #[serde(skip_deserializing)]
    pub tensor_coefficient_index_map: Vec<(usize, usize)>,
}

// Implement Setup functions
impl Amplitude {
    pub fn from_file(filename: &str) -> HashMap<String, Amplitude> {
        let f = File::open(filename).expect("Could not open amplitude file");
        let amplitudes: Vec<Amplitude> = serde_yaml::from_reader(f).unwrap();
        amplitudes
            .into_iter()
            .map(|a| (a.name.clone(), a))
            .collect()
    }

    pub fn process(&mut self, settings: &GeneralSettings) {
        //Compute the needed polarizations
        self.polarizations = Vec::new();
        for (pol, p) in self.pols_type.iter().zip(self.ps.iter()) {
            self.polarizations.push(if p[0] < 0.0 {
                match pol {
                    Polarizations::AMinus | Polarizations::APlus => {
                        compute_polarization(*p, pol.complement()).unwrap()
                    }
                    _ => compute_polarization(-*p, pol.complement())
                        .unwrap()
                        .iter()
                        .map(|x| -x)
                        .collect(),
                }
            } else {
                compute_polarization(*p, *pol).unwrap()
            })
        }

        if settings.debug > 2 {
            for (i, (_pol, p)) in self.pols_type.iter().zip(self.ps.iter()).enumerate() {
                println!("  | p{} = {:?}", i, p);
            }
            for (i, (pol, p)) in self
                .pols_type
                .iter()
                .zip(self.polarizations.iter())
                .enumerate()
            {
                print!("  | (p{}, {:?}) \t= [", i, pol);
                for x in p.iter() {
                    print!(" {:+.5e}", x);
                }
                println!("]");
            }
        }

        // Convert coefficient to Complex<f64>
        self.diagrams[0].tensor_coefficients = Vec::new();
        for dn in 0..self.diagrams.len() {
            let mut diag = &mut self.diagrams[dn];
            diag.tensor_coefficients = Vec::new();
            for cn in 0..diag.tensor_coefficients_split.len() {
                diag.tensor_coefficients.push(Complex::new(
                    diag.tensor_coefficients_split[cn][0],
                    diag.tensor_coefficients_split[dn][1],
                ));
            }
        }

        //Attach the polarizations to vectors (based on the amplitude)
        match self.amp_type.as_ref() {
            "qqbarphotonsNLO" => {
                // Convert the input coming from vectorsf64 to allow complex vectors
                self.vectors = Vec::new();
                for slashed in self.vectors_f64.iter() {
                    self.vectors.push(SlashedMomentum {
                        loop_mom: slashed.loop_mom.clone(),
                        const_mom: slashed.const_mom.map(|x| Complex::new(x, 0.0)),
                    });
                }
                //Attach polarizaitons
                for p in self.polarizations[2..].iter() {
                    self.vectors.push(SlashedMomentum {
                        loop_mom: Vec::default(),
                        const_mom: LorentzVector::from_slice(&p),
                    });
                }
                // Create Coefficients
                self.chain_to_coefficients().unwrap();
                for diag in self.diagrams.iter_mut() {
                    diag.numerator = LTDNumerator::new(self.n_loops, &diag.tensor_coefficients);
                }

                // Print vector elements
                if settings.debug > 2 {
                    for v in self.vectors.iter() {
                        println!("{}", v);
                    }
                }
                //Compute the integrated counterterms
                self.integrated_ct();
            }
            "DiHiggsTopologyLO" => {
                self.int_ct = vec![Complex::new(0.0, 0.0)];
                for diag in self.diagrams.iter_mut() {
                    diag.numerator = LTDNumerator::new(self.n_loops, &diag.tensor_coefficients);
                }
            }
            _ => panic!("Unknown amplitude type: {}", self.amp_type),
        };

        // TODO: Evaluate with new numerator
        self.construct_numerator();
        println!("DONE");
    }

    pub fn construct_numerator(&mut self) {
        let n_loops = self.n_loops;
        // Determine max_rank
        let mut max_rank = 0;
        let mut check_size = 1;
        let mut level_size = 1;
        for diag in self.diagrams.iter() {
            while check_size < diag.tensor_coefficients_split.len() {
                level_size = (level_size * (n_loops * 4 + max_rank)) / (max_rank + 1);
                max_rank += 1;
                check_size += level_size;
            }
            // Check if the sizes match
            if max_rank != 0 {
                assert_eq!(diag.tensor_coefficients_split.len(), check_size);
            }
        }
        let sorted_linear: Vec<Vec<usize>> = (0..max_rank + 1)
            .map(|rank| (0..4 * n_loops).combinations_with_replacement(rank))
            .flatten()
            .collect();
        self.tensor_coefficient_index_map = sorted_linear
            .iter()
            .map(|l| {
                if l.len() == 1 {
                    (0, l[l.len() - 1]) // note: these indices will never be used
                } else {
                    (
                        sorted_linear
                            .iter()
                            .position(|x| x[..] == l[..l.len() - 1])
                            .unwrap(),
                        l[l.len() - 1],
                    )
                }
            })
            .collect();
    }
}

// Implement Evaluation Functions
impl Amplitude {
    pub fn compute_chain<T: FloatLike>(
        &self,
        indices: &[i8],
        vbar: &[Complex<T>],
        u: &[Complex<T>],
        vectors: &[LorentzVector<Complex<T>>],
    ) -> Result<Complex<T>, &'static str> {
        GammaChain::new(vbar, u, indices, vectors)
            .unwrap()
            .compute_chain()
    }

    pub fn chain_to_coefficients(&mut self) -> Result<(), &'static str> {
        // Create the polarization vecotors
        let vbar: ArrayVec<[Complex<f64>; 4]> = self.polarizations[1]
            .iter()
            .map(|x| Complex::new(x.re, x.im))
            .collect();
        let u: ArrayVec<[Complex<f64>; 4]> = self.polarizations[0]
            .iter()
            .map(|x| Complex::new(x.re, x.im))
            .collect();

        let e_mu = vec![
            LorentzVector::from_args(
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ),
            LorentzVector::from_args(
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ),
            LorentzVector::from_args(
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ),
            LorentzVector::from_args(
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ),
        ];

        // Evaluate the vectors with loop momentum 0
        let mut vectors: ArrayVec<[LorentzVector<Complex<f64>>; 32]> = self
            .vectors
            .iter()
            .map(|v| v.evaluate(&[LorentzVector::default()]))
            .collect();
        // Check that there are enough slots available in self.vector
        // to allow an easy evaluation of the coefficients
        let mut count = 0;
        for v in self.vectors.iter() {
            if v.loop_mom.len() > 0 {
                count += 1;
            }
        }
        assert!(
            count
                >= self
                    .diagrams
                    .iter()
                    .map(|d| d.positions.len())
                    .max()
                    .unwrap()
        );

        // Go through all the diagrams to translate the the gamma_chain into a list of coefficients
        for diag_and_cut in self.diagrams.iter_mut() {
            // Global factor
            let factor = Complex::new(diag_and_cut.factor_f64.re, diag_and_cut.factor_f64.im);
            // Get chain info
            let position = diag_and_cut.positions.to_vec();
            //Find all unique permutations coming form each elements
            //of combinations_with_repetitions
            diag_and_cut.tensor_coefficients = Vec::new();

            // If 4 set the loop momentum to zero
            for (coeff_n, coeff_mus) in (0..=4)
                .combinations_with_replacement(position.len())
                .enumerate()
            {
                diag_and_cut.tensor_coefficients.push(Complex::default());
                // Sum all over all the contributions for this coefficient
                for permuted_mus in coeff_mus.iter().permutations(coeff_mus.len()).unique() {
                    // Change the values store in vectors to evaluate the right coefficient
                    let mut index = 0;
                    let mut chain = diag_and_cut.chain.to_vec();
                    for (&mu, &pos) in permuted_mus.iter().zip_eq(position.iter()) {
                        //Update Vector
                        let old_pos = chain[pos as usize] as usize - 1;
                        vectors[index] = if *mu != 0 {
                            match self.vectors[old_pos].loop_mom[0] {
                                -1 => -e_mu[*mu - 1],
                                1 => e_mu[*mu - 1],
                                _ => panic!(
                                    "Unknown option while changing from chain to coefficients"
                                ),
                            }
                        //self.vectors[old_pos].evaluate(&[e_mu[*mu - 1]])
                        } else {
                            self.vectors[old_pos].evaluate(&[LorentzVector::default()])
                        };
                        //Update chain
                        index += 1;
                        chain[pos as usize] = index as i8;
                    }
                    diag_and_cut.tensor_coefficients[coeff_n] += factor
                        * GammaChain::new(
                            vbar.as_slice(),
                            u.as_slice(),
                            &chain,
                            vectors.as_slice(),
                        )
                        .unwrap()
                        .compute_chain()?;
                }
            }
        }
        Ok(())
    }

    pub fn compute_chain_amplitude<T: FloatLike>(
        &self,
        propagators: &HashMap<(usize, usize), Complex<T>>,
        vectors: &[LorentzVector<Complex<T>>],
        cut_2energy: Complex<T>,
        cut_ll_id: Vec<(usize, usize)>,
        settings: &GeneralSettings,
        e_cm_sq: T,
    ) -> Result<Complex<T>, &'static str> {
        let mut res: Complex<T> = Complex::default();

        // Select which set of diagrams to consider
        // sets1: with UV approximation for all the diagrams
        // sets0: regular amplitude
        let diaglist = if cut_2energy.norm() > e_cm_sq * Into::<T>::into(settings.mu_uv_sq_re_im[0])
        {
            self.sets[1].clone()
        } else {
            self.sets[0].clone()
        };
        if settings.debug >= 2 {
            println!("Set of diagrams: {:?}", diaglist);
        }

        //TODO: avoid converting them all the time
        let vbar: ArrayVec<[Complex<T>; 4]> = self.polarizations[1]
            .iter()
            .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
            .collect();
        let u: ArrayVec<[Complex<T>; 4]> = self.polarizations[0]
            .iter()
            .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
            .collect();

        for diag_and_cut in self.diagrams.iter() {
            if cut_ll_id
                .iter()
                .all(|cut_id| diag_and_cut.denominators.iter().any(|v| *v == *cut_id))
                & diaglist.iter().any(|v| v == &diag_and_cut.name)
            {
                let res0 = res;
                let chain = diag_and_cut.chain.to_vec();
                //Compute denominator
                let mut diag_den = Complex::new(T::one(), T::zero());
                for cut in diag_and_cut.denominators.iter() {
                    diag_den *= propagators[cut];
                }
                //Compute numerator
                let factor = Complex::new(
                    T::from_f64(diag_and_cut.factor_f64.re).unwrap(),
                    T::from_f64(diag_and_cut.factor_f64.im).unwrap(),
                );

                // TODO: extend the derivatie to mulit-loop
                let cut_id = cut_ll_id[0];
                let pow_position = diag_and_cut
                    .denominators
                    .iter()
                    .position(|&x| x == cut_id)
                    .unwrap();
                let diag_num = self.uv_residue(
                    &chain,
                    &diag_and_cut.positions,
                    vectors,
                    vbar.as_slice(),
                    u.as_slice(),
                    diag_and_cut.pows[pow_position] - 1,
                    cut_2energy,
                    factor,
                );
                //println!("{} Numerator = {}", diag_and_cut.name, diag_num);
                res += match diag_and_cut.ct {
                    true => -diag_num * utils::finv(diag_den),
                    false => diag_num * utils::finv(diag_den),
                };
                //println!("{} = {}", diag_and_cut.name, res - res0);
                if settings.debug > 2 {
                    println!(
                        "  |cut[{:?}]: {} \t=> {:e}",
                        cut_id,
                        diag_and_cut.name,
                        res - res0
                    );
                }
            }
        }
        return Ok(res);
    }

    pub fn compute_coefficient_amplitude<T: FloatLike>(
        &self,
        topo: &Topology,
        loop_momenta: &mut [LorentzVector<Complex<T>>],
        cut_2energy: Complex<T>,
        cut_ll_id: Vec<(usize, usize)>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &'static str> {
        let mut res: Complex<T> = Complex::default();

        // Select which set of diagrams to consider
        // sets1: with UV approximation for all the diagrams
        // sets0: regular amplitude
        let diaglist = if cut_2energy.norm()
            > Into::<T>::into(topo.e_cm_squared * topo.settings.general.mu_uv_sq_re_im[0])
        {
            self.sets[1].clone()
        } else {
            self.sets[0].clone()
        };
        if topo.settings.general.debug >= 2 {
            println!("Set of diagrams: {:?}", diaglist);
        }
        for diag_and_cut in self.diagrams.iter() {
            //Filter to get only diagrams that are cut
            if diaglist.iter().any(|v| v == &diag_and_cut.name) {
                let res0 = res;
                //Whole evaluation is taken care by Topology::evaluate_cut once the correct
                //powers for the propagators are set and the corresponding numerator is selected
                for prop_pow in cache.propagator_powers.iter_mut() {
                    *prop_pow = 0;
                }
                for (ll_id, &pow) in diag_and_cut
                    .denominators
                    .iter()
                    .zip_eq(diag_and_cut.pows.iter())
                {
                    cache.propagator_powers[topo.loop_lines[ll_id.0].propagators[ll_id.1].id] = pow;
                }
                //Compute diagram on cut
                res += match diag_and_cut.ct {
                    true => -topo.evaluate_cut(
                        loop_momenta,
                        &diag_and_cut.numerator,
                        cut,
                        mat,
                        cache,
                        false,
                    )?,
                    false => topo.evaluate_cut(
                        loop_momenta,
                        &diag_and_cut.numerator,
                        cut,
                        mat,
                        cache,
                        false,
                    )?,
                };
                if topo.settings.general.debug >= 2 {
                    println!(
                        "  |cut[{:?}]: {} \t=> {:e}",
                        cut_ll_id[0],
                        diag_and_cut.name,
                        res - res0
                    );
                }
            }
        }
        return Ok(res);
    }

    pub fn uv_residue<T: FloatLike>(
        &self,
        chain: &[i8],
        loop_momentum_positions: &[i8],
        vectors: &[LorentzVector<Complex<T>>],
        vbar: &[Complex<T>],
        u: &[Complex<T>],
        order: usize,
        cut_2energy: Complex<T>,
        factor: Complex<T>,
    ) -> Complex<T> {
        //TODO: use try_into
        let gamma_0_pos = vectors.len() as i8;
        let mut coeff = 1;
        let mut norm = 1;
        //The sign variable depends only on order because we choose a negative k
        let sign = Complex::new(T::one(), T::zero());
        let mut res = Complex::default();
        for i in (0..order + 1).rev() {
            //println!("coeff: {:?}", coeff);
            //println!(
            //    "[{},{}, {:?}]",
            //    i,
            //    2 * order - i,
            //    utils::powi(sign, order - i).re
            //);
            res += self
                .compute_chain_residue(
                    vbar,
                    u,
                    chain,
                    vectors,
                    loop_momentum_positions,
                    gamma_0_pos,
                    i,
                )
                .unwrap()
                * utils::powi(sign, order - i)
                * T::from_usize(coeff).unwrap()
                * utils::finv(utils::powi(cut_2energy, 2 * order - i));
            //Update coefficients
            if i != 0 {
                coeff = (coeff * i * (2 * order - i + 1)) / (order - i + 1);
                norm *= i;
            }
        }
        if order % 2 == 0 {
            factor * res / T::from_usize(norm).unwrap()
        } else {
            -factor * res / T::from_usize(norm).unwrap()
        }
    }
    //Compute the n-th derivative w.r.t. k0 giving the corresponding gamma chain list
    //Needs the gamma0 position inside vectors.
    pub fn compute_chain_residue<T: FloatLike>(
        &self,
        vbar: &[Complex<T>],
        u: &[Complex<T>],
        indices: &[i8],
        vectors: &[LorentzVector<Complex<T>>],
        loop_momentum_positions: &[i8],
        gamma_0_pos: i8,
        order: usize,
    ) -> Result<Complex<T>, &'static str> {
        if order == 0 {
            return GammaChain::new(vbar, u, indices, vectors)
                .unwrap()
                .compute_chain();
        }
        let mut res = Complex::default();
        let factor = T::from_usize(order).unwrap();
        let mut left_indices = loop_momentum_positions.to_vec();
        for (i, pos) in loop_momentum_positions.iter().enumerate() {
            if i > order {
                break;
            }
            left_indices.remove(0);
            let mut diff_indices = indices.to_vec();
            diff_indices[*pos as usize] = gamma_0_pos;
            res += self
                .compute_chain_residue(
                    vbar,
                    u,
                    &diff_indices,
                    vectors,
                    &left_indices,
                    gamma_0_pos,
                    order - 1,
                )
                .unwrap()
                * factor;
        }
        Ok(res)
    }

    pub fn evaluate<T: FloatLike>(
        &self,
        propagators: &HashMap<(usize, usize), Complex<T>>,
        loop_momenta: &[LorentzVector<Complex<T>>],
        cut_2energy: Complex<T>,
        cut_ll_id: Vec<(usize, usize)>,
        e_cm_sq: T,
        settings: &GeneralSettings,
    ) -> Result<Complex<T>, &'static str> {
        let gamma_0 = LorentzVector::from_args(
            Complex::new(T::one(), T::zero()),
            Complex::new(T::zero(), T::zero()),
            Complex::new(T::zero(), T::zero()),
            Complex::new(T::zero(), T::zero()),
        );

        //Collinear to a1
        //let loop_momenta = [(-self.ps[0]-self.ps[1]+self.ps[2].map(|x| 0.5* x)).map(|x| Complex::new(Into::<T>::into(x),T::zero()))];
        //println!("loop_mom = {}",loop_momenta[0]);
        //Get the incormation about the slashed vectors
        let mut vectors: ArrayVec<[LorentzVector<Complex<T>>; 32]> = self
            .vectors
            .iter()
            .map(|v| v.evaluate(&loop_momenta.to_vec()))
            .collect();
        vectors.push(gamma_0);
        self.compute_chain_amplitude(
            propagators,
            &vectors.to_vec(),
            cut_2energy,
            cut_ll_id,
            settings,
            e_cm_sq,
        )
    }
    pub fn evaluate_with_coefficients<T: FloatLike>(
        &self,
        topo: &Topology,
        loop_momenta: &mut [LorentzVector<Complex<T>>],
        cut_2energy: Complex<T>,
        cut_ll_id: Vec<(usize, usize)>,
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &'static str> {
        // Update tensor loop dependent part
        self.update_numerator_momentum(loop_momenta, cache);

        self.compute_coefficient_amplitude(
            topo,
            loop_momenta,
            cut_2energy,
            cut_ll_id,
            cut,
            mat,
            cache,
        )
    }
    //In debug mode the integrated counterterms can be seen
    pub fn integrated_ct(&mut self) {
        match self.amp_type.as_ref() {
            "qqbarphotonsNLO" => {
                //Get the information about the slashed vectors
                let loop_mom = vec![LorentzVector::default()];
                let vectors: ArrayVec<[LorentzVector<Complex<f64>>; 32]> =
                    self.vectors.iter().map(|v| v.evaluate(&loop_mom)).collect();

                //Spinors
                let vbar = self.polarizations[1].clone();
                let u = self.polarizations[0].clone();

                // Compute Born
                let mut res: Complex<f64> = Complex::default();
                for diag_and_cut in self.diagrams.iter() {
                    if diag_and_cut.name == "born" {
                        //Compute gamma chain
                        let chain = GammaChain::new(
                            vbar.as_slice(),
                            u.as_slice(),
                            &diag_and_cut.chain,
                            vectors.as_slice(),
                        )
                        .unwrap()
                        .compute_chain()
                        .unwrap();
                        res = diag_and_cut.factor_f64 * chain;
                    }
                }
                self.born = res;

                // Get right expression from the logarithm
                let mu_r = self.mu_r_sq;
                let s12 = (self.ps[0] + self.ps[1]).square();
                let ln = if s12 < 0.0 {
                    Complex::new((-mu_r * mu_r / s12).ln(), 0.0)
                } else {
                    Complex::new((mu_r * mu_r / s12).ln(), std::f64::consts::PI)
                };
                // Get integrated CT poles and finite part
                let epm2 = -self.born * 2.0 * Parameters::C_F / (4.0 * std::f64::consts::PI)
                    * Parameters::alpha_s;
                let epm1 = epm2 * (ln + 1.5);
                let ep0 = epm2 * (ln * (ln + 3.) * 0.5 + 4.);
                //Store the integrated counterterms
                self.int_ct = vec![ep0, epm1, epm2];
            }
            _ => panic!("Unknown integrated counterterm for {}.", self.amp_type),
        }
    }
    //Print out the result for the integrated counterterms
    pub fn print_integrated_ct(&self) {
        match self.amp_type.as_ref() {
            "qqbarphotonsNLO" => {
                // Print polarizations
                for (i, (pol, p)) in self
                    .pols_type
                    .iter()
                    .zip(self.polarizations.iter())
                    .enumerate()
                {
                    //println!("  | p{}        \t= {:.5e} ", i, self.ps[i]);
                    print!("  | (p{}, {:?}) \t= [", i, pol);
                    for x in p.iter() {
                        print!(" {:+.5e}", x);
                    }
                    println!("]");
                }
                // Print CT
                println!("\n  | born    \t: {:.5e}", self.born);
                for (i, val) in self.int_ct.iter().enumerate() {
                    match i {
                        0 => println!("  | Finite  \t: {:+.5e}", val),
                        1 => println!("  | 1/ep    \t: {:+.5e}", val),
                        _ => println!("  | 1/ep**{} \t: {:+.5e}", i, val),
                    };
                }
            }
            _ => panic!("Unknown integrated counterterm for {}.", self.amp_type),
        }
    }

    /* cache numerator coefficients */

    pub fn update_numerator_momentum<T: FloatLike>(
        &self,
        loop_momenta: &[LorentzVector<Complex<T>>],
        cache: &mut LTDCache<T>,
    ) -> () {
        if cache.numerator_momentum_cache.len() < self.tensor_coefficient_index_map.len() + 1 {
            println!(
                "Resizing numerator_momentum_cache {} -> {}",
                cache.numerator_momentum_cache.len(),
                self.tensor_coefficient_index_map.len() + 1
            );
            cache.numerator_momentum_cache.resize(
                self.tensor_coefficient_index_map.len() + 1,
                Complex::new(T::zero(), T::zero()),
            );
        }
        for (i, &(cache_index, vec_index)) in self.tensor_coefficient_index_map.iter().enumerate() {
            cache.numerator_momentum_cache[i + 1] = if cache_index == 0 {
                loop_momenta[vec_index / 4][vec_index % 4]
            } else {
                cache.numerator_momentum_cache[cache_index]
                    * loop_momenta[vec_index / 4][vec_index % 4]
            };
        }
    }
    pub fn evaluate_numerator<T: FloatLike>(
        &self,
        coefficients: &[Complex<f64>],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        let mut coefficient = Complex::new(
            Into::<T>::into(coefficients[0].re),
            Into::<T>::into(coefficients[0].im),
        );
        let mut result = coefficient;
        for (i, &c) in coefficients.iter().skip(1).enumerate() {
            coefficient = Complex::new(Into::<T>::into(c.re), Into::<T>::into(c.im));
            result += cache.numerator_momentum_cache[i + 1] * coefficient;
        }
        result
    }
}

impl Topology {
    pub fn evaluate_amplitude_cut<T: FloatLike>(
        &self,
        k_def: &mut [LorentzVector<Complex<T>>],
        cut: &Vec<Cut>,
        mat: &Vec<i8>,
        cache: &mut LTDCache<T>,
    ) -> Result<Complex<T>, &'static str> {
        // Ensure that an amplitude is defined
        assert!(self.settings.general.use_amplitude);

        // Update cache information about the propagators
        // The cut propagator will be replaced with the corresponding energy residue

        // Update loop momenta
        self.set_loop_momentum_energies(k_def, cut, mat, cache);
        // compute propagators
        for (n, ll) in self.loop_lines.iter().enumerate() {
            // Build loopline loop momentum part from signature
            let mut ll_ks: LorentzVector<Complex<T>> = LorentzVector::new();
            for (loop_mom_n, sign) in ll.signature.iter().enumerate() {
                ll_ks += match sign {
                    1 => k_def[loop_mom_n],
                    -1 => -k_def[loop_mom_n],
                    _ => panic!("Unknown loopline signature option {:?}", ll.signature),
                };
            }

            for (m, p) in ll.propagators.iter().enumerate() {
                cache.propagators.insert(
                    (n, m),
                    utils::powi(ll_ks.t + T::from_f64(p.q.t).unwrap(), 2)
                        - cache.complex_prop_spatial[p.id],
                );
                //println!("Prop[{:?}] = {:?}", (n, m), props[&(n, m)])
            }
        }
        // compute residue energy
        let mut cut_2energy = Complex::new(T::one(), T::zero());
        let mut cut_id;
        let mut cut_ll_id = Vec::new();
        for (n, (ll_cut, ll)) in cut.iter().zip(self.loop_lines.iter()).enumerate() {
            // get the loop line result from the cache
            match ll_cut {
                Cut::PositiveCut(j) | Cut::NegativeCut(j) => {
                    cut_id = ll.propagators[*j].id;
                    cut_ll_id.push((n, *j));
                    cut_2energy = cache.complex_cut_energies[cut_id] * Into::<T>::into(2.);
                    //Replace the cut propagator by its 2E
                    //println!("old {:?} = {:?}", (n, j), props.get(&(n, *j)));
                    cache.propagators.insert((n, *j), cut_2energy);
                    //println!("new {:?} = {:?}", (n, j), props.get(&(n, *j)));
                }
                _ => continue,
            };
        }
        // Evaluate Amplitude
        match self.amplitude.amp_type.as_ref() {
            "qqbarphotonsNLO" => {
                if true {
                    self.amplitude.evaluate_with_coefficients(
                        &self,
                        k_def,
                        cut_2energy,
                        cut_ll_id,
                        cut,
                        mat,
                        cache,
                    )
                } else {
                    //Rotate back PS for numerator
                    let rot_matrix = self.rotation_matrix;
                    let mut l = k_def[0];
                    let old_x = l.x;
                    let old_y = l.y;
                    let old_z = l.z;
                    l.x = old_x * T::from_f64(rot_matrix[0][0]).unwrap()
                        + old_y * T::from_f64(rot_matrix[1][0]).unwrap()
                        + old_z * T::from_f64(rot_matrix[2][0]).unwrap();
                    l.y = old_x * T::from_f64(rot_matrix[0][1]).unwrap()
                        + old_y * T::from_f64(rot_matrix[1][1]).unwrap()
                        + old_z * T::from_f64(rot_matrix[2][1]).unwrap();
                    l.z = old_x * T::from_f64(rot_matrix[0][2]).unwrap()
                        + old_y * T::from_f64(rot_matrix[1][2]).unwrap()
                        + old_z * T::from_f64(rot_matrix[2][2]).unwrap();
                    let l_def = vec![l];
                    //Evaluate Amplitude
                    self.amplitude.evaluate(
                        &cache.propagators,
                        &l_def,
                        cut_2energy,
                        cut_ll_id,
                        Into::<T>::into(self.e_cm_squared),
                        &self.settings.general,
                    )
                }
            }
            "DiHiggsTopologyLO" => {
                //Evaluate Amplitude
                self.amplitude.evaluate_with_coefficients(
                    &self,
                    k_def,
                    cut_2energy,
                    cut_ll_id,
                    cut,
                    mat,
                    cache,
                )
            }
            _ => {
                panic!("Unknown amplitude type: {}", self.amplitude.amp_type);
            }
        }
    }
}

#[allow(unused_variables)]
fn compute_polarization<T: FloatLike>(
    p: LorentzVector<f64>,
    polarization: Polarizations,
) -> Result<ArrayVec<[Complex<T>; 4]>, &'static str> {
    //Only for massless case
    //Spherical coordinates of the spatial part of p
    let rho = p.spatial_squared().sqrt();
    let theta = (p[3] / rho).acos();
    let phi = if p[2].abs() < p[3].abs() * 1e-10 && p[1].abs() < p[3].abs() * 1e-10 && p[3] > 0.0 {
        0.0
    } else if p[2].abs() < p[3].abs() * 1e-10 && p[1].abs() < p[3].abs() * 1e-10 && p[3] < 0.0 {
        std::f64::consts::PI
    } else {
        //p.arctan2(p[2], p[1])
        p[2].atan2(p[1])
    };
    let phase = (Complex::new(0.0, 1.0) * phi).exp();

    // define gamma matrices in Weyl representation,
    let gamma0 = [
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
    ];
    let gamma5 = [
        [-1., 0., 0., 0.],
        [0., -1., 0., 0.],
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
    ];

    //HELAS convention for V_+- to U_-+ conversion
    // E > 0 : v_plus = -u_minus, v_minus = -u_plus
    // E < 0 : v_plus = u_minus, v_minus = u_plus
    let factor = match polarization {
        Polarizations::VPlus
        | Polarizations::VBarPlus
        | Polarizations::VMinus
        | Polarizations::VBarMinus
            if p[0] > 0. =>
        {
            -1.0
        }
        _ => 1.0,
    };

    let pol = match polarization {
        Polarizations::UPlus
        | Polarizations::UBarPlus
        | Polarizations::VMinus
        | Polarizations::VBarMinus => {
            // spinors as documented in HELAS reference, rewritten in Dirac representation and in polar coordinates
            let mut u_plus: Vec<Complex<f64>> = if p[0] > 0.0 {
                vec![
                    Complex::default(),
                    Complex::default(),
                    Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                    (1.0 - (theta).cos()).sqrt() * phase,
                ]
            } else {
                vec![
                    Complex::new(0.0, (1.0 + (theta).cos()).sqrt()),
                    Complex::new(0.0, 1.0) * (1.0 - (theta).cos()).sqrt() * phase,
                    Complex::default(),
                    Complex::default(),
                ]
            }
            .iter()
            .map(|x| factor * x * (rho).sqrt())
            .collect();
            // this line adjust sign to HELAS convention
            if p[2].abs() < p[3].abs() * 1e-10 && p[1].abs() < p[3].abs() * 1e-10 && p[3] < 0.0 {
                u_plus = u_plus.iter().map(|x| -x).collect();
            }
            match polarization {
                Polarizations::UPlus | Polarizations::VMinus => u_plus,
                _ => {
                    let mut u_bar_plus = Vec::new();
                    for g in gamma0.iter() {
                        let mut res = Complex::default();
                        for (i, x) in u_plus.iter().enumerate() {
                            res += x.conj() * g[i];
                        }
                        u_bar_plus.push(res);
                    }
                    u_bar_plus
                }
            }
        }
        Polarizations::UMinus
        | Polarizations::UBarMinus
        | Polarizations::VPlus
        | Polarizations::VBarPlus => {
            // spinors as documented in HELAS reference, rewritten in Dirac representation and in polar coordinates
            // according to HELAS convention, v_plus = -u_minus, v_minus = -u_plus
            let mut u_minus: Vec<Complex<f64>> = if p[0] > 0.0 {
                vec![
                    -(1.0 - (theta).cos()).sqrt() / phase,
                    Complex::new((1.0 + (theta).cos()).sqrt(), 0.0),
                    Complex::default(),
                    Complex::default(),
                ]
            } else {
                vec![
                    Complex::default(),
                    Complex::default(),
                    -Complex::new(0.0, 1.0) * (1.0 - (theta).cos()).sqrt() / phase,
                    Complex::new(0.0, (1.0 + (theta).cos()).sqrt()),
                ]
            }
            .iter()
            .map(|x| factor * x * (rho).sqrt())
            .collect();
            // this line adjust sign to HELAS convention
            if p[2].abs() < p[3].abs() * 1e-10 && p[1].abs() < p[3].abs() * 1e-10 && p[3] < 0.0 {
                u_minus = u_minus.iter().map(|x| -x).collect();
            }
            match polarization {
                Polarizations::UMinus | Polarizations::VPlus => u_minus,
                _ => {
                    let mut u_bar_minus = Vec::new();
                    for g in gamma0.iter() {
                        let mut res = Complex::default();
                        for (i, x) in u_minus.iter().enumerate() {
                            res += x.conj() * g[i];
                        }
                        u_bar_minus.push(res);
                    }
                    u_bar_minus
                }
            }
        }
        Polarizations::APlus => {
            // photon polarization vectors
            let a1 = vec![
                0.,
                (theta).cos() * (phi).cos(),
                (theta).cos() * (phi).sin(),
                -(theta).sin(),
            ];
            let a2 = vec![0., -(phi).sin(), (phi).cos(), 0.];

            //helicity eigenstates of photons
            let ii = Complex::new(0.0, 1.0);
            let norm = 1.0 / (2.0 as f64).sqrt();
            let mut a_plus = Vec::new();
            for i in 0..a1.len() {
                a_plus.push(norm * (-a1[i] - ii * a2[i]));
            }

            a_plus
        }
        Polarizations::AMinus => {
            // photon polarization vectors
            let a1 = vec![
                0.,
                (theta).cos() * (phi).cos(),
                (theta).cos() * (phi).sin(),
                -(theta).sin(),
            ];
            let a2 = vec![0., -(phi).sin(), (phi).cos(), 0.];
            //helicity eigenstates of photons
            let ii = Complex::new(0.0, 1.0);
            let norm = 1.0 / (2.0 as f64).sqrt();
            let mut a_minus = Vec::new();
            for i in 0..a1.len() {
                a_minus.push(norm * (a1[i] - ii * a2[i]));
            }

            a_minus
        }
    };
    Ok(pol
        .iter()
        .map(|x| Complex::new(T::from_f64(x.re).unwrap(), T::from_f64(x.im).unwrap()))
        .collect())
}
//END compute_polarization
