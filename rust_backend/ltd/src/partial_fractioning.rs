use itertools::Itertools;
use num::Complex;
use topologies::{LTDCache, LTDNumerator, LoopLine};
use FloatLike;

#[derive(Default, Debug, Clone)]
pub struct PartialFractioning {
    pub partial_fractioning_element: Vec<PartialFractioningMonomial>,
    pub n_props_deg: usize,
    pub numerator_rank: usize,
    //splits: Vec<usize>,
}

#[derive(Default, Debug, Clone)]
pub struct PFCache {
    splits: Vec<usize>,
    ellipsoids_product: Vec<(usize, usize)>,
    numerator: Vec<usize>,
    numerator_size: usize,
    numerator_index_map: Vec<usize>,
}

impl PFCache {
    pub fn new(n_prop: usize) -> PFCache {
        if n_prop == 0 {
            PFCache {
                splits: vec![],
                ellipsoids_product: vec![],
                numerator: vec![],
                numerator_size: 0,
                numerator_index_map: vec![],
            }
        } else {
            // Use this to avoid redundant allocation
            PFCache {
                splits: vec![0; n_prop - 1],
                ellipsoids_product: vec![(0, 0); n_prop - 1],
                numerator: vec![0; n_prop],
                numerator_size: 0,
                numerator_index_map: vec![0; n_prop],
            }
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningMonomial {
    pub ellipsoids_product: Vec<(usize, usize)>,
    pub numerator: Vec<usize>,
}

/// Loop through all subset of lenth n of the set of
/// all integers from 0 to set_size-1
fn get_next_subset(subset: &mut [usize], set_size: usize, bottom_up: bool) -> bool {
    let n = subset.len();
    // Move subset forward
    if bottom_up {
        for i in 0..n {
            if i + 1 == n || subset[i] + 1 < subset[i + 1] {
                subset[i] += 1;
                for j in 0..i {
                    subset[j] = j;
                }
                return subset[n - 1] < set_size;
            }
        }
    } else {
        // Move subset backward
        for i in 0..n {
            if subset[i] > i {
                subset[i] -= 1;
                for j in 1..=i {
                    subset[i - j] = subset[i] - j
                }
                return subset[n - 1] < set_size;
            }
        }
    }
    return false; // TODO: Use a panic
}

fn get_next_set_pair(subset1: &mut [usize], subset2: &mut [usize], set_size: usize) -> bool {
    // Move subset1 forward and subset2 backward
    // such that if subset1 U subset2 == set remains true
    return get_next_subset(subset1, set_size, true) && get_next_subset(subset2, set_size, false);
}

impl PartialFractioningMonomial {
    pub fn new(ellipsoids: &[(usize, usize)], num: &[usize]) -> PartialFractioningMonomial {
        let numerator = num.iter().map(|x| *x).collect();
        let ellipsoids_product = ellipsoids.iter().map(|x| *x).collect();

        return PartialFractioningMonomial {
            ellipsoids_product: ellipsoids_product,
            numerator: numerator,
        };
    }

    pub fn evaluate<T: FloatLike>(
        &self,
        numerator: &LTDNumerator,
        loop_line: &[LoopLine], // make sure to only pass loop lines with a loop momentum in it
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let den = PartialFractioningMonomial::evaluate_ellipsoids_product(
            &self.ellipsoids_product,
            ltd_cache,
        );
        let num = PartialFractioningMonomial::evaluate_numerator(
            self.numerator.as_slice(),
            numerator,
            loop_line,
            ltd_cache,
        );

        return den.inv() * num;
    }

    pub fn evaluate_ellipsoids_product<T: FloatLike>(
        ellipsoids_product: &[(usize, usize)],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let mut res = Complex::new(T::one(), T::zero());
        for (id1, id2) in ellipsoids_product.iter() {
            // If surface not present skip it
            if ltd_cache.propagator_powers[*id1] == 0 || ltd_cache.propagator_powers[*id2] == 0 {
                panic!("Ellipsoid do not exists!");
            }
            // Check if single power propagators
            //if ltd_cache.propagator_powers[*id1] != 1 || ltd_cache.propagator_powers[*id2] != 1 {
            //    panic!("Propagators have power higher then 1 and partial fractioning still has to be defined!");
            //}
            res *= ltd_cache.complex_ellipsoids[*id1][*id2];
        }
        return res;
    }

    pub fn evaluate_numerator<T: FloatLike>(
        pf_numerator: &[usize],
        ltd_numerator: &LTDNumerator,
        ll: &[LoopLine],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        //return Complex::new(T::one(), T::zero());
        if ltd_numerator.max_rank + 1 < pf_numerator.len() {
            Complex::default()
        } else {
            PartialFractioningMonomial::evaluate_numerator2(
                pf_numerator,
                ltd_numerator.max_rank + 1 - pf_numerator.len(),
                0,
                ltd_numerator,
                ll,
                ltd_cache,
            )
        }
    }

    fn evaluate_numerator2<T: FloatLike>(
        pf_numerator: &[usize],
        rank: usize,
        offset: usize,
        ltd_numerator: &LTDNumerator,
        ll: &[LoopLine],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let mut res = Complex::default();

        let prop = ll[0]
            .propagators
            .iter()
            .filter(|p| p.id == pf_numerator[0])
            .next()
            .unwrap();
        let x = ltd_cache.complex_cut_energies[prop.id] - Into::<T>::into(prop.q.t);

        if pf_numerator.len() == 1 {
            for n in offset..=rank {
                res += ltd_cache.reduced_coefficient_lb[0]
                    //[ltd_numerator.max_rank - n]// Only at 1-Loop
                    [ltd_numerator.reduced_powers_to_position(&[n as u8,0,0,0])]
                    * x.powi((n - offset) as i32);
            }
            return res;
        }
        // TODO: Multiplication by the ltd_cache.... is complete bollocks
        // Waiting to have the right thing
        for n in offset..=rank {
            res -= PartialFractioningMonomial::evaluate_numerator2(
                &pf_numerator[1..pf_numerator.len()],
                rank + 1,
                n + 1,
                ltd_numerator,
                ll,
                ltd_cache,
            ) * x.powi((n - offset) as i32);
        }

        return res;
    }
}

impl PartialFractioning {
    // Create the partial fractioninig expression for a one-loop LTD integrand
    pub fn new(n_prop: usize, numerator_rank: usize) -> PartialFractioning {
        {
            // Create container
            let mut pf_expr = PartialFractioning {
                partial_fractioning_element: Vec::new(),
                n_props_deg: n_prop,
                numerator_rank: numerator_rank,
            };

            // Use this to avoid redundant allocation
            let mut pf_cache = PFCache::new(n_prop);

            let mut lh = vec![0; n_prop];
            let mut le = vec![0; n_prop];

            // Loop trough all the generating elemnets
            for ne in 0..n_prop {
                let nh = n_prop - ne - 1;
                for j in 0..=nh {
                    lh[j] = j;
                }
                for j in 0..ne {
                    le[j] = j + nh + 1;
                }
                //println!("{:?}[:{}] {:?}[:{}]", lh, nh + 1, le, ne);
                pf_cache.numerator_size = 0;
                pf_expr.element_partial_fractioning(&lh[0..=nh], &le[0..ne], &mut pf_cache);
                while get_next_set_pair(&mut lh[0..=nh], &mut le[0..ne], n_prop) {
                    //println!(" -> {:?}[:{}] {:?}[:{}]", lh, nh + 1, le, ne);
                    pf_cache.numerator_size = 0;
                    pf_expr.element_partial_fractioning(&lh[0..=nh], &le[0..ne], &mut pf_cache);
                }
            }
            return pf_expr;
        }
    }

    /// Add a new element PartialFractioningMonomial to the partial fractioning container
    pub fn add(&mut self, ellipsoids: &[(usize, usize)], num: &[usize]) {
        self.partial_fractioning_element
            .push(PartialFractioningMonomial::new(ellipsoids, num));
    }

    #[allow(dead_code)]
    /// Perform the partial fractioning on an element with
    ///  - cut        :  h_index[0]
    ///  - hyperboloid: h_index[1..]
    ///  - ellipsoids : e_index
    pub fn element_partial_fractioning(
        &mut self,
        h_index: &[usize],
        e_index: &[usize],
        pf_cache: &mut PFCache,
    ) -> bool {
        // Get nh and ne from the input so that we can import all the indices in
        // a single slice
        let n_e = e_index.len();
        // If all the surfaces are hyperboloids we get a pure numerator contribution
        if n_e == 0 {
            if self.numerator_rank + 1 < h_index.len() {
                return false;
            }
            for (&j, num) in h_index
                .iter()
                .zip_eq(pf_cache.numerator[0..self.n_props_deg].iter_mut())
            {
                *num = j;
            }
            self.add(
                &pf_cache.ellipsoids_product[0..0],
                &pf_cache.numerator[0..self.n_props_deg],
            );
            return true;
        }
        let n_h = h_index.len() - 1;
        let k = h_index[0];

        // Update the numerator structure
        pf_cache.numerator[pf_cache.numerator_size] = k;
        pf_cache.numerator_size += 1;

        //If all the surfaces are ellipsoids we already are in the final form
        if n_h == 0 {
            for (&j, prod) in e_index.iter().zip_eq(
                pf_cache.ellipsoids_product[0..self.n_props_deg - pf_cache.numerator_size]
                    .iter_mut(),
            ) {
                *prod = (k, j);
            }
            self.add(
                &pf_cache.ellipsoids_product[0..self.n_props_deg - pf_cache.numerator_size],
                &pf_cache.numerator[0..pf_cache.numerator_size],
            );
            return true;
        }

        // generate all the pure ellipsoids expression coming from this combination of
        // hyperboloids and ellipsoids surfaces
        let mut first_split = true;
        for (n, v) in pf_cache.splits.iter_mut().enumerate() {
            *v = n;
        }
        while first_split || get_next_subset(&mut pf_cache.splits[0..n_h], n_e + n_h - 1, true) {
            if first_split {
                first_split = false
            };
            for i in 0..=pf_cache.splits[0] {
                pf_cache.ellipsoids_product[i] = (k, e_index[i]);
            }
            for n in 0..n_h {
                if n + 1 != n_h {
                    for i in pf_cache.splits[n]..pf_cache.splits[n + 1] {
                        pf_cache.ellipsoids_product[i + 1] = (h_index[n + 1], e_index[i - n]);
                    }
                } else {
                    for i in pf_cache.splits[n]..n_e + n_h - 1 {
                        pf_cache.ellipsoids_product[i + 1] = (h_index[n + 1], e_index[i - n]);
                    }
                }
            }
            self.add(
                &pf_cache.ellipsoids_product[0..self.n_props_deg - pf_cache.numerator_size],
                &pf_cache.numerator[0..pf_cache.numerator_size],
            );
        }

        // Recoursively reduce a lower rank topology due to the numberator
        // Because at each iteration the rank of the numerator is lowerd by one
        // we can truncate the recursion earlier in many cases
        if self.numerator_rank + 1 == pf_cache.numerator_size {
            return false;
        } else {
            return self.element_partial_fractioning(&h_index[1..=n_h], e_index, pf_cache);
        }
    }

    /// Evaluate the partial fractioned expression using the information contained in LTDCache
    ///  - propagator_powers
    ///  - reduced_coefficient_lb
    ///  - complex_ellipsoids
    pub fn evaluate<T: FloatLike>(
        &self,
        ltd_numerator: &LTDNumerator,
        ll: &[LoopLine],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        // make sure that numerator.evaluate_reduced_in_lb has been called before this function
        // is not done here to avoid multiple calls in the case of amplitudes
        let mut result: na::Complex<T> = Complex::default();
        // Compute the overall factor coming from partial fractioning
        let mut norm: na::Complex<T> = Complex::new(-T::one(), T::zero());
        let mut min_index = self.n_props_deg;
        for p in ll[0].propagators.iter() {
            // Build the map to read the indices in PartialFractioningMonomial
            // in terms of propagator's id
            for _ in 0..cache.propagator_powers[p.id] {
                min_index -= 1;
                cache.pf_cache.numerator_index_map[min_index] = p.id;
            }
            norm *= (cache.complex_cut_energies[p.id] * Into::<T>::into(-2.0))
                .powi(cache.propagator_powers[p.id] as i32);
        }
        let mut skip: bool;
        let mut ellipsoid_count;
        for mono in self.partial_fractioning_element.iter() {
            // check if the all the ellipsoids exists otherwise move to the next entry
            skip = false;
            ellipsoid_count = 0;
            for (i1, i2) in mono.ellipsoids_product.iter() {
                if *i1 < min_index || *i2 < min_index {
                    skip = true;
                    break;
                }
                cache.pf_cache.ellipsoids_product[ellipsoid_count] = (
                    cache.pf_cache.numerator_index_map[*i1],
                    cache.pf_cache.numerator_index_map[*i2],
                );
                ellipsoid_count += 1;
            }
            if skip {
                continue;
            }

            cache.pf_cache.numerator_size = 0;
            for i in mono.numerator.iter() {
                if *i >= min_index
                    && cache.propagator_powers[cache.pf_cache.numerator_index_map[*i]] != 0
                {
                    cache.pf_cache.numerator[cache.pf_cache.numerator_size] =
                        cache.pf_cache.numerator_index_map[*i];
                    cache.pf_cache.numerator_size += 1;
                }
            }
            //            println!(
            //                "USE:: {:?} \n\t-> {:?} num{:?}",
            //                mono,
            //                &pf_cache.ellipsoids_product[0..ellipsoid_count],
            //                &pf_cache.numerator[0..pf_cache.numerator_size]
            //            );
            result += PartialFractioningMonomial::evaluate_ellipsoids_product(
                &cache.pf_cache.ellipsoids_product[0..ellipsoid_count],
                cache,
            )
            .inv()
                * PartialFractioningMonomial::evaluate_numerator(
                    &cache.pf_cache.numerator[0..cache.pf_cache.numerator_size],
                    ltd_numerator,
                    ll,
                    cache,
                );
        }
        return result / norm;
    }
}
/* MULTI LOOPS  */

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningDen {
    pub lambdas: Vec<f64>,
    pub energies: Vec<f64>,
    pub shifts: Vec<f64>,
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningBlock {
    pub factor: f64,
    pub denominators: Vec<PartialFractioningDen>,
    // TODO: numerator
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningMultiLoops {
    pub partial_fractioning_element: Vec<PartialFractioningBlock>,
    pub loop_lines: Vec<LoopLine>,
    pub ll_n_props_deg: Vec<usize>,
    pub n_loops: usize,
    //splits: Vec<usize>,
}
impl PartialFractioningBlock {
    pub fn evaluate<T: FloatLike>(
        &self,
        numerator: &LTDNumerator,
        loop_lines: &[LoopLine],
        min_index: usize,
        map_id: &[(usize, usize)],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        println!("===============================================");
        if self.denominators.len() == 0 {
            return Complex::default();
        }
        let den = self.evaluate_denominators(loop_lines, min_index, map_id, ltd_cache);
        println!("den: {}", den);
        //let num = PartialFractioningBlock::evaluate_numerator()

        return den.inv() * Into::<T>::into(self.factor); // * num;
    }

    pub fn evaluate_denominators<T: FloatLike>(
        &self,
        loop_lines: &[LoopLine],
        min_index: usize,
        map_id: &[(usize, usize)],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let mut res = Complex::new(T::one(), T::zero());
        let mut den_res: Complex<T> = Complex::default();
        for (n, den) in self.denominators.iter().enumerate() {
            println!("  :: sub [{}]", n);
            den_res = Complex::default();
            for (idx, (v_e, v_s)) in den
                .energies
                .iter()
                .zip(den.shifts.iter())
                .enumerate()
                .skip(min_index)
            {
                if *v_e < 1e-10 {
                    continue;
                }
                let id = ltd_cache.pf_cache.numerator_index_map[idx];
                println!("\t| {}, {}, {}", id, v_e, v_s);
                if ltd_cache.propagator_powers[id] == 0 {
                    panic!("Ellipsoid do not exists!");
                }
                den_res += ltd_cache.complex_cut_energies[id] * Into::<T>::into(*v_e);
                den_res +=
                    Into::<T>::into(loop_lines[map_id[id].0].propagators[map_id[id].1].q.t * v_s);
            }
            println!("\t> den_res = {}", den_res);

            res *= den_res;
        }
        return res;
    }
}

impl PartialFractioningMultiLoops {
    pub fn new(loop_lines: &Vec<LoopLine>, n_loops: usize) -> PartialFractioningMultiLoops {
        let mut pf_expr = PartialFractioningMultiLoops {
            partial_fractioning_element: Vec::new(),
            loop_lines: loop_lines.clone(),
            ll_n_props_deg: loop_lines
                .iter()
                .map(|x| x.propagators.iter().map(|x| x.power).sum())
                .collect(),
            n_loops: n_loops,
        };
        pf_expr.partial_fractioning_nl();
        return pf_expr;
    }

    fn partial_fractioning_nl(&mut self) {
        // ll_n_props_deg: number of propagator per loop line
        // signatures: signature of each loop line
        // sigmas: sign coming from closing the contour

        //id_to_ll = []
        let mut id_to_ll = Vec::new();
        for (n, &n_props_deg) in self.ll_n_props_deg.iter().enumerate() {
            for _ in 0..n_props_deg {
                id_to_ll.push(n)
            }
        }
        println!("id_to_ll: {:?}", id_to_ll);
        let n_props_deg = self.ll_n_props_deg.iter().sum();

        // Take all combination of plus and minus energies from the first E_fractioning
        //    0: positive energy
        //    1: negative energy
        for choose in (0..n_props_deg).map(|_| 0..=1).multi_cartesian_product() {
            //println!("choose: {:?}", choose);
            let mut product = PartialFractioningBlock {
                factor: 1.0,
                denominators: Vec::new(),
            };
            for (n, &which) in choose.iter().enumerate() {
                // the positive enegy component comes with a - sign
                product.factor *= -(1 - 2 * which) as f64;

                let mut den = PartialFractioningDen {
                    lambdas: self.loop_lines[id_to_ll[n]]
                        .signature
                        .iter()
                        .map(|x| *x as f64)
                        .collect(),
                    energies: vec![0.0; n_props_deg],
                    shifts: vec![0.0; n_props_deg],
                };
                den.energies[n] = (1 - 2 * which) as f64;
                den.shifts[n] = 1.0;
                product.denominators.push(den);
            }

            self.pf_product(&mut product, 0);
            //        for n, res in enumerate(pf_product(product, n_loops)):
            //            print("PF RESULT :: %d" %n)
            //            for d in res[1]:
            //                    print("\t", d)
        }
        //return pf_res;
    }

    /// append new result
    fn add(&mut self, product: &mut PartialFractioningBlock) {
        // TODO: Check that all the lambdas are now zero
        self.partial_fractioning_element
            .push(PartialFractioningBlock {
                factor: product.factor,
                denominators: product.denominators.iter().map(|x| x.clone()).collect(),
            })
    }

    // Perform partial fractioning on a single product for an arbitrary number of loops
    fn pf_product(&mut self, product: &mut PartialFractioningBlock, residue_n: usize) -> bool {
        if residue_n == self.n_loops {
            self.add(product);
            return true;
            //return [(global_factor, [{k: v for k, v in x.items() if k != 'lambdas'} for x in product], numerator)]
        }
        //println!("========================================");
        //println!("             LOOP {}", residue_n);
        //println!("========================================");
        let mut indices = Vec::new();
        let mut h_index = Vec::new();
        let mut e_index = Vec::new();
        let mut spectators = Vec::new();

        for (n, den) in product.denominators.iter_mut().enumerate() {
            let norm = den.lambdas[residue_n];
            if norm.abs() > 1e-10 {
                // Element depneds on the loop momenta

                // Extract the factor in front of the loop momenta in order to
                // to have a consistent unit factor before taking the residue
                product.factor /= norm;
                for (en, sn) in den.energies.iter_mut().zip(den.shifts.iter_mut()) {
                    *en /= norm;
                    *sn /= norm;
                }
                for ln in den.lambdas.iter_mut() {
                    *ln /= norm;
                }

                indices.push(n);
                if den.energies.iter().all(|x| *x >= 0.0) {
                    e_index.push(n);
                } else {
                    h_index.push(n);
                }
            } else {
                spectators.push(n)
                //left_product += [den]
            }
        }
        //println!("  idx: {:?}, hdx: {:?}, edx: {:?}", indices, h_index, e_index);
        // Correction coming from applying partial fractioning to remove hyperboloids
        product.factor *= (-1.0 as f64).powi(1 + h_index.len() as i32);
        if h_index.len() == 0 {
            // no pole
            return true;
        } else {
            // apply residue and partial fractioning

            // TODO: use some sort of cache
            // Create container
            let mut pf_expr = PartialFractioning {
                partial_fractioning_element: Vec::new(),
                n_props_deg: h_index.len() + e_index.len(),
                numerator_rank: 0,
            };
            // Use this to avoid redundant allocation
            let mut pf_cache = PFCache::new(h_index.len() + e_index.len());
            pf_cache.numerator_size = 0;
            pf_expr.element_partial_fractioning(&h_index, &e_index, &mut pf_cache);

            if pf_expr.partial_fractioning_element.len() == 0 {
                return true;
            }

            for mapping in pf_expr.partial_fractioning_element.iter() {
                let mut new_product = PartialFractioningBlock {
                    factor: product.factor,
                    denominators: Vec::new(),
                };
                for i in spectators.iter() {
                    new_product
                        .denominators
                        .push(product.denominators[*i].clone());
                }
                //print("\t   > ", mapping)
                // Add partial fractioned elements
                for pair in mapping.ellipsoids_product.iter() {
                    PartialFractioningMultiLoops::pf_mapper(
                        &mut new_product,
                        &product.denominators[pair.0],
                        &product.denominators[pair.1],
                        residue_n,
                    );
                }

                // Take next residue
                self.pf_product(&mut new_product, residue_n + 1);
            }
        }

        return true;
    }

    /// takes two denominator, one is the one over which the residue is taken (giver)
    /// the other one is one where we have to replace the value of the loop momentum (receiver)
    /// The result is appended to the PartialFractioningBlock given as first argument
    fn pf_mapper(
        product: &mut PartialFractioningBlock,
        den_giver: &PartialFractioningDen,
        den_receiver: &PartialFractioningDen,
        residue_n: usize,
    ) {
        let factor = den_receiver.lambdas[residue_n] / den_giver.lambdas[residue_n];
        product.denominators.push(PartialFractioningDen {
            lambdas: den_receiver
                .lambdas
                .iter()
                .zip(den_giver.lambdas.iter())
                .map(|(x1, x2)| x1 - factor * x2)
                .collect(),
            energies: den_receiver
                .energies
                .iter()
                .zip(den_giver.energies.iter())
                .map(|(x1, x2)| x1 - factor * x2)
                .collect(),
            shifts: den_receiver
                .shifts
                .iter()
                .zip(den_giver.shifts.iter())
                .map(|(x1, x2)| x1 - factor * x2)
                .collect(),
        });
    }
    /// Evaluate the partial fractioned expression using the information contained in LTDCache
    ///  - propagator_powers
    ///  - reduced_coefficinet_lb
    ///  - complex_ellipsoids
    pub fn evaluate<T: FloatLike>(
        &self,
        ltd_numerator: &LTDNumerator,
        loop_lines: &[LoopLine],
        map_id: &[(usize, usize)],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        // WARNING: make sure that numerator.evaluate_reduced_in_lb has been called before this function
        // is not done here to avoid multiple calls in the case of amplitudes
        let mut result: na::Complex<T> = Complex::default();
        // Compute the overall factor coming from partial fractioning
        let mut norm: na::Complex<T> = Complex::new(T::one(), T::zero());
        let n_props_deg = self.ll_n_props_deg.iter().sum();
        dbg!(&n_props_deg);
        let mut min_index = n_props_deg;
        for ll in loop_lines.iter().rev() {
            //for ll in loop_lines.iter() {
            for p in ll.propagators.iter() {
                // Build the map to read the indices in PartialFractioningMonomial
                // in terms of propagator's id
                for _ in 0..cache.propagator_powers[p.id] {
                    min_index -= 1;
                    cache.pf_cache.numerator_index_map[min_index] = p.id;
                }
                norm *= (cache.complex_cut_energies[p.id] * Into::<T>::into(2.0))
                    .powi(cache.propagator_powers[p.id] as i32);
            }
        }
        dbg!(&min_index);
        let mut store: bool;
        let mut skip: bool;
        let mut new_block = self.partial_fractioning_element[0].clone();
        for block in self.partial_fractioning_element.iter() {
            // check if the all the ellipsoids exists otherwise move to the next entry
            skip = false;

            new_block.factor = block.factor;
            new_block.denominators.clear();
            for den in block.denominators.iter() {
                store = true;
                for (v_e, v_s) in den.energies[..min_index]
                    .iter()
                    .zip(den.shifts[..min_index].iter())
                {
                    // TODO: decide if we need to include abs()
                    // If we have fewer propagators then the one computed before we
                    // we can still recycle the original partial fractioning by discarding
                    // all the superfluous denominators
                    if *v_e > 0.0 && *v_s < 0.0 {
                        skip = true;
                        break;
                    } else if *v_e > 0.0 && *v_s > 0.0 {
                        store = store && false
                    } else {
                        store = store && true
                    }
                }
                if skip {
                    break;
                }

                if store {
                    new_block.denominators.push(den.clone());
                }
            }
            // COMING-SOON: numerator
            println! {"---"};
            println! {"ll_n_props_deg: {:?}", self.ll_n_props_deg};
            println! {"map: {:?}", cache.pf_cache.numerator_index_map};
            for den in new_block.denominators.iter() {
                println! {"{:?}", den.energies};
                println! {"{:?}", den.shifts};
            }
            if skip || new_block.denominators.len() + min_index != n_props_deg - self.n_loops {
                continue;
            } else {
                result += new_block.evaluate(ltd_numerator, loop_lines, min_index, map_id, cache);
            }
        }
        return result / norm;
    }
}
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use topologies::Propagators;
    use vector::LorentzVector;

    #[allow(non_snake_case, dead_code)]
    #[test]
    fn test_partial_fractioning() {
        let start = std::time::Instant::now();
        let pf_expr = PartialFractioning::new(10, 9);

        println!("{:?}", start.elapsed());
        println!(
            "{} -> size: {}",
            10,
            pf_expr.partial_fractioning_element.len()
        );
        assert_eq!(pf_expr.partial_fractioning_element.len(), 92378);
        //for expr in pf_expr.partial_fractioning_element.iter() {
        //    println!("{:?} :: num({:?})", expr.ellipsoids_product, expr.numerator);
        //}
    }
    #[test]
    fn test_partial_fractioning_nl() {
        let propagator = Propagators {
            id: 0, // global id
            m_squared: 0.0,
            q: LorentzVector::default(),
            parametric_shift: (Vec::new(), Vec::new()),
            signature: Vec::new(),
            power: 1,
        };
        // Define looplines
        // Replicate one loop
        let mut looplines_1loop = Vec::new();
        looplines_1loop.push(LoopLine {
            start_node: 0,
            end_node: 0,
            signature: vec![1],
            propagators: vec![propagator.clone(); 10],
        });
        // Replicate two loops
        let signatures = [vec![1, 0], vec![1, -1], vec![0, 1]];
        let multiplicity = [6, 1, 4];
        let mut looplines_2loop = Vec::new();
        for (s, m) in signatures.iter().zip(multiplicity.iter()) {
            looplines_2loop.push(LoopLine {
                start_node: 0,
                end_node: 0,
                signature: s.clone(),
                propagators: vec![propagator.clone(); m.clone()],
            });
        }
        // Replicate 2x2 FISHNET
        let signatures = [
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![1, -1, 0, 0],
            vec![-1, 0, -1, 0],
            vec![0, -1, 0, -1],
            vec![0, 0, 1, 0],
            vec![0, 0, -1, 1],
            vec![0, 0, 0, -1],
        ];
        let multiplicity = [1, 1, 1, 1, 1, 1, 1, 1];
        let mut looplines_2x2_fishnet = Vec::new();
        for (s, m) in signatures.iter().zip(multiplicity.iter()) {
            looplines_2x2_fishnet.push(LoopLine {
                start_node: 0,
                end_node: 0,
                signature: s.clone(),
                propagators: vec![propagator.clone(); m.clone()],
            });
        }

        // TEST 1 LOOP
        let start = std::time::Instant::now();
        let mut pf_expr = PartialFractioningMultiLoops::new(&looplines_1loop, 1);
        pf_expr.partial_fractioning_nl();
        print!("1LOOP -> #{}", pf_expr.partial_fractioning_element.len());
        println!(" :: {:?}", start.elapsed());
        assert_eq!(pf_expr.partial_fractioning_element.len(), 48620);

        // TEST 2 LOOPS
        let start = std::time::Instant::now();
        let mut pf_expr = PartialFractioningMultiLoops::new(&looplines_2loop, 2);
        pf_expr.partial_fractioning_nl();
        print!("2LOOPS -> #{}", pf_expr.partial_fractioning_element.len());
        println!(" :: {:?}", start.elapsed());
        assert_eq!(pf_expr.partial_fractioning_element.len(), 48620);
        // TEST 2x2 FISHNET
        let start = std::time::Instant::now();
        let mut pf_expr = PartialFractioningMultiLoops::new(&looplines_2x2_fishnet, 4);
        pf_expr.partial_fractioning_nl();
        print!("FISHNET - >#{}", pf_expr.partial_fractioning_element.len());
        println!(" :: {:?}", start.elapsed());
        assert_eq!(pf_expr.partial_fractioning_element.len(), 98);
    }
}
