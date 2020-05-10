use itertools::Itertools;
use num::Complex;
use topologies::{LTDCache, LTDNumerator, LoopLine, Propagators, Topology};
use vector::LorentzVector;
use FloatLike;
use {MAX_LOOP, MAX_PROP};

#[derive(Default, Debug, Clone)]
pub struct PartialFractioning {
    pub partial_fractioning_element: Vec<PartialFractioningMonomial>,
    pub n_props_deg: usize,
    pub numerator_rank: usize,
    //splits: Vec<usize>,
}

#[derive(Default, Debug, Clone)]
pub struct PFCache<T: FloatLike> {
    splits: Vec<usize>,
    ellipsoids_product: Vec<(usize, usize)>,
    den_short: PartialFractioningDen<T>,
    block_fiter: Vec<bool>,
    numerator: Vec<usize>,
    numerator_size: usize,
    numerator_index_map: Vec<usize>,
}

impl<T: FloatLike> PFCache<T> {
    pub fn new(n_props_deg: usize) -> PFCache<T> {
        if n_props_deg == 0 {
            PFCache {
                splits: vec![],
                ellipsoids_product: vec![],
                den_short: PartialFractioningDen {
                    indices: vec![],
                    energies_and_shifts: vec![],
                    size: 0,
                },
                block_fiter: vec![],
                numerator: vec![],
                numerator_size: 0,
                numerator_index_map: vec![],
            }
        } else {
            // Use this to avoid redundant allocation
            PFCache {
                splits: vec![0; n_props_deg - 1],
                ellipsoids_product: vec![(0, 0); n_props_deg - 1],
                den_short: PartialFractioningDen {
                    indices: vec![0; 2 * MAX_LOOP],
                    energies_and_shifts: vec![(0.0, 0.0); 2 * MAX_LOOP],
                    size: 0,
                },
                block_fiter: vec![false; MAX_PROP],
                numerator: vec![0; n_props_deg],
                numerator_size: 0,
                numerator_index_map: vec![0; n_props_deg],
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
        numerator: &LTDNumerator<T>,
        loop_line: &[LoopLine<T>], // make sure to only pass loop lines with a loop momentum in it
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
        ltd_numerator: &LTDNumerator<T>,
        ll: &[LoopLine<T>],
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
        ltd_numerator: &LTDNumerator<T>,
        ll: &[LoopLine<T>],
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
    pub fn element_partial_fractioning<T: FloatLike>(
        &mut self,
        h_index: &[usize],
        e_index: &[usize],
        pf_cache: &mut PFCache<T>,
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
        for (n, v) in pf_cache.splits[0..self.n_props_deg - 1]
            .iter_mut()
            .enumerate()
        {
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
        ltd_numerator: &LTDNumerator<T>,
        ll: &[LoopLine<T>],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        println!("1L map: {:?}", cache.pf_cache.numerator_index_map);
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
        dbg!(&norm);
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
            if mono.numerator.len() == 1 {
                println!(
                    //"USE:: {:?} num {:?} -> {:?} num {:?}",
                    "USE:: {:?} -> {:?}",
                    mono.ellipsoids_product,
                    //mono.numerator,
                    &cache.pf_cache.ellipsoids_product[0..ellipsoid_count],
                    //&cache.pf_cache.numerator[0..cache.pf_cache.numerator_size]
                );
                println!(
                    //"\t den: {:?}, num: {:?}",
                    "\t den: {:?}",
                    PartialFractioningMonomial::evaluate_ellipsoids_product(
                        &cache.pf_cache.ellipsoids_product[0..ellipsoid_count],
                        cache,
                    )
                    .inv(),
                    //                    PartialFractioningMonomial::evaluate_numerator(
                    //                        &cache.pf_cache.numerator[0..cache.pf_cache.numerator_size],
                    //                        ltd_numerator,
                    //                        ll,
                    //                        cache,
                    //                    )
                );
            }
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
pub struct PartialFractioningDen<T: FloatLike> {
    pub indices: Vec<u8>,
    pub energies_and_shifts: Vec<(T, T)>,
    pub size: usize,
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningBlock<T: FloatLike> {
    pub factor: T,
    pub denominators: Vec<PartialFractioningDen<T>>,
    // TODO: numerator
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningDen2<T: FloatLike> {
    pub lambdas: Vec<T>,
    pub energies: Vec<T>,
    pub shifts: Vec<T>,
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningBlock2<T: FloatLike> {
    pub factor: T,
    pub denominators: Vec<PartialFractioningDen2<T>>,
    // TODO: numerator
}

#[derive(Default, Debug, Clone)]
pub struct PartialFractioningMultiLoops<T: FloatLike> {
    pub partial_fractioning_element: Vec<PartialFractioningBlock<T>>,
    pub loop_lines: Vec<LoopLine<T>>,
    pub ll_n_props_deg: Vec<usize>,
    pub n_loops: usize,
    //splits: Vec<usize>,
}
impl<T: FloatLike> PartialFractioningDen<T> {
    pub fn evaluate(
        &self,
        loop_lines: &[LoopLine<T>],
        min_index: usize,
        map_id: &[(usize, usize)],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let mut den_res = Complex::default();
        for (idx, (v_e, v_s)) in self.indices[..self.size]
            .iter()
            .zip_eq(self.energies_and_shifts[..self.size].iter())
            .skip(min_index)
        {
            let id = ltd_cache.pf_cache.numerator_index_map[*idx as usize];
            if ltd_cache.propagator_powers[id] == 0 {
                panic!("Ellipsoid do not exists!");
            }
            den_res += ltd_cache.complex_cut_energies[id] * Into::<T>::into(*v_e);
            den_res +=
                Into::<T>::into(loop_lines[map_id[id].0].propagators[map_id[id].1].q.t * v_s);
        }
        den_res.inv()
    }
}
impl<T: FloatLike> PartialFractioningBlock<T> {
    pub fn evaluate(
        &self,
        numerator: &LTDNumerator<T>,
        loop_lines: &[LoopLine<T>],
        min_index: usize,
        map_id: &[(usize, usize)],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        //println!("===============================================");
        if self.denominators.len() == 0 {
            return Complex::default();
        }
        let den_inv = self.evaluate_dens(loop_lines, min_index, map_id, ltd_cache);
        //println!("den: {}", den);
        //let num = PartialFractioningBlock::evaluate_numerator()

        return den_inv * Into::<T>::into(self.factor); // * num;
    }

    pub fn evaluate_dens(
        &self,
        loop_lines: &[LoopLine<T>],
        min_index: usize,
        map_id: &[(usize, usize)],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let mut res = Complex::new(T::one(), T::zero());
        for den in self.denominators.iter() {
            //println!("  :: sub [{}]", n);
            res *= den.evaluate(loop_lines, min_index, map_id, ltd_cache);
        }
        return res;
    }
}

impl<T: FloatLike> PartialFractioningMultiLoops<T> {
    pub fn new(loop_lines: &Vec<LoopLine<T>>) -> PartialFractioningMultiLoops<T> {
        let mut pf_expr = PartialFractioningMultiLoops {
            partial_fractioning_element: Vec::new(),
            loop_lines: loop_lines.clone(),
            ll_n_props_deg: loop_lines
                .iter()
                .map(|x| x.propagators.iter().map(|x| x.power).sum())
                .collect(),
            n_loops: loop_lines[0].signature.len(),
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
        //println!("id_to_ll: {:?}", id_to_ll);
        let n_props_deg = self.ll_n_props_deg.iter().sum();
        // create cache
        let mut pf_cache = PFCache::new(n_props_deg);

        // Take all combination of plus and minus energies from the first E_fractioning
        //    1.0: positive energy
        //   -1.0: negative energy
        for choose in (0..n_props_deg)
            .map(|_| [-1.0, 1.0].iter())
            .multi_cartesian_product()
        {
            //println!("choose: {:?}", choose);
            let mut product = PartialFractioningBlock2 {
                factor: 1.0,
                denominators: Vec::new(),
            };
            for (n, &which) in choose.iter().enumerate() {
                // the positive enegy component comes with a - sign
                product.factor *= -which;

                let mut den = PartialFractioningDen2 {
                    lambdas: self.loop_lines[id_to_ll[n]]
                        .signature
                        .iter()
                        .map(|x| *x as T)
                        .collect(),
                    energies: vec![0.0; n_props_deg],
                    shifts: vec![0.0; n_props_deg],
                };
                den.energies[n] = *which;
                den.shifts[n] = 1.0;
                product.denominators.push(den);
            }

            self.pf_product(&mut product, 0, &mut pf_cache);
        }
    }

    /// append new result
    fn add(&mut self, product: &mut PartialFractioningBlock2<T>, pf_cache: &mut PFCache<T>) {
        // TODO: Check that all the lambdas are now zero
        let mut dens_short = Vec::new();
        for den in product.denominators.iter() {
            //            println!("energies: {:?}", den.energies);
            //            println!("shifts  : {:?}", den.shifts);

            // Reset the vectors before filling them again
            pf_cache.den_short.size = 0;
            for (idx, (v_e, v_s)) in den.energies.iter().zip(den.shifts.iter()).enumerate() {
                if *v_e == 0.0 {
                    assert_eq!(v_s, v_e);
                    continue;
                }
                pf_cache.den_short.indices[pf_cache.den_short.size] = idx as u8;
                pf_cache.den_short.energies_and_shifts[pf_cache.den_short.size] = (*v_e, *v_s);
                pf_cache.den_short.size += 1;
            }
            //println!("indices   : {:?}", pf_cache.den_short.indices);
            //println!("en_and_sh : {:?}", pf_cache.den_short.energies_and_shifts);
            // Store the value in dens_short
            dens_short.push(pf_cache.den_short.clone());
        }
        self.partial_fractioning_element
            .push(PartialFractioningBlock {
                factor: product.factor,
                denominators: dens_short,
            });
    }

    // Perform partial fractioning on a single product for an arbitrary number of loops
    fn pf_product(
        &mut self,
        product: &mut PartialFractioningBlock2<T>,
        residue_n: usize,
        pf_cache: &mut PFCache<T>,
    ) -> bool {
        if residue_n == self.n_loops {
            self.add(product, pf_cache);
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
            if norm.abs() > 0.0 {
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
        product.factor *= (-1.0 as T).powi(1 + h_index.len() as i32);
        if h_index.len() == 0 {
            // no pole
            return true;
        }
        // apply residue and partial fractioning

        // TODO: use some sort of cache
        // Create container
        let mut pf_expr = PartialFractioning {
            partial_fractioning_element: Vec::new(),
            n_props_deg: h_index.len() + e_index.len(),
            numerator_rank: 0,
        };
        // Use this to avoid redundant allocation
        //       let mut pf_cache2 = PFCache::new(h_index.len() + e_index.len());
        pf_cache.numerator_size = 0;

        pf_expr.element_partial_fractioning(&h_index, &e_index, pf_cache);

        if pf_expr.partial_fractioning_element.len() == 0 {
            return true;
        }

        for mapping in pf_expr.partial_fractioning_element.iter() {
            let mut new_product = PartialFractioningBlock2 {
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
            self.pf_product(&mut new_product, residue_n + 1, pf_cache);

            //TODO: RESET HERE
        }
        true
    }

    /// takes two denominator, one is the one over which the residue is taken (giver)
    /// the other one is one where we have to replace the value of the loop momentum (receiver)
    /// The result is appended to the PartialFractioningBlock given as first argument
    fn pf_mapper(
        product: &mut PartialFractioningBlock2<T>,
        den_giver: &PartialFractioningDen2<T>,
        den_receiver: &PartialFractioningDen2<T>,
        residue_n: usize,
    ) {
        let factor = den_receiver.lambdas[residue_n] / den_giver.lambdas[residue_n];
        product.denominators.push(PartialFractioningDen2 {
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
    pub fn evaluate(
        &self,
        ltd_numerator: &LTDNumerator<T>,
        loop_lines: &[LoopLine<T>],
        map_id: &[(usize, usize)],
        cache: &mut LTDCache<T>,
    ) -> Complex<T> {
        println!("NL map: {:?}", cache.pf_cache.numerator_index_map);
        // WARNING: make sure that numerator.evaluate_reduced_in_lb has been called before this function
        // is not done here to avoid multiple calls in the case of amplitudes
        let mut result: na::Complex<T> = Complex::default();
        // Compute the overall factor coming from partial fractioning
        let mut norm: na::Complex<T> = Complex::new(T::one(), T::zero());
        let n_props_deg = self.ll_n_props_deg.iter().sum();
        dbg!(&n_props_deg, &self.ll_n_props_deg);
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
        dbg!(&min_index, &norm);
        let mut skip: bool;

        for block in self.partial_fractioning_element.iter() {
            // check if the all the ellipsoids exists otherwise move to the next entry
            let mut block_size = 0;
            skip = false;
            let mut block_res = Complex::new(Into::<T>::into(block.factor), T::zero());

            //new_block.dens.clear();
            for (den_n, den) in block.denominators.iter().enumerate() {
                // clear the short denominator
                cache.pf_cache.den_short.size = 0;
                // Check if we need to evaluate it
                cache.pf_cache.block_fiter[den_n] = true;
                for (idx, (v_e, v_s)) in den.indices[..den.size]
                    .iter()
                    .zip(den.energies_and_shifts[..den.size].iter())
                {
                    if (*idx as usize) >= min_index {
                        continue;
                    }
                    // TODO: decide if we need to include abs()
                    // If we have fewer propagators then the one computed before we
                    // we can still recycle the original partial fractioning by discarding
                    // all the superfluous denominators
                    if *v_e > 0.0 && *v_s < 0.0 {
                        skip = true;
                        break;
                    } else if *v_e > 0.0 && *v_s > 0.0 {
                        cache.pf_cache.block_fiter[den_n] = false;
                        break;
                    }
                }
                if skip {
                    break;
                }
                if cache.pf_cache.block_fiter[den_n] {
                    block_size += 1;
                }
            }
            // COMING-SOON: numerator
            //for (store, den) in cache.pf_cache.block_fiter[..n_props]
            //    .iter()
            //    .zip(block.dens.iter())
            //{
            //    if *store {
            //        println!("{:?}", den.indices);
            //        println!("{:?}", den.energies_and_shifts);
            //        println!("{:?}\n", den.size);
            //    }
            //}

            if skip || block_size + min_index != n_props_deg - self.n_loops {
                //println!("SKIP");
                continue;
            } else {
                for (store, den) in cache.pf_cache.block_fiter[..n_props_deg]
                    .iter()
                    .zip(block.denominators.iter())
                {
                    if *store {
                        block_res *= den.evaluate(loop_lines, min_index, map_id, cache);
                    }
                }
                //                print!("USE:: [");
                //                for den in block.dens.iter() {
                //                    let d1 = den.energies_and_shifts[0];
                //                    let e1 = den.indices[0];
                //                    let e2 = den.indices[1];
                //                    if d1.1 < 0.0 {
                //                        print!("({},{}), ", e1, e2);
                //                    } else {
                //                        print!("({},{}), ", e2, e1);
                //                    }
                //                }
                //                print!("] -> [");
                //                for den in block.dens.iter() {
                //                    let d1 = den.energies_and_shifts[0];
                //                    let e1 = cache.pf_cache.numerator_index_map[den.indices[0] as usize];
                //                    let e2 = cache.pf_cache.numerator_index_map[den.indices[1] as usize];
                //                    if d1.1 < 0.0 {
                //                        print!("({},{}), ", e1, e2);
                //                    } else {
                //                        print!("({},{}), ", e2, e1);
                //                    }
                //                }
                //                println!("]");
                //
                //                println!("\t den: {:?}", block_res);
                result += block_res;
            }
        }
        return result / norm;
    }
}
