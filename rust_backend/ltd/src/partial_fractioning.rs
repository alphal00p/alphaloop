use itertools::Itertools;
use num::Complex;
use topologies::{LTDCache, LTDNumerator, LoopLine, Topology};
use FloatLike;

#[derive(Default, Debug, Clone)]
pub struct PartialFractioning {
    pub partial_fractioning_element: Vec<PartialFractioningMonomial>,
    pub n_props: usize,
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
        topo: &Topology,
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        let den = PartialFractioningMonomial::evaluate_ellipsoids_product(
            &self.ellipsoids_product,
            ltd_cache,
        );
        let num = PartialFractioningMonomial::evaluate_numerator(
            self.numerator.as_slice(),
            numerator,
            &topo.loop_lines,
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
        if ltd_numerator.max_rank < pf_numerator.len() - 1 {
            return Complex::default();
        }
        return PartialFractioningMonomial::evaluate_numerator2(
            pf_numerator,
            ltd_numerator.max_rank + 1 - pf_numerator.len(),
            0,
            ltd_numerator,
            ll,
            ltd_cache,
        );
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
        let mut x;
        if pf_numerator.len() == 1 {
            for n in offset..=rank {
                x = ltd_cache.complex_cut_energies[ll[0].propagators[pf_numerator[0]].id]
                    - Into::<T>::into(ll[0].propagators[pf_numerator[0]].q.t);
                res += ltd_cache.reduced_coefficient_lb
                    //[numerator.powers_to_position[&[n as u8, 0, 0, 0]]]
                    [ltd_numerator.max_rank - n]// Only at 1-Loop
                    * x.powi((n - offset) as i32);
            }
            return res;
        }
        // TODO: Multiplication by the ltd_cache.... is complete bollocks
        // Waiting to have the right thing
        for n in offset..=rank {
            x = ltd_cache.complex_cut_energies[ll[0].propagators[pf_numerator[0]].id]
                - Into::<T>::into(ll[0].propagators[pf_numerator[0]].q.t);
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
                n_props: n_prop,
                numerator_rank: numerator_rank,
            };

            // Use this to avoid redundant allocation
            let mut cache = PFCache::new(n_prop);

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
                cache.numerator_size = 0;
                pf_expr.element_partial_fractioning(&lh[0..=nh], &le[0..ne], &mut cache);
                while get_next_set_pair(&mut lh[0..=nh], &mut le[0..ne], n_prop) {
                    //println!(" -> {:?}[:{}] {:?}[:{}]", lh, nh + 1, le, ne);
                    cache.numerator_size = 0;
                    pf_expr.element_partial_fractioning(&lh[0..=nh], &le[0..ne], &mut cache);
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
    fn element_partial_fractioning(
        &mut self,
        h_index: &[usize],
        e_index: &[usize],
        cache: &mut PFCache,
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
                .zip_eq(cache.numerator[0..self.n_props].iter_mut())
            {
                *num = j;
            }
            self.add(
                &cache.ellipsoids_product[0..0],
                &cache.numerator[0..self.n_props],
            );
            return true;
        }
        let n_h = h_index.len() - 1;
        let k = h_index[0];

        // Update the numerator structure
        cache.numerator[cache.numerator_size] = k;
        cache.numerator_size += 1;

        //If all the surfaces are ellipsoids we already are in the final form
        if n_h == 0 {
            for (&j, prod) in e_index
                .iter()
                .zip_eq(cache.ellipsoids_product[0..self.n_props - cache.numerator_size].iter_mut())
            {
                *prod = (k, j);
            }
            self.add(
                &cache.ellipsoids_product[0..self.n_props - cache.numerator_size],
                &cache.numerator[0..cache.numerator_size],
            );
            return true;
        }

        // generate all the pure ellipsoids expression coming from this combination of
        // hyperboloids and ellipsoids surfaces
        let mut first_split = true;
        for (n, v) in cache.splits.iter_mut().enumerate() {
            *v = n;
        }
        while first_split || get_next_subset(&mut cache.splits[0..n_h], n_e + n_h - 1, true) {
            if first_split {
                first_split = false
            };
            for i in 0..=cache.splits[0] {
                cache.ellipsoids_product[i] = (k, e_index[i]);
            }
            for n in 0..n_h {
                if n + 1 != n_h {
                    for i in cache.splits[n]..cache.splits[n + 1] {
                        cache.ellipsoids_product[i + 1] = (h_index[n + 1], e_index[i - n]);
                    }
                } else {
                    for i in cache.splits[n]..n_e + n_h - 1 {
                        cache.ellipsoids_product[i + 1] = (h_index[n + 1], e_index[i - n]);
                    }
                }
            }
            self.add(
                &cache.ellipsoids_product[0..self.n_props - cache.numerator_size],
                &cache.numerator[0..cache.numerator_size],
            );
        }

        // Recoursively reduce a lower rank topology due to the numberator
        // Because at each iteration the rank of the numerator is lowerd by one
        // we can truncate the recursion earlier in many cases
        if self.numerator_rank + 1 == cache.numerator_size {
            return false;
        } else {
            return self.element_partial_fractioning(&h_index[1..=n_h], e_index, cache);
        }
    }

    /// Evaluate the partial fractioned expression using the information contained in LTDCache
    ///  - propagator_powers
    ///  - reduced_coefficinet_lb
    ///  - complex_ellipsoids
    pub fn evaluate<T: FloatLike>(
        &self,
        ltd_numerator: &LTDNumerator,
        ll: &[LoopLine],
        ltd_cache: &LTDCache<T>,
    ) -> Complex<T> {
        // make sure that numerator.evaluate_reduced_in_lb has been called before this function
        // is not done here to avoid multiple calls in the case of amplitudes
        let mut pf_cache = PFCache::new(self.n_props);
        let mut result: na::Complex<T> = Complex::default();
        // Compute the overall factor coming from partial fractioning
        let mut norm: na::Complex<T> = Complex::new(-T::one(), T::zero());
        let mut min_index = self.n_props;
        for p in ll[0].propagators.iter() {
            // Build the map to read the indices in PartialFractioningMonomial
            // in terms of propagator's id
            for _ in 0..ltd_cache.propagator_powers[p.id] {
                min_index -= 1;
                pf_cache.numerator_index_map[min_index] = p.id;
            }
            norm *= (ltd_cache.complex_cut_energies[p.id] * Into::<T>::into(-2.0))
                .powi(ltd_cache.propagator_powers[p.id] as i32);
        }
        let mut skip: bool;
        let mut ellipsoid_count;
        for mono in self.partial_fractioning_element.iter() {
            // check if the all the ellipsoids exists otherwise move to the next entry
            skip = false;
            ellipsoid_count = 0;
            for (i1, i2) in mono.ellipsoids_product.iter() {
                if *i1 < min_index|| *i2 < min_index {
                    skip = true;
                    break;
                }
                pf_cache.ellipsoids_product[ellipsoid_count] = (
                    pf_cache.numerator_index_map[*i1],
                    pf_cache.numerator_index_map[*i2],
                );
                ellipsoid_count += 1;
            }
            if skip {
                continue;
            }

            pf_cache.numerator_size = 0;
            for i in mono.numerator.iter() {
                if *i >= min_index && ltd_cache.propagator_powers[pf_cache.numerator_index_map[*i]] != 0
                {
                    pf_cache.numerator[pf_cache.numerator_size] = pf_cache.numerator_index_map[*i];
                    pf_cache.numerator_size += 1;
                }
            }
//            println!(
//                "USE:: {:?} \n\t-> {:?} num{:?}",
//                mono,
//                &pf_cache.ellipsoids_product[0..ellipsoid_count],
//                &pf_cache.numerator[0..pf_cache.numerator_size]
//            );
            result += PartialFractioningMonomial::evaluate_ellipsoids_product(
                &pf_cache.ellipsoids_product[0..ellipsoid_count],
                ltd_cache,
            )
            .inv()
                * PartialFractioningMonomial::evaluate_numerator(
                    &pf_cache.numerator[0..pf_cache.numerator_size],
                    ltd_numerator,
                    ll,
                    ltd_cache,
                );
        }
        return result / norm;
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

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

}
