use num::Complex;
use vector::LorentzVector;
use FloatLike;

pub struct GammaChain<'a, T: FloatLike> {
    vbar: &'a [Complex<T>],
    u: &'a [Complex<T>],
    vectors: &'a [LorentzVector<Complex<T>>],
    indices: &'a [i8],
    repeated_indices_values: Vec<i8>,
    repeated_indices: Vec<i8>,
}

impl<'a, T: FloatLike> GammaChain<'a, T> {
    pub fn new(
        vbar: &'a [Complex<T>],
        u: &'a [Complex<T>],
        indices: &'a [i8],
        vectors: &'a [LorentzVector<Complex<T>>],
    ) -> Result<GammaChain<'a, T>, &'static str> {
        //Inintialize the value of all repeated indices values to -1
        //i.e not used
        let mut repeated_indices_values = vec![-1; indices.len()];

        //List all repeated repeated indices
        let mut repeated_indices = Vec::new();
        for i in 0..indices.len() {
            if indices[i] < 0 && !repeated_indices.contains(&indices[i]) {
                //Add it to the list
                repeated_indices.push(indices[i]);
                repeated_indices_values[(indices[i].abs() - 1) as usize] = 0;
            }
        }
        //Resize to the needed dimension
        repeated_indices_values.truncate(repeated_indices.len());

        //TODO: add check to see if all the repeated indices use consecutive negative numbers
        //TODO: add check to see if there are no 0 indices
        return Ok(GammaChain {
            vbar: vbar,
            u: u,
            vectors: vectors,
            indices: indices,
            repeated_indices_values: repeated_indices_values,
            repeated_indices: repeated_indices,
        });
    }

    pub fn gamma_multiply(
        flow: &[Complex<T>],
        gamma_index: i8,
        factor: Complex<T>,
        res_flow: &mut [Complex<T>],
    ) -> Result<(), &'static str> {
        //  Variables declaration
        let ii = Complex::new(T::zero(), T::one());

        // Implementation
        match gamma_index {
            0 => {
                res_flow[0] += flow[2] * factor;
                res_flow[1] += flow[3] * factor;
                res_flow[2] += flow[0] * factor;
                res_flow[3] += flow[1] * factor;
            }
            1 => {
                res_flow[0] -= flow[3] * factor;
                res_flow[1] -= flow[2] * factor;
                res_flow[2] += flow[1] * factor;
                res_flow[3] += flow[0] * factor;
            }
            2 => {
                res_flow[0] -= ii * flow[3] * factor;
                res_flow[1] += ii * flow[2] * factor;
                res_flow[2] += ii * flow[1] * factor;
                res_flow[3] -= ii * flow[0] * factor;
            }
            3 => {
                res_flow[0] -= flow[2] * factor;
                res_flow[1] += flow[3] * factor;
                res_flow[2] += flow[0] * factor;
                res_flow[3] -= flow[1] * factor;
            }
            5 => {
                res_flow[0] -= flow[0] * factor;
                res_flow[1] -= flow[1] * factor;
                res_flow[2] += flow[2] * factor;
                res_flow[3] += flow[3] * factor;
            }
            _ => return Err("Unknown"),
        }
        return Ok(());
    }

    fn compute_chain_with_dummy_indices_fixed(&mut self) -> Result<Complex<T>, &'static str> {
        //Variable declaration
        let one = Complex::new(T::one(), T::zero());
        let ii = Complex::new(T::zero(), T::one());
        let zero = Complex::new(T::zero(), T::zero());
        //Follow the spinor flow
        let mut flow: [Complex<T>; 4] = [self.vbar[0], self.vbar[1], self.vbar[2], self.vbar[3]];

        for index in self.indices {
            let mut next_flow = [zero; 4];
            let pos = (index.abs() - 1) as usize;
            if *index > 0 {
                //Contraction with an external vector
                GammaChain::gamma_multiply(&flow, 0, self.vectors[pos][0], &mut next_flow).unwrap();
                for j in 1..4 {
                    GammaChain::gamma_multiply(
                        &flow,
                        j as i8,
                        -self.vectors[pos][j],
                        &mut next_flow,
                    )
                    .unwrap();
                }
            } else {
                // Internal Summation
                if self.repeated_indices_values[pos] == 0 {
                    GammaChain::gamma_multiply(
                        &flow,
                        self.repeated_indices_values[pos],
                        one,
                        &mut next_flow,
                    )
                    .unwrap();
                } else {
                    GammaChain::gamma_multiply(
                        &flow,
                        self.repeated_indices_values[pos],
                        ii,
                        &mut next_flow,
                    )
                    .unwrap();
                }
            }
            // Update running flow
            flow = next_flow;
        }

        //Build and return result
        let mut res = zero;
        for i in 0..4 {
            res += flow[i] * self.u[i];
        }
        return Ok(res);
    }

    fn compute_chain_element(
        &mut self,
        repeated_index_pointer: usize,
    ) -> Result<Complex<T>, &'static str> {
        let mut res = Complex::new(T::zero(), T::zero());
        // If the repeated_indicex_pointer is beyond the number of repeated indices then return
        if repeated_index_pointer >= self.repeated_indices.len() {
            res += self.compute_chain_with_dummy_indices_fixed().unwrap();
        } else {
            //Then perfrom the summation over the dummy index
            for i in 0..4 {
                let pos = (self.repeated_indices[repeated_index_pointer].abs() - 1) as usize;
                self.repeated_indices_values[pos] = i;
                res += self
                    .compute_chain_element(repeated_index_pointer + 1)
                    .unwrap();
            }
        }

        return Ok(res);
    }

    pub fn compute_chain(&mut self) -> Result<Complex<T>, &'static str> {
        let res = self.compute_chain_element(0).unwrap();
        return Ok(res);
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[allow(non_snake_case, dead_code)]
    #[test]
    fn test_gamma() {
        //Initialize some vectors contracted with matrix 1 and 2, respectively
        let mut vectors = Vec::new();
        vectors.push(LorentzVector::from_args(
            Complex::new(1.0, 0.0),
            Complex::new(1.1, 0.0),
            Complex::new(1.2, 0.0),
            Complex::new(1.3, 0.0),
        ));
        vectors.push(LorentzVector::from_args(
            Complex::new(2.0, 0.0),
            Complex::new(2.1, 0.0),
            Complex::new(2.2, 0.0),
            Complex::new(2.3, 0.0),
        ));

        let zero = Complex::default();
        let ii = Complex::new(0.0, 1.0) / (2.0 as f64).sqrt();
        let one = Complex::new(1.0, 0.0) / (2.0 as f64).sqrt();

        //Some beatiful contraction of dummy indices together with some slashed momenta
        let indices = vec![1, -1, -2, 2, -1, -2, -3, -4, -5, -6, -5, -4, -6, -3];

        //Define vbar and u
        let vbar = vec![
            Complex::new(3.0, 0.0),
            Complex::new(3.1, 0.0),
            Complex::new(3.2, 0.0),
            Complex::new(3.3, 0.0),
        ];
        let u = vec![
            Complex::new(4.0, 0.0),
            Complex::new(4.1, 0.0),
            Complex::new(4.2, 0.0),
            Complex::new(4.3, 0.0),
        ];
        //Change to Weyl Representation
        //T = 1/Sqrt(2) ( Id x Id + i \sigma_2 x I )
        let T = [
            [one, zero, one, zero],
            [zero, one, zero, one],
            [-one, zero, one, zero],
            [zero, -one, zero, one],
        ];

        let Tt = [
            [one, zero, -one, zero],
            [zero, one, zero, -one],
            [one, zero, one, zero],
            [zero, one, zero, one],
        ];
        let mut vbarW = Vec::new();
        for t in Tt.iter() {
            let mut res = Complex::new(0.0, 0.0);
            for (x, y) in vbar.iter().zip(t) {
                res += x * y;
            }
            vbarW.push(res);
        }

        let mut uW = Vec::new();
        for t in Tt.iter() {
            let mut res = Complex::new(0.0, 0.0);
            for (x, y) in u.iter().zip(t) {
                res += x * y;
            }
            uW.push(res);
        }
        println!("{:?}", uW);
        println!("{:?}", vbarW);
        //Initialize gamma chain
        let mut chain = GammaChain::new(&vbarW, &uW, &indices, &vectors).unwrap();

        //Solve
        let result = Complex::new(-78354.8416, 1312.256);
        let chain_result = chain.compute_chain().unwrap();
        println!("Computed result: {:?}", chain_result);
        println!("Exact result: {:?}", result);

        assert!((chain_result - result).norm() < 1e-10);
    }
}
