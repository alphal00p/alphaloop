use arrayvec::ArrayVec;

type Complex = num::Complex<f64>;
const MAX_DIM: usize = 128;

#[inline]
/// Invert with better precision
pub fn finv(c: Complex) -> Complex {
    let norm = c.norm();
    c.conj() / norm / norm
}

/// Calculate the determinant of any complex-valued input matrix using LU-decomposition.
/// Original C-code by W. Gong and D.E. Soper.
pub fn determinant(bb: &ArrayVec<[Complex; MAX_DIM]>) -> Complex {
    // Define matrix related variables.
    let dimension = (bb.len() as f64).sqrt() as usize;
    let mut determinant = Complex::new(1.0, 0.0);
    let mut indx = [0; MAX_DIM];
    let mut d = 1; // initialize parity parameter

    // Inintialize the matrix to be decomposed with the transferred matrix b.
    let mut aa = bb.clone();

    // Define parameters used in decomposition.
    let mut imax = 0;
    let mut flag = 1;
    let mut dumc;
    let mut sum;

    let mut aamax;
    let mut dumr;
    let mut vv = [0.0; MAX_DIM];

    // Get the implicit scaling information.
    for i in 0..dimension {
        aamax = 0.0;
        for j in 0..dimension {
            let r = aa[i * dimension + j].norm_sqr();
            if r > aamax {
                aamax = r;
            }
        }
        // Set a flag to check if the determinant is zero.
        if aamax == 0.0 {
            flag = 0;
        }
        // Save the scaling.
        vv[i] = 1.0 / aamax;
    }
    if flag == 1 {
        for j in 0..dimension {
            for i in 0..j {
                sum = aa[i * dimension + j];
                for k in 0..i {
                    sum = sum - aa[i * dimension + k] * aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
            }
            //Initialize for the search for largest pivot element.
            aamax = 0.0;
            for i in j..dimension {
                sum = aa[i * dimension + j];
                for k in 0..j {
                    sum = sum - aa[i * dimension + k] * aa[k * dimension + j];
                }
                aa[i * dimension + j] = sum;
                // Figure of merit for the pivot.
                dumr = vv[i] * sum.norm_sqr();
                // Is it better than the best so far?
                if dumr >= aamax {
                    imax = i;
                    aamax = dumr;
                }
            }
            // See if we need to interchange rows.
            if j != imax {
                for k in 0..dimension {
                    dumc = aa[imax * dimension + k];
                    aa[imax * dimension + k] = aa[j * dimension + k];
                    aa[j * dimension + k] = dumc;
                }
                // Change the parity of d.
                d = -d;
                // Interchange the scale factor.
                vv[imax] = vv[j];
            }
            indx[j] = imax;
            if j + 1 != dimension {
                dumc = 1.0 / aa[j * dimension + j];
                for i in j + 1..dimension {
                    aa[i * dimension + j] = aa[i * dimension + j] * dumc;
                }
            }
        }
    }
    // Calculate the determinant using the decomposed matrix.
    if flag == 0 {
        determinant = Complex::new(0.0, 0.0);
    } else {
        // Multiply the diagonal elements.
        for diagonal in 0..dimension {
            determinant = determinant * aa[diagonal * dimension + diagonal];
        }
        determinant = d as f64 * determinant;
    }
    determinant
}
