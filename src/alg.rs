use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, Neg, Sub, SubAssign},
};

use crate::isclose::IsClose;
use crate::{abs::Abs, LUDecomposition};

use crate::{argmax, Matrix};

impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Neg<Output = T>
        + Default
        + AddAssign
        + SubAssign
        + DivAssign
        + Abs
        + IsClose
        + Debug
        + PartialOrd
        + From<f64>,
    f64: From<T>,
{
    pub fn into_row_echelon(&mut self) {
        let mut h = 0; /* Initialization of the pivot row */
        let mut k = 0; /* Initialization of the pivot column */

        while h < M && k < N {
            /* Find the k-th pivot: */

            let i_max = h + argmax(
                self.view::<M, 1>([h..M, k..k])
                    .iter()
                    .take(M - h) // Takes only the values within the valid range
                    .map(|el| el.abs()),
            )
            .unwrap();

            //let i_max = argmax (i = h..m, A[[i, k]].abs());
            if self[[i_max, k]] == 0.0.into() {
                /* No pivot in this column, pass to next column */
            } else {
                self.swap_rows(h, i_max);
                /* Do for all rows below pivot: */
                for i in h + 1..M {
                    let f = self[[i, k]] / self[[h, k]];
                    /* Fill with zeros the lower part of pivot column: */
                    self[[i, k]] = 0.0.into();
                    /* Do for all remaining elements in current row: */
                    for j in (k + 1)..N {
                        let el = self[[h, j]];
                        self[[i, j]] -= el * f;
                    }
                }
                /* Increase pivot row and column */
                h += 1;
            }
            k += 1;
        }
    }

    pub fn into_reduced_row_echelon(&mut self) {
        self.into_row_echelon();

        let mut h = M - 1; // Initialization of the pivot row
        let mut k = 0; // Initialization of the pivot column

        // Find first non-zero element in the last row
        while self[[h, k]] == 0.0.into() {
            k += 1;
        }

        loop {
            // Make leading coefficient 1
            let f = self[[h, k]];
            for i in k..N {
                self[[h, i]] /= f;
            }

            // Don't reduce above rows if on the first row
            if h == 0 {
                break;
            }

            // Set each element in each row above to zero
            for i in 0..h {
                // For each row above
                let f = self[[i, k]] / self[[h, k]];
                for j in k..N {
                    // For each element in the row that isn't to the left of the coefficient
                    let el = self[[h, j]];
                    self[[i, j]] -= el * f;
                }
            }
            h -= 1;
            k -= 1;
            while self[[h, k]] == 0.0.into() && k > 0 {
                k -= 1;
            }
        }
    }

    pub fn inv(&self) -> Matrix<T, M, N>
    where
        [(); N + N]: Sized,
    {
        // Stack identity matrix to the right
        let mut augmented = self.hstack(Self::identity());

        // Convert augmented matrix into reduced row echelon form
        augmented.into_reduced_row_echelon();

        // Remove identity matrix from the left
        augmented.view([0..M, N..2 * N]).into()
    }

    pub fn lu(&self) -> LUDecomposition<T, M, N> {
        let mut mtx = self.clone();

        let mut perms = Matrix::identity();
        let mut lower = Matrix::zeros();
        let mut upper = Matrix::zeros();

        for col in 0..N {
            let i_max = col
                + argmax(
                    self.view::<M, 1>([col..M, col..col])
                        .iter()
                        .take(M - col) // Takes only the values within the valid range
                        .map(|el| el.abs()),
                )
                .unwrap();

            mtx.swap_rows(col, i_max);
            perms.swap_rows(col, i_max);
        }

        // Doolittle's algorithm
        for i in 0..M {
            for j in i..N {
                upper[[i, j]] = mtx[[i, j]];
                for k in 0..i {
                    let el = upper[[k, j]];
                    upper[[i, j]] -= lower[[i, k]] * el;
                }
            }
            for j in i..N {
                lower[[j, i]] = mtx[[j, i]];
                for k in 0..i {
                    let el = lower[[j, k]];
                    lower[[j, i]] -= el * upper[[k, i]];
                }
                lower[[j, i]] /= upper[[i, i]];
            }
        }

        LUDecomposition {
            l: lower,
            u: upper,
            p: perms,
        }
    }

    pub fn det(&self) -> T {
        if M == 2 && N == 2 {
            self[[0, 0]] * self[[1, 1]] - self[[0, 1]] * self[[1, 0]]
        } else {
            let decomp = self.lu();

            let num_perms = M - decomp.p.diag().count() - 1;
            let det_pinv: T = ((-1_f64).powf(num_perms as f64)).into();
            dbg!(num_perms);

            let lower_det = decomp.l.diag().prod();
            let upper_det = decomp.u.diag().prod();

            dbg!(det_pinv, lower_det, upper_det);

            det_pinv * lower_det * upper_det
        }
    }
}
