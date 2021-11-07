#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::{
    cmp,
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, Neg, Range, Sub, SubAssign},
};

use abs::Abs;
use matrixview::MatrixView;

mod abs;
mod matrixview;

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize> {
    pub values: [[T; N]; M],
}

pub struct LUDecomposition<T, const M: usize, const N: usize> {
    pub l: Matrix<T, M, M>, //  Not a typo, see https://en.wikipedia.org/wiki/LU_decomposition#Rectangular_matrices
    pub u: Matrix<T, M, N>,
    pub p: Matrix<T, M, N>,
}

pub type Vector<T, const M: usize> = Matrix<T, M, 1>;

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
        + Debug
        + PartialOrd
        + From<f32>,
    f64: From<T>,
{
    pub fn new(array: [[T; N]; M]) -> Self {
        Self { values: array }
    }

    pub fn zeros() -> Self {
        Self::new([[0.0.into(); N]; M])
    }

    pub fn ones() -> Self {
        Self::new([[1.0.into(); N]; M])
    }

    pub fn identity() -> Self {
        let mut mtx = Self::zeros();

        let min_axis = cmp::min(M, N);
        for i in 0..min_axis {
            mtx[[i, i]] = 1.0.into();
        }

        mtx
    }

    pub fn mul<const O: usize>(&self, rhs: &Matrix<T, N, O>) -> Matrix<T, M, O> {
        let mut res = [[Default::default(); O]; M];

        // For each element in the result
        #[allow(clippy::needless_range_loop)]
        for i in 0..M {
            for j in 0..O {
                // For each value contributing to the resulting element
                for k in 0..N {
                    res[i][j] += self.values[i][k] * rhs.values[k][j];
                }
            }
        }

        Matrix::<T, M, O> { values: res }
    }

    pub fn isclose(&self, other: &Matrix<T, M, N>) -> bool {
        let rt: T = 1e-05.into();
        let at: T = 1e-08.into();

        // For each element in the matrices
        for i in 0..M {
            for j in 0..N {
                // Compare difference to the threshold
                let threshold = rt * self[[i, j]].abs() + at;
                if (self[[i, j]] - other[[i, j]]).abs() > threshold {
                    return false;
                }
            }
        }
        true
    }

    pub fn view<const R: usize, const C: usize>(
        &self,
        bounds: [Range<usize>; 2],
    ) -> MatrixView<T, M, N, R, C> {
        MatrixView::new(self, bounds)
    }

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

    fn swap_rows(&mut self, row1: usize, row2: usize) {
        if row1 == row2 {
            return;
        }
        self.values.swap(row1, row2);
    }

    pub fn hstack<const C2: usize>(&self, other: Matrix<T, M, C2>) -> Matrix<T, M, { N + C2 }> {
        let mut mtx: Matrix<T, M, { N + C2 }> = self.view([0..M, 0..N]).into();

        // For each element in the other matrix
        for row in 0..M {
            for col in 0..C2 {
                // Copy the element from the other matrix to the new matrix
                mtx[[row, col + N]] = other[[row, col]];
            }
        }
        mtx
    }

    pub fn vstack<const R2: usize>(&self, other: Matrix<T, R2, N>) -> Matrix<T, { M + R2 }, N> {
        let mut mtx: Matrix<T, { M + R2 }, N> = self.view([0..M, 0..N]).into();

        // For each element in the other matrix
        for row in 0..R2 {
            for col in 0..N {
                // Copy the element from the other matrix to the new matrix
                mtx[[row + M, col]] = other[[row, col]];
            }
        }
        mtx
    }

    pub fn lu(&self) -> LUDecomposition<T, M, N> {
        let mut mtx = self.clone();

        let mut perms = Matrix::identity();
        let mut lower = Matrix::zeros();
        let mut upper = Matrix::zeros();

        for col in 0..N {
            let mut next_nonzero_row = col;
            for row in col..M {
                if mtx[[row, col]] != 0.0.into() {
                    next_nonzero_row = row;
                    break;
                }
            }

            mtx.swap_rows(col, next_nonzero_row);
            perms.swap_rows(col, next_nonzero_row);
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

    pub fn diag(&self) -> Matrix<T, M, 1> {
        let mut vec = Vector::zeros();
        for i in 0..M {
            vec[i] = self[[i, i]];
        }
        vec
    }

    pub fn count(&self) -> usize {
        let mut count = 0;
        for i in 0..M {
            for j in 0..N {
                if self[[i, j]] != 0.0.into() {
                    count += 1;
                }
            }
        }
        count
    }
}

fn argmax<T: PartialOrd>(slice: impl Iterator<Item = T>) -> Option<usize> {
    slice
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value0.partial_cmp(value1).expect("nan found in argmax"))
        .map(|(idx, _)| idx)
}

impl<T, const M: usize, const N: usize> Add for Matrix<T, M, N>
where
    T: Copy + Add<Output = T> + Default,
{
    type Output = Matrix<T, M, N>;

    fn add(self, other: Self) -> Self {
        let mut res = [[Default::default(); N]; M];

        // For each element
        #[allow(clippy::needless_range_loop)]
        for row in 0..M {
            for col in 0..N {
                // Add each element
                res[row][col] = self.values[row][col] + other.values[row][col];
            }
        }

        Matrix { values: res }
    }
}

impl<T, const M: usize, const N: usize> Mul for Matrix<T, M, N>
where
    T: Copy + Mul<Output = T> + Default,
{
    type Output = Matrix<T, M, N>;

    fn mul(self, rhs: Self) -> Matrix<T, M, N> {
        let mut res = [[Default::default(); N]; M];

        // For each element
        #[allow(clippy::needless_range_loop)]
        for row in 0..M {
            for col in 0..N {
                // Multiply the elements together
                res[row][col] = self.values[row][col] * rhs.values[row][col];
            }
        }

        Matrix { values: res }
    }
}

impl<T, const M: usize, const N: usize> Index<[usize; 2]> for Matrix<T, M, N> {
    type Output = T;

    fn index(&self, idx: [usize; 2]) -> &Self::Output {
        &self.values[idx[0]][idx[1]]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<[usize; 2]> for Matrix<T, M, N> {
    fn index_mut(&mut self, idx: [usize; 2]) -> &mut Self::Output {
        &mut self.values[idx[0]][idx[1]]
    }
}

impl<T, const M: usize, const N: usize, const R: usize, const C: usize>
    From<MatrixView<'_, T, M, N, R, C>> for Matrix<T, R, C>
where
    T: Sub<Output = T>
        + Copy
        + Mul<Output = T>
        + Add<Output = T>
        + Neg<Output = T>
        + Div<Output = T>
        + Default
        + Abs
        + Debug
        + AddAssign
        + SubAssign
        + DivAssign
        + PartialOrd
        + From<f32>,
    f64: From<T>,
{
    fn from(view: MatrixView<T, M, N, R, C>) -> Self {
        // Create a new empty matrix
        let mut mtx = Matrix::new([[0.0.into(); C]; R]);

        // For each element in the view
        for row in view.bounds[0].clone() {
            for col in view.bounds[1].clone() {
                // Copy the element from the view to the new matrix
                mtx[[row - view.bounds[0].start, col - view.bounds[1].start]] =
                    view.mtx[[row, col]];
            }
        }
        mtx
    }
}
