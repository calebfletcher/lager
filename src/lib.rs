#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::{
    cmp,
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, Neg, Range, Sub, SubAssign},
};

use abs::Abs;
use isclose::IsClose;
use matrixview::MatrixView;

mod abs;
mod alg;
pub mod isclose;
mod matrixview;

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize> {
    pub values: [[T; N]; M],
}

#[derive(Debug, Clone)]
pub struct LUDecomposition<T, const M: usize, const N: usize> {
    pub l: Matrix<T, M, M>, //  Not a typo, see https://en.wikipedia.org/wiki/LU_decomposition#Rectangular_matrices
    pub u: Matrix<T, M, N>,
    pub p: Matrix<T, M, N>,
}

pub struct MatrixIter<'a, T, const M: usize, const N: usize> {
    mtx: &'a Matrix<T, M, N>,
    i: usize,
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
        + IsClose
        + Debug
        + PartialOrd
        + From<f64>,
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

    pub fn t(&self) -> Matrix<T, N, M> {
        let mut mtx = Matrix::zeros();

        for row in 0..M {
            for col in 0..N {
                mtx[[col, row]] = self[[row, col]];
            }
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
        // For each element in the matrices
        for i in 0..M {
            for j in 0..N {
                // Compare difference to the threshold
                if !self[[i, j]].isclose(other[[i, j]]) {
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

    pub fn iter(&self) -> MatrixIter<T, M, N> {
        MatrixIter { mtx: self, i: 0 }
    }

    pub fn prod(&self) -> T {
        self.iter().fold(1.0.into(), |acc: T, x| acc * *x)
    }

    pub fn sum(&self) -> T {
        self.iter().fold(0.0.into(), |acc: T, x| acc + *x)
    }

    pub fn is_upper_triangular(&self) -> bool {
        for row in 0..M {
            for col in 0..row {
                dbg!(row, col);
                if !self[[row, col]].isclose(0.0.into()) {
                    return false;
                }
            }
        }
        true
    }

    pub fn is_lower_triangular(&self) -> bool {
        for row in 0..M {
            for col in (row + 1)..N {
                dbg!(row, col);
                if !self[[row, col]].isclose(0.0.into()) {
                    return false;
                }
            }
        }
        true
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

impl<T, const M: usize> Index<usize> for Vector<T, M> {
    type Output = T;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.values[idx][0]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<[usize; 2]> for Matrix<T, M, N> {
    fn index_mut(&mut self, idx: [usize; 2]) -> &mut Self::Output {
        &mut self.values[idx[0]][idx[1]]
    }
}

impl<T, const M: usize> IndexMut<usize> for Vector<T, M> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.values[idx][0]
    }
}

impl<'a, T, const M: usize, const N: usize> Iterator for MatrixIter<'a, T, M, N> {
    type Item = &'a T;
    fn next(&mut self) -> Option<Self::Item> {
        let current = self.i;

        if self.i == M * N {
            return None;
        }

        self.i += 1;

        Some(&self.mtx[[current / N, current % N]])
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
        + IsClose
        + Debug
        + AddAssign
        + SubAssign
        + DivAssign
        + PartialOrd
        + From<f64>,
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
