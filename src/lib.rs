use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Range, Sub};

use abs::Abs;
use matrixview::MatrixView;

mod abs;
mod matrixview;

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize> {
    pub values: [[T; N]; M],
}

pub type Vector<T, const M: usize> = Matrix<T, M, 1>;

impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Neg<Output = T>
        + Default
        + AddAssign
        + Abs
        + PartialOrd
        + From<f32>,
    f64: From<T>,
{
    pub fn new(array: [[T; N]; M]) -> Self {
        Self { values: array }
    }

    pub fn mul<const O: usize>(&self, rhs: &Matrix<T, N, O>) -> Matrix<T, M, O> {
        let mut res = [[Default::default(); O]; M];

        // For each element in the result
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

impl<T, const M: usize, const N: usize> Add for Matrix<T, M, N>
where
    T: Copy + Add<Output = T> + Default,
{
    type Output = Matrix<T, M, N>;

    fn add(self, other: Self) -> Self {
        let mut res = [[Default::default(); N]; M];

        // For each element
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
        + Default
        + Abs
        + AddAssign
        + PartialOrd
        + From<f32>,
    f64: From<T>,
{
    fn from(view: MatrixView<T, M, N, R, C>) -> Self {
        // Crete a new empty matrix
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

