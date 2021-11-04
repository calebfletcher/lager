use std::ops::{Index, Range};

use crate::Matrix;

pub struct MatrixView<'a, T, const M: usize, const N: usize, const R: usize, const C: usize> {
    pub mtx: &'a Matrix<T, M, N>,
    pub bounds: [Range<usize>; 2],
}

impl<T, const M: usize, const N: usize, const R: usize, const C: usize>
    MatrixView<'_, T, M, N, R, C>
{
    pub fn new(mtx: &Matrix<T, M, N>, bounds: [Range<usize>; 2]) -> MatrixView<T, M, N, R, C> {
        MatrixView { mtx, bounds }
    }
}

impl<'a, T, const M: usize, const N: usize, const R: usize, const C: usize> Index<[usize; 2]>
    for MatrixView<'a, T, M, N, R, C>
{
    type Output = T;

    fn index(&self, idx: [usize; 2]) -> &Self::Output {
        &self.mtx.values[self.bounds[0].start + idx[0]][self.bounds[1].start + idx[1]]
    }
}
