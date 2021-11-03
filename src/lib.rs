use std::ops::{Add, AddAssign, Mul};

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize>
where
    T: Copy + Mul<Output = T> + Default,
{
    pub values: [[T; N]; M],
}

impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where
    T: Copy + Mul<Output = T> + Default + AddAssign,
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
}

impl<T, const M: usize, const N: usize> Add for Matrix<T, M, N>
where
    T: Copy + Add<Output = T> + Default + Mul<Output = T>,
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

impl<T, const M: usize, const N: usize> Mul for &Matrix<T, M, N>
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

#[cfg(test)]
mod tests {
    use crate::Matrix;

    #[test]
    fn multiply_identity_matrices() {
        let i3 = Matrix::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        let res = i3.mul(&i3);

        assert_eq!(res, i3);
    }

    #[test]
    fn multiply_single_matrices() {
        let mtx1 = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let mtx2 = Matrix::new([[1.0, 2.0, 1.0], [2.0, 4.0, 6.0], [7.0, 2.0, 5.0]]);
        let res = mtx1.mul(&mtx2);

        assert_eq!(
            res,
            Matrix::new([[26.0, 16.0, 28.0], [56.0, 40.0, 64.0], [86.0, 64.0, 100.0]])
        );
    }

    #[test]
    fn multiply_matrices_different_shapes() {
        let mtx1 = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let mtx2 = Matrix::new([[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]]);
        let res = mtx1.mul(&mtx2);

        assert_eq!(res, Matrix::new([[140.0, 146.0], [320.0, 335.0]]));
    }
}
