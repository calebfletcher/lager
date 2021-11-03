use std::ops::{Add, Mul};

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix {
    pub values: [[f32; 3]; 3],
}

impl Matrix {
    pub fn new(array: [[f32; 3]; 3]) -> Self {
        Self { values: array }
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        let mut res = [[0.0; 3]; 3];

        // For each element in the result
        for i in 0..3 {
            for j in 0..3 {
                // For each value contributing to the resulting element
                for k in 0..3 {
                    res[i][j] += self.values[i][k] * rhs.values[k][j];
                }
            }
        }

        Self { values: res }
    }
}

impl Add for Matrix {
    type Output = Matrix;

    fn add(self, other: Self) -> Matrix {
        let mut res = [[0.0; 3]; 3];

        // For each element
        for row in 0..3 {
            for col in 0..3 {
                // Add each element
                res[row][col] = self.values[row][col] + other.values[row][col];
            }
        }

        Matrix { values: res }
    }
}

impl Mul for &Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Self) -> Matrix {
        let mut res = [[0.0; 3]; 3];

        // For each element
        for row in 0..3 {
            for col in 0..3 {
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
}
