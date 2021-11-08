use lager::{isclose::IsClose, Matrix};

#[test]
fn index_matrices() {
    let mtx = Matrix::new([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]);

    assert_eq!(mtx[[0, 0]], 1.);
    assert_eq!(mtx[[0, 2]], 3.);
    assert_eq!(mtx[[2, 0]], 7.);
    assert_eq!(mtx[[2, 2]], 9.);
}

#[test]
fn check_matrices_are_close() {
    let mtx1 = Matrix::new([[4.4, 7.3, 2.5], [2.4, 6.4, 2.2], [9.5, 3.5, 4.6]]);
    let mtx2 = Matrix::new([
        [4.40004, 7.299998, 2.50002],
        [2.39998, 6.400001, 2.19998],
        [9.50006, 3.499997, 4.60004],
    ]);

    assert!(mtx1.isclose(&mtx2));
}

#[test]
fn multiply_3x3_identity_matrices() {
    let i3 = Matrix::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);

    let res = i3.mul(&i3);

    assert!(res.isclose(&i3));
}

#[test]
fn multiply_2x2_matrices_elementwise() {
    let mtx1 = Matrix::new([[1., 2.], [3., 4.]]);
    let mtx2 = Matrix::new([[4., 6.], [7., 5.]]);

    let res = mtx1 * mtx2;

    assert!(res.isclose(&Matrix::new([[4., 12.], [21., 20.]])));
}

#[test]
fn multiply_3x3_matrices() {
    let mtx1 = Matrix::new([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]);
    let mtx2 = Matrix::new([[1., 2., 1.], [2., 4., 6.], [7., 2., 5.]]);

    let res = mtx1.mul(&mtx2);
    let expected = Matrix::new([[26., 16., 28.], [56., 40., 64.], [86., 64., 100.]]);

    assert!(res.isclose(&expected));
}

#[test]
fn multiply_matrices_different_shapes() {
    let mtx1 = Matrix::new([[1., 2., 3.], [4., 5., 6.]]);
    let mtx2 = Matrix::new([[10., 11.], [20., 21.], [30., 31.]]);

    let res = mtx1.mul(&mtx2);
    let expected = Matrix::new([[140., 146.], [320., 335.]]);

    assert!(res.isclose(&expected));
}

#[test]
fn get_view_from_matrix() {
    let mtx = Matrix::new([[1., 2., 3.], [4., 5., 6.]]);
    let slice = mtx.view([0..2, 1..3]);
    let slice_mtx: Matrix<f64, 2, 2> = slice.into();

    let expected = Matrix::new([[2., 3.], [5., 6.]]);
    assert!(slice_mtx.isclose(&expected));
}

#[test]
fn view_iterator_column() {
    let mtx = Matrix::new([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]);
    let slice = mtx.view::<3, 1>([0..3, 1..2]);

    let mut iter = slice.iter();

    assert_eq!(Some(&2.), iter.next());
    assert_eq!(Some(&5.), iter.next());
    assert_eq!(Some(&8.), iter.next());
    assert_eq!(None, iter.next());
}

#[test]
fn view_iterator_matrix() {
    let mtx = Matrix::new([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]]);
    let slice = mtx.view::<2, 2>([1..3, 0..2]);

    let mut iter = slice.iter();

    assert_eq!(Some(&4.), iter.next());
    assert_eq!(Some(&5.), iter.next());
    assert_eq!(Some(&7.), iter.next());
    assert_eq!(Some(&8.), iter.next());
    assert_eq!(None, iter.next());
}

#[test]
fn row_echelon() {
    let mut mtx = Matrix::new([
        [2., 4., 1., 1., 0., 0.],
        [-1., 1., -1., 0., 1., 0.],
        [1., 4., 0., 0., 0., 1.],
    ]);

    mtx.into_row_echelon();
    let expected = Matrix::new([
        [2., 4., 1., 1., 0., 0.],
        [0., 3., -0.5, 0.5, 1., 0.],
        [0., 0., -1. / 6., -5. / 6., -2. / 3., 1.],
    ]);

    assert!(mtx.isclose(&expected));
}

#[test]
fn reduced_row_echelon() {
    let mut mtx = Matrix::new([
        [2., 4., 1., 1., 0., 0.],
        [-1., 1., -1., 0., 1., 0.],
        [1., 4., 0., 0., 0., 1.],
    ]);

    mtx.into_reduced_row_echelon();
    let expected = Matrix::new([
        [1., 0., 0., -4., -4., 5.],
        [0., 1., 0., 1., 1., -1.],
        [0., 0., 1., 5., 4., -6.],
    ]);

    assert!(mtx.isclose(&expected));
}

#[test]
fn matrix_inverse() {
    let mtx = Matrix::new([[2., 4., 1.], [-1., 1., -1.], [1., 4., 0.]]);

    let result = mtx.inv();
    let expected = Matrix::new([[-4., -4., 5.], [1., 1., -1.], [5., 4., -6.]]);

    assert!(result.isclose(&expected));
}

#[test]
fn hstack() {
    let mtx1: Matrix<f64, 3, 4> =
        Matrix::new([[2., 4., 1., 7.], [-1., 1., -1., 2.], [1., 4., 0., 9.]]);
    let mtx2: Matrix<f64, 3, 3> = Matrix::new([[-4., -4., 5.], [1., 1., -1.], [5., 4., -6.]]);

    let result: Matrix<f64, 3, 7> = mtx1.hstack(mtx2);
    let expected = Matrix::new([
        [2., 4., 1., 7., -4., -4., 5.],
        [-1., 1., -1., 2., 1., 1., -1.],
        [1., 4., 0., 9., 5., 4., -6.],
    ]);

    assert!(result.isclose(&expected));
}

#[test]
fn vstack() {
    let mtx1: Matrix<f64, 3, 3> = Matrix::new([[2., 4., 1.], [-1., 1., -1.], [1., 4., 0.]]);
    let mtx2: Matrix<f64, 4, 3> =
        Matrix::new([[-4., -4., 5.], [1., 1., -1.], [5., 4., -6.], [7., 2., 9.]]);

    let result: Matrix<f64, 7, 3> = mtx1.vstack(mtx2);
    let expected = Matrix::new([
        [2., 4., 1.],
        [-1., 1., -1.],
        [1., 4., 0.],
        [-4., -4., 5.],
        [1., 1., -1.],
        [5., 4., -6.],
        [7., 2., 9.],
    ]);

    assert!(result.isclose(&expected));
}

#[test]
fn zeros() {
    let result = Matrix::zeros();
    let expected = Matrix::new([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]]);

    assert!(result.isclose(&expected));
}

#[test]
fn ones() {
    let result = Matrix::ones();
    let expected = Matrix::new([[1., 1., 1.], [1., 1., 1.], [1., 1., 1.]]);

    assert!(result.isclose(&expected));
}

#[test]
fn identity_3x3() {
    let result = Matrix::identity();
    let expected = Matrix::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);

    assert!(result.isclose(&expected));
}

#[test]
fn lu_decomposition() {
    let mtx = Matrix::new([[2., -1., -2.], [-4., 6., 3.], [-4., -2., 8.]]);

    let expected = LUDecomposition {
        l: Matrix::new([[1., 0., 0.], [-2., 1., 0.], [-2., -1., 1.]]),
        u: Matrix::new([[2., -1., -2.], [0., 4., -1.], [0., 0., 3.]]),
        p: Matrix::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
    };

    let result = mtx.lu();

    assert!(result.l.isclose(&expected.l));
    assert!(result.u.isclose(&expected.u));
    assert!(result.p.isclose(&expected.p));
}

#[test]
fn lu_decomposition_unknown_result() {
    // This test gives a result that doesn't match other implementations,
    // but still seems to give a valid result. Probably something to do with
    // the choice of pivots.

    let mtx = Matrix::new([[10., -7., 0.], [-3., 2., 6.], [5., -1., 5.]]);

    // What other tools give as the output, not used here
    // let expected = LUDecomposition {
    //     l: Matrix::new([[1., 0., 0.], [0.5, 1., 0.], [-0.3, -0.04, 1.]]),
    //     u: Matrix::new([[10., -7., 0.], [0., 2.5, 5.], [0., 0., 6.2]]),
    //     p: Matrix::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
    // };

    let result = mtx.lu();

    assert!(mtx.isclose(&result.p.mul(&result.l.mul(&result.u))));
}

#[test]
fn lu_decomposition_with_permuations() {
    // This test requires permutations otherwise NaNs result in the
    // decomposition
    let mtx = Matrix::new([
        [0., 2., 8., 6.],
        [0., 0., 1., 2.],
        [0., 1., 0., 1.],
        [3., 7., 1., 0.],
    ]);

    let result = mtx.clone().lu();

    assert!(mtx.isclose(&result.p.mul(&result.l.mul(&result.u))));
}

#[test]
fn count_3x3() {
    let mtx = Matrix::new([[2., 4., 1.], [0., 1., -1.], [1., 4., 0.]]);
    assert_eq!(mtx.count(), 7);
}

#[test]
fn diag_3x3() {
    let mtx = Matrix::new([[2., 4., 1.], [0., 1., -1.], [1., 4., 0.]]);
    let result = mtx.diag();

    let expected = Matrix::new([[2.], [1.], [0.]]);

    assert!(expected.isclose(&result))
}
