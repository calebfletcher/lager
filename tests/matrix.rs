use lager::Matrix;

#[test]
fn index_matrices() {
    let mtx = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);

    assert_eq!(mtx[[0, 0]], 1.0);
    assert_eq!(mtx[[0, 2]], 3.0);
    assert_eq!(mtx[[2, 0]], 7.0);
    assert_eq!(mtx[[2, 2]], 9.0);
}

#[test]
fn matrices_are_close() {
    let mtx1 = Matrix::new([[4.4, 7.3, 2.5], [2.4, 6.4, 2.2], [9.5, 3.5, 4.6]]);
    let mtx2 = Matrix::new([
        [4.40004, 7.299998, 2.50002],
        [2.39998, 6.400001, 2.19998],
        [9.50006, 3.499997, 4.60004],
    ]);
    assert!(mtx1.isclose(&mtx2));
}

#[test]
fn multiply_identity_matrices() {
    let i3 = Matrix::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let res = i3.mul(&i3);

    assert_eq!(res, i3);
}

#[test]
fn multiply_2x2_matrices_elementwise() {
    let mtx1 = Matrix::new([[1.0, 2.0], [3.0, 4.0]]);
    let mtx2 = Matrix::new([[4.0, 6.0], [7.0, 5.0]]);
    let res = mtx1 * mtx2;

    assert!(res.isclose(&Matrix::new([[4.0, 12.0], [21.0, 20.0]])));
}

#[test]
fn multiply_3x3_matrices() {
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
