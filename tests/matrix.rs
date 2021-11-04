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
    let i3 = Matrix::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);

    let res = i3.mul(&i3);

    assert!(res.isclose(&i3));
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
    let expected = Matrix::new([[26.0, 16.0, 28.0], [56.0, 40.0, 64.0], [86.0, 64.0, 100.0]]);

    assert!(res.isclose(&expected));
}

#[test]
fn multiply_matrices_different_shapes() {
    let mtx1 = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    let mtx2 = Matrix::new([[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]]);

    let res = mtx1.mul(&mtx2);
    let expected = Matrix::new([[140.0, 146.0], [320.0, 335.0]]);

    assert!(res.isclose(&expected));
}

#[test]
fn get_view_from_matrix() {
    let mtx = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
    let slice = mtx.view([0..2, 1..3]);
    let slice_mtx: Matrix<f64, 2, 2> = slice.into();

    let expected = Matrix::new([[2.0, 3.0], [5.0, 6.0]]);
    assert!(slice_mtx.isclose(&expected));
}

#[test]
fn view_iterator_column() {
    let mtx = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    let slice = mtx.view::<3, 1>([0..3, 1..2]);

    let mut iter = slice.iter();

    assert_eq!(Some(&2.0), iter.next());
    assert_eq!(Some(&5.0), iter.next());
    assert_eq!(Some(&8.0), iter.next());
    assert_eq!(None, iter.next());
}

#[test]
fn view_iterator_matrix() {
    let mtx = Matrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    let slice = mtx.view::<2, 2>([1..3, 0..2]);

    let mut iter = slice.iter();

    assert_eq!(Some(&4.0), iter.next());
    assert_eq!(Some(&5.0), iter.next());
    assert_eq!(Some(&7.0), iter.next());
    assert_eq!(Some(&8.0), iter.next());
    assert_eq!(None, iter.next());
}
