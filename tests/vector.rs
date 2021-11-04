use lager::{Matrix, Vector};

#[test]
fn add_3x1_vectors_elementwise() {
    let vec1: Vector<f64, 3> = Vector::new([[1.0], [2.0], [3.0]]);
    let vec2: Vector<f64, 3> = Vector::new([[11.0], [21.0], [31.0]]);
    let res = vec1 + vec2;

    let expected: Vector<f64, 3> = Vector::new([[12.0], [23.0], [34.0]]);
    assert!(res.isclose(&expected));
}

#[test]
fn multiply_3x1_vectors_elementwise() {
    let vec1: Vector<f64, 3> = Vector::new([[1.0], [2.0], [3.0]]);
    let vec2: Vector<f64, 3> = Vector::new([[11.0], [21.0], [31.0]]);
    let res: Vector<f64, 3> = vec1 * vec2;

    let expected: Vector<f64, 3> = Vector::new([[11.0], [42.0], [93.0]]);
    assert!(res.isclose(&expected));
}

#[test]
fn multiply_3x3_matrix_with_3x1_vector() {
    let mtx = Matrix::new([[1.0, 0.0, -2.0], [0.0, 3.0, -1.0], [1.0, 2.0, 1.0]]);
    let vec: Vector<f64, 3> = Vector::new([[3.0], [-1.0], [4.0]]);

    let res: Vector<f64, 3> = mtx.mul(&vec);

    let expected: Vector<f64, 3> = Vector::new([[-5.0], [-7.0], [5.0]]);
    assert!(res.isclose(&expected));
}
