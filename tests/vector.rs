use lager::{Matrix, Vector};

#[test]
fn add_vector() {
    let vec1: Vector<f32, 3> = Vector::new([[1.0], [2.0], [3.0]]);
    let vec2: Vector<f32, 3> = Vector::new([[11.0], [21.0], [31.0]]);
    let res = vec1 + vec2;

    let expected: Vector<f32, 3> = Vector::new([[12.0], [23.0], [34.0]]);
    assert!(res.isclose(&expected));
}

#[test]
fn mult_vector_elementwise() {
    let vec1: Vector<f32, 3> = Vector::new([[1.0], [2.0], [3.0]]);
    let vec2: Vector<f32, 3> = Vector::new([[11.0], [21.0], [31.0]]);
    let res: Vector<f32, 3> = vec1 * vec2;

    let expected: Vector<f32, 3> = Vector::new([[11.0], [42.0], [93.0]]);
    assert!(res.isclose(&expected));
}
