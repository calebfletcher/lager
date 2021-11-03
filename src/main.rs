fn main() {
    let mtx = Matrix::new();
    dbg!("{}", mtx);
}

#[derive(Debug)]
struct Matrix {
    values: [[f32; 3]; 3],
}

impl Matrix {
    fn new() -> Self {
        Self {
            values: [[0.0; 3]; 3],
        }
    }
}
