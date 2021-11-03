use criterion::{black_box, criterion_group, criterion_main, Criterion};
use matrix::Matrix;

fn add(a: Matrix, b: Matrix) -> Matrix {
    a + b
}

fn mul(a: Matrix, b: Matrix) -> Matrix {
    a.mul(&b)
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("mat<3,3> add", |b| {
        b.iter(|| {
            add(
                black_box(Matrix::new([
                    [4.4, 7.3, 2.5],
                    [2.4, 6.4, 2.2],
                    [9.5, 3.5, 4.6],
                ])),
                black_box(Matrix::new([
                    [4.4, 7.3, 2.5],
                    [2.4, 6.4, 2.2],
                    [9.5, 3.5, 4.6],
                ])),
            )
        })
    });

    c.bench_function("mat<3,3> mult", |b| {
        b.iter(|| {
            mul(
                black_box(Matrix::new([
                    [4.4, 7.3, 2.5],
                    [2.4, 6.4, 2.2],
                    [9.5, 3.5, 4.6],
                ])),
                black_box(Matrix::new([
                    [4.4, 7.3, 2.5],
                    [2.4, 6.4, 2.2],
                    [9.5, 3.5, 4.6],
                ])),
            )
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
