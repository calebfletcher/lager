use std::ops::{Add, Mul};

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lager::Matrix;

fn add<T: Copy + Default + Mul<Output = T> + Add<Output = T>, const M: usize, const N: usize>(
    a: Matrix<T, M, N>,
    b: Matrix<T, M, N>,
) -> Matrix<T, M, N> {
    a + b
}

fn mul<T: Copy + Default + Mul<Output = T>, const M: usize, const N: usize>(
    a: Matrix<T, M, N>,
    b: Matrix<T, M, N>,
) -> Matrix<T, M, N> {
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
