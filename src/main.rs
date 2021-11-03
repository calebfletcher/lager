#![feature(bench_black_box)]
use lager::Matrix;
use std::hint::black_box;

fn main() {
    let mtx1 = black_box(Matrix::new([
        [4.4, 7.3, 2.5],
        [2.4, 6.4, 2.2],
        [9.5, 3.5, 4.6],
    ]));
    let mtx2 = black_box(Matrix::new([
        [4.4, 7.3, 2.5],
        [2.4, 6.4, 2.2],
        [9.5, 3.5, 4.6],
    ]));
    let _mtx3 = black_box(Matrix::new([[4.4, 7.3, 2.5], [2.4, 6.4, 2.2]]));
    let mtx4 = mtx1 + mtx2;
    dbg!("{}", mtx4.values);
}
