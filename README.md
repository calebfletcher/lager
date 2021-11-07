# lager
A lightweight linear algebra implementation in Rust.

## Purpose

Lager aims to be a lightweight implementation of various linear algebra techniques which is optimised for small matrices and allows a greater level of compile time checking of matrix math.

## Advantages
 - Compile time checking of matrix dimensions for operations like matrix multiplication and addition.
 - All matrices stored on the stack to increase performance. Note that due to this, if you want to use very large arrays you may want to look at other crates which store matrices on the heap like ndarray.
 - Straightforward numpy-inspired API.

## Disadvantages
 - Uses some nightly features currently, and so does not compile on stable:
   - `generic_const_exprs`
 - Does not have a huge number of features or calls to BLAS/LAPACK implementations