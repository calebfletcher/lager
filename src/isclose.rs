use std::fmt::Debug;
use std::ops::{Add, Mul, Neg, Sub};

use crate::abs::Abs;

pub trait IsClose {
    fn isclose(&self, other: Self) -> bool
    where
        Self: Sized
            + From<f64>
            + Abs
            + PartialOrd
            + Copy
            + Debug
            + Neg<Output = Self>
            + Add<Output = Self>
            + Sub<Output = Self>
            + Mul<Output = Self>,
    {
        let rt: Self = 1e-05.into();
        let at: Self = 1e-08.into();

        let threshold = rt * self.abs() + at;
        (*self - other).abs() < threshold
    }
}

impl IsClose for f32 {}
impl IsClose for f64 {}
