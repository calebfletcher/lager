use std::ops::Neg;

pub trait Abs {
    fn abs(&self) -> Self
    where
        Self: Sized + PartialOrd + Neg<Output = Self> + From<f32> + Copy,
    {
        if *self >= 0.0.into() {
            *self
        } else {
            -*self
        }
    }
}

impl Abs for f32 {}
impl Abs for f64 {}
impl Abs for u8 {}
impl Abs for i8 {}
impl Abs for u16 {}
impl Abs for i16 {}
impl Abs for u32 {}
impl Abs for i32 {}
impl Abs for u64 {}
impl Abs for i64 {}
