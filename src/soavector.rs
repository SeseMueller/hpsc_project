use crate::force::LjFloat;

#[derive(Debug, Clone, Copy)]
pub struct SoAVector<T, const N: usize> {
    pub x: [T; N],
    pub y: [T; N],
    pub z: [T; N],
}

impl<T, const N: usize> SoAVector<T, N>
where
    T: LjFloat,
{
    pub fn new() -> Self {
        SoAVector {
            x: [T::zero(); N],
            y: [T::zero(); N],
            z: [T::zero(); N],
        }
    }
}

#[derive(Debug, Clone)]
pub struct SoAVectorDyn<T> {
    pub x: Vec<T>,
    pub y: Vec<T>,
    pub z: Vec<T>,
}

impl<T> SoAVectorDyn<T>
where
    T: LjFloat,
{
    pub fn new(num_particles: usize) -> Self {
        SoAVectorDyn {
            x: vec![T::zero(); num_particles],
            y: vec![T::zero(); num_particles],
            z: vec![T::zero(); num_particles],
        }
    }
}
