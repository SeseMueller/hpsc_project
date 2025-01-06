pub struct SoAVector<T, const N: usize> {
    pub x: [T; N],
    pub y: [T; N],
    pub z: [T; N],
}