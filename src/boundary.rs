/// The dimension of the boundary condition; Either X, Y, or Z.
pub enum Dimension {
    X,
    Y,
    Z,
}

/// Represents a boundary condition that can be applied to a container.
pub struct BoundaryCondition<T> {
    pub dimension: Dimension,
    pub min: T,
    pub max: T,
}

/// Represents a boundary condition that can be applied to a container.
/// Only requires the container to be able to run the "apply_boundary" method.
/// (We only use reflective boundary conditions in this project)
pub trait Boundary<T> {
    /// Formulas:
    /// In position:
    /// - if `xi < min`, `xi = 2 * min - xi`
    /// - if `xi > max`, `xi = 2 * max - xi`
    ///
    /// In velocity:
    /// - if `xi < min`, `xi = -xi`
    /// - if `xi > max`, `xi = -xi`
    fn apply_boundary(&mut self, boundary: &BoundaryCondition<T>);
}
