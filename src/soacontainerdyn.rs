use crate::{
    boundary::{self, Boundary, BoundaryCondition},
    force::LjFloat,
    soavector::SoAVectorDyn,
};

/// An SoA container that can dynamically grow and shrink.
/// Functionally equivalent to `SoAContainer`, but with a different implementation
/// allowing for dynamic resizing, useful for linked cell implementations.
/// Note that because of the lack of a proper effect system in Rust,
/// it is not possible to declare a trait over the static and dynamic SoA containers
/// that would allow for the same code to be used for both while still having the
/// performance benefits of the static version.
/// A way to do it would be to put the accesses to the field behind a function call
/// which could then be manipulated and read through a trait over both SoAVectors.
/// However, Static dispatch could then not optimize the code as well as it could with
/// direct field access.
#[derive(Debug, Clone)]
pub struct SoAContainerDyn<T> {
    pub position: SoAVectorDyn<T>,
    pub velocity: SoAVectorDyn<T>,
    pub old_force: SoAVectorDyn<T>,
    pub current_force: SoAVectorDyn<T>,
}

impl<T> SoAContainerDyn<T>
where
    T: LjFloat,
{
    /// Creates a new `SoAContainerDyn` with the given number of particles.
    pub fn new(num_particles: usize) -> Self {
        let mut container = SoAContainerDyn {
            position: SoAVectorDyn::new(num_particles),
            velocity: SoAVectorDyn::new(num_particles),
            old_force: SoAVectorDyn::new(num_particles),
            current_force: SoAVectorDyn::new(num_particles),
        };
        for i in 0..num_particles {
            container.position.x[i] = T::zero();
            container.position.y[i] = T::zero();
            container.position.z[i] = T::zero();
        }
        container
    }

    /// Distributes the particles in a 3D lattice.
    #[allow(dead_code)]
    pub fn init_3d_lattice(num_particles: usize) -> Self {
        // The size of the lattice such that the particles are evenly distributed
        let lattice_size = num_particles as f64;
        let lattice_size = lattice_size.cbrt().ceil() as usize;

        let mut container = SoAContainerDyn {
            position: SoAVectorDyn::new(num_particles),
            velocity: SoAVectorDyn::new(num_particles),
            old_force: SoAVectorDyn::new(num_particles),
            current_force: SoAVectorDyn::new(num_particles),
        };

        let mut i = 0;
        for x in 0..lattice_size {
            for y in 0..lattice_size {
                for z in 0..lattice_size {
                    if i >= num_particles {
                        break;
                    }
                    container.position.x[i] = T::from_f64(x as f64).unwrap();
                    container.position.y[i] = T::from_f64(y as f64).unwrap();
                    container.position.z[i] = T::from_f64(z as f64).unwrap();
                    i += 1;
                }
            }
        }
        container
    }

    /// Swaps the old forces with the current forces and clear the (now old) force.
    pub fn flush_forces(&mut self) {
        std::mem::swap(&mut self.old_force, &mut self.current_force);
        self.current_force = SoAVectorDyn::new(self.current_force.x.len());
    }
}

impl<T> Boundary<T> for SoAContainerDyn<T>
where
    T: LjFloat,
{
    fn apply_boundary(&mut self, boundary: &BoundaryCondition<T>) {
        let (velocities, positions) = match boundary.dimension {
            boundary::Dimension::X => (&mut self.velocity.x, &mut self.position.x),
            boundary::Dimension::Y => (&mut self.velocity.y, &mut self.position.y),
            boundary::Dimension::Z => (&mut self.velocity.z, &mut self.position.z),
        };

        for i in 0..positions.len() {
            if positions[i] < boundary.min {
                // For double tunneling or more, we take the modulus of the position
                if positions[i] < boundary.min - boundary.max {
                    positions[i] = boundary.max + (positions[i] % boundary.max);
                // Don't reflect, just move to the other side
                } else {
                    positions[i] = boundary.min + (boundary.min - positions[i]); // Reflect: x = 2 * min - x
                    velocities[i] = -velocities[i]; // Reflect: v = -v
                }
            } else if positions[i] > boundary.max {
                // For double tunneling or more, we take the modulus of the position
                if positions[i] > boundary.max + boundary.max {
                    positions[i] = boundary.min + (positions[i] % boundary.max);
                // Don't reflect, just move to the other side
                } else {
                    positions[i] = boundary.max + (boundary.max - positions[i]); // Reflect: x = 2 * max - x
                    velocities[i] = -velocities[i]; // Reflect: v = -v
                }
            }
        }
    }
}
