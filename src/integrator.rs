use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::{force::LjFloat, linkedcell::LinkedCell, soacontainer::SoAContainer, soacontainerdyn::SoAContainerDyn};

/// A Velocity Verlet integrator for all containers of this project.
pub struct VelocityStörmerVerlet<T> {
    pub dt: T,
    pub half_dt: T,
    pub half_dt_sq: T,
}

impl<T> VelocityStörmerVerlet<T> where T: LjFloat {
    pub fn new(dt: T) -> Self {
        Self {
            dt,
            half_dt: dt / (T::one() + T::one()),
            half_dt_sq: dt * dt / (T::one() + T::one()),
        }
    }
}
impl<T> VelocityStörmerVerlet<T> where T: LjFloat + std::ops::AddAssign {

    /// Update the position of the particles in the container for a fixed number of particles.
    #[inline(never)] // To be able to see the function in the profiler
    #[allow(dead_code)]
    pub fn update_position< const N: usize> (&self, container: &mut SoAContainer<T, N>) {
        for i in 0..N {
            container.position.x[i] += container.velocity.x[i] * self.dt + container.current_force.x[i] * self.half_dt_sq;
            container.position.y[i] += container.velocity.y[i] * self.dt + container.current_force.y[i] * self.half_dt_sq;
            container.position.z[i] += container.velocity.z[i] * self.dt + container.current_force.z[i] * self.half_dt_sq;
        }
    }
    /// Note that neither update_position nor update_velocity can be trivially parallelized because of the ownership structure of the SoAContainer.
    #[allow(dead_code)]
    #[inline(never)] // To be able to see the function in the profiler
    pub fn update_velocity< const N: usize> (&self, container: &mut SoAContainer<T, N>) {
        for i in 0..N {
            container.velocity.x[i] += (container.current_force.x[i] + container.old_force.x[i]) * self.half_dt;
            container.velocity.y[i] += (container.current_force.y[i] + container.old_force.y[i]) * self.half_dt;
            container.velocity.z[i] += (container.current_force.z[i] + container.old_force.z[i]) * self.half_dt;
        }
    }

    /// Update the position of the particles in the container for a dynamic number of particles.
    #[allow(dead_code)]
    #[inline(never)] // To be able to see the function in the profiler
    pub fn update_position_dyn (&self, container: &mut SoAContainerDyn<T>) {
        for i in 0..container.position.x.len() {
            container.position.x[i] += container.velocity.x[i] * self.dt + container.current_force.x[i] * self.half_dt_sq;
            container.position.y[i] += container.velocity.y[i] * self.dt + container.current_force.y[i] * self.half_dt_sq;
            container.position.z[i] += container.velocity.z[i] * self.dt + container.current_force.z[i] * self.half_dt_sq;
        }
    }

    /// Update the velocity of the particles in the container for a dynamic number of particles.
    #[inline(never)] // To be able to see the function in the profiler
    #[allow(dead_code)]
    pub fn update_velocity_dyn (&self, container: &mut SoAContainerDyn<T>) {
        for i in 0..container.position.x.len() {
            container.velocity.x[i] += (container.current_force.x[i] + container.old_force.x[i]) * self.half_dt;
            container.velocity.y[i] += (container.current_force.y[i] + container.old_force.y[i]) * self.half_dt;
            container.velocity.z[i] += (container.current_force.z[i] + container.old_force.z[i]) * self.half_dt;
        }
    }


    /// Update the position of all particles in a linked cell structure.
    #[allow(dead_code)]
    pub fn update_position_linked<const D0: usize, const D1: usize, const D2: usize>(&self, linked_cell: &mut LinkedCell<T, D0, D1, D2>) {
        for i in 0..D0 {
            for j in 0..D1 {
                for k in 0..D2 {
                    for l in 0..linked_cell.cells[i][j][k].position.x.len() {
                        linked_cell.cells[i][j][k].position.x[l] += linked_cell.cells[i][j][k].velocity.x[l] * self.dt + linked_cell.cells[i][j][k].current_force.x[l] * self.half_dt_sq;
                        linked_cell.cells[i][j][k].position.y[l] += linked_cell.cells[i][j][k].velocity.y[l] * self.dt + linked_cell.cells[i][j][k].current_force.y[l] * self.half_dt_sq;
                        linked_cell.cells[i][j][k].position.z[l] += linked_cell.cells[i][j][k].velocity.z[l] * self.dt + linked_cell.cells[i][j][k].current_force.z[l] * self.half_dt_sq;
                    }
                }
            }
        }
    }

    /// Update the position of all particles in a linked cell structure in a parallel manner.
    pub fn update_position_linked_par<const D0: usize, const D1: usize, const D2: usize>(&self, linked_cell: &mut LinkedCell<T, D0, D1, D2>) 
    where T: Send + Sync
    {
        linked_cell.cells.par_iter_mut().for_each(|plane| {
            for j in 0..D1 {
                for k in 0..D2 {
                    for l in 0..plane[j][k].position.x.len() {
                        plane[j][k].position.x[l] += plane[j][k].velocity.x[l] * self.dt + plane[j][k].current_force.x[l] * self.half_dt_sq;
                        plane[j][k].position.y[l] += plane[j][k].velocity.y[l] * self.dt + plane[j][k].current_force.y[l] * self.half_dt_sq;
                        plane[j][k].position.z[l] += plane[j][k].velocity.z[l] * self.dt + plane[j][k].current_force.z[l] * self.half_dt_sq;
                    }
                }
            }
        });
    }

    /// Update the velocity of all particles in a linked cell structure.
    #[allow(dead_code)]
    pub fn update_velocity_linked<const D0: usize, const D1: usize, const D2: usize>(&self, linked_cell: &mut LinkedCell<T, D0, D1, D2>) {
        for i in 0..D0 {
            for j in 0..D1 {
                for k in 0..D2 {
                    for l in 0..linked_cell.cells[i][j][k].position.x.len() {
                        linked_cell.cells[i][j][k].velocity.x[l] += (linked_cell.cells[i][j][k].current_force.x[l] + linked_cell.cells[i][j][k].old_force.x[l]) * self.half_dt;
                        linked_cell.cells[i][j][k].velocity.y[l] += (linked_cell.cells[i][j][k].current_force.y[l] + linked_cell.cells[i][j][k].old_force.y[l]) * self.half_dt;
                        linked_cell.cells[i][j][k].velocity.z[l] += (linked_cell.cells[i][j][k].current_force.z[l] + linked_cell.cells[i][j][k].old_force.z[l]) * self.half_dt;
                    }
                }
            }
        }
    }

    /// Update the velocity of all particles in a linked cell structure in a parallel manner.
    pub fn update_velocity_linked_par<const D0: usize, const D1: usize, const D2: usize>(&self, linked_cell: &mut LinkedCell<T, D0, D1, D2>) 
    where T: Send + Sync
    {
        linked_cell.cells.par_iter_mut().for_each(|plane| {
            for j in 0..D1 {
                for k in 0..D2 {
                    for l in 0..plane[j][k].position.x.len() {
                        plane[j][k].velocity.x[l] += (plane[j][k].current_force.x[l] + plane[j][k].old_force.x[l]) * self.half_dt;
                        plane[j][k].velocity.y[l] += (plane[j][k].current_force.y[l] + plane[j][k].old_force.y[l]) * self.half_dt;
                        plane[j][k].velocity.z[l] += (plane[j][k].current_force.z[l] + plane[j][k].old_force.z[l]) * self.half_dt;
                    }
                }
            }
        });
    }
}