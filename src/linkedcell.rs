use boxarray::boxarray_;

use crate::{boundary::Boundary, soacontainerdyn::SoAContainerDyn, soavector::SoAVectorDyn};

/// Simple implementation of a linked cell utilizing SoAcontainer under the hood.
/// It has a fixed number of cells in each dimension, and a cutoff distance, which is used to determine the size of the cells.
pub struct LinkedCell<T, const D0: usize, const D1: usize, const D2: usize> {
    pub cutoff: T,
    pub cells: Box<[[[SoAContainerDyn<T>; D2]; D1]; D0]>, // Hell yeah, 3D array of const generics!
    pub boundaries: Vec<crate::boundary::BoundaryCondition<T>>,
}

impl<T, const D0: usize, const D1: usize, const D2: usize> LinkedCell<T, D0, D1, D2>
where
    T: crate::force::LjFloat,
{
    // This implementation is not possible because SoAContainerDyn is not Copy.
    // While I did try to implement a fast implementation, I was limited by the way multidimensional arrays work with uninitialized values.
    // /// Initialize an empty LinkedCell with the given cutoff distance.
    // pub fn new(cutoff: T) -> Self {
    //     // This doesn't work because SoAContainerDyn is not Copy
    //     // let mut cells = [[[SoAContainerDyn::new(0); D2]; D1]; D0];
    //     // Instead, we can do some uninit magic.

    //     #[allow(clippy::uninit_assumed_init)]
    //     let mut cells: [[[SoAContainerDyn<T>; D2]; D1]; D0] =
    //         unsafe { std::mem::MaybeUninit::uninit().assume_init() };

    //     // Note that the following also requires the Copy trait, which SoAContainerDyn doesn't have.
    //     // let mut cells: [[[SoAContainerDyn<T>; D2]; D1]; D0] = unsafe { [[[std::mem::MaybeUninit::uninit().assume_init(); D2]; D1]; D0] };

    //     for i in 0..D0 {
    //         for j in 0..D1 {
    //             for k in 0..D2 {
    //                 cells[i][j][k] = SoAContainerDyn::new(0);
    //             }
    //         }
    //     }
    //     LinkedCell {
    //         cutoff,
    //         cells,
    //         boundaries: Vec::new(),
    //     }
    // }

    /// Initialize an empty LinkedCell with the given cutoff distance.
    pub fn new(cutoff: T) -> Self {
        let constructor = |_| SoAContainerDyn::<T>::new(0);
        let cells: Box<[[[SoAContainerDyn<T>; D2]; D1]; D0]> = boxarray_(constructor); // This is a bit more readable, but still not perfect. Also pretty cursed.
        LinkedCell {
            cutoff,
            cells,
            boundaries: Vec::new(),
        }
    }

    /// Naively redistributes the particles in the LinkedCell to the correct cells.
    /// First collects all particles in the cells, then puts them back in the correct cells.
    pub fn redistribute_particles_slow(&mut self) {
        let mut positions: SoAVectorDyn<T> = SoAVectorDyn::new(0);
        let mut velocities: SoAVectorDyn<T> = SoAVectorDyn::new(0);
        let mut old_forces: SoAVectorDyn<T> = SoAVectorDyn::new(0);
        let mut current_forces: SoAVectorDyn<T> = SoAVectorDyn::new(0);

        // Collect all particles in the cells
        for i in 0..D0 {
            for j in 0..D1 {
                for k in 0..D2 {
                    positions.x.extend(self.cells[i][j][k].position.x.iter());
                    positions.y.extend(self.cells[i][j][k].position.y.iter());
                    positions.z.extend(self.cells[i][j][k].position.z.iter());

                    velocities.x.extend(self.cells[i][j][k].velocity.x.iter());
                    velocities.y.extend(self.cells[i][j][k].velocity.y.iter());
                    velocities.z.extend(self.cells[i][j][k].velocity.z.iter());

                    old_forces.x.extend(self.cells[i][j][k].old_force.x.iter());
                    old_forces.y.extend(self.cells[i][j][k].old_force.y.iter());
                    old_forces.z.extend(self.cells[i][j][k].old_force.z.iter());

                    current_forces
                        .x
                        .extend(self.cells[i][j][k].current_force.x.iter());
                    current_forces
                        .y
                        .extend(self.cells[i][j][k].current_force.y.iter());
                    current_forces
                        .z
                        .extend(self.cells[i][j][k].current_force.z.iter());
                }
            }
        }

        // Redistribute the particles to the correct cells
        // Note that the data doesn't need to be cleared, as the content is overwritten.
        for i in 0..D0 {
            let x_min = T::from_usize(i).unwrap() * self.cutoff;
            let x_max = T::from_usize(i + 1).unwrap() * self.cutoff;
            for j in 0..D1 {
                let y_min = T::from_usize(j).unwrap() * self.cutoff;
                let y_max = T::from_usize(j + 1).unwrap() * self.cutoff;
                for k in 0..D2 {
                    let z_min = T::from_usize(k).unwrap() * self.cutoff;
                    let z_max = T::from_usize(k + 1).unwrap() * self.cutoff;

                    // Now we need to filter the particles that belong to this cell
                    let mut new_positions = SoAVectorDyn::new(0);
                    let mut new_velocities = SoAVectorDyn::new(0);
                    let mut new_old_forces = SoAVectorDyn::new(0);
                    let mut new_current_forces = SoAVectorDyn::new(0);
                    for l in 0..positions.x.len() {
                        if positions.x[l] >= x_min
                            && positions.x[l] < x_max
                            && positions.y[l] >= y_min
                            && positions.y[l] < y_max
                            && positions.z[l] >= z_min
                            && positions.z[l] < z_max
                        {
                            new_positions.x.push(positions.x[l]);
                            new_positions.y.push(positions.y[l]);
                            new_positions.z.push(positions.z[l]);

                            new_velocities.x.push(velocities.x[l]);
                            new_velocities.y.push(velocities.y[l]);
                            new_velocities.z.push(velocities.z[l]);

                            new_old_forces.x.push(old_forces.x[l]);
                            new_old_forces.y.push(old_forces.y[l]);
                            new_old_forces.z.push(old_forces.z[l]);

                            new_current_forces.x.push(current_forces.x[l]);
                            new_current_forces.y.push(current_forces.y[l]);
                            new_current_forces.z.push(current_forces.z[l]);
                        }
                    }

                    self.cells[i][j][k].position = new_positions;
                    self.cells[i][j][k].velocity = new_velocities;
                    self.cells[i][j][k].old_force = new_old_forces;
                    self.cells[i][j][k].current_force = new_current_forces;
                    // Note that because the content is overwritten, particles can be lost if they are outside the linked cell.
                }
            }
        }
    }

    /// Redistributes the particles in the LinkedCell to the correct cells.
    /// First collects all particles in the cells that don't belong, then puts them back in the correct cells.
    pub fn redistribute_particles(&mut self) {
        let mut positions: SoAVectorDyn<T> = SoAVectorDyn::new(0);
        let mut velocities: SoAVectorDyn<T> = SoAVectorDyn::new(0);
        let mut old_forces: SoAVectorDyn<T> = SoAVectorDyn::new(0);
        let mut current_forces: SoAVectorDyn<T> = SoAVectorDyn::new(0);

        /// Helper function that takes in a SoAVectorDyn and an index and deletes the element at that index.
        /// It is then copied to the end of the second vector
        fn move_to_local<T>(
            vec: &mut SoAVectorDyn<T>,
            index: usize,
            local_vec: &mut SoAVectorDyn<T>,
        ) {
            // It has am x,y,z component, so we need to remove all three components
            let x = vec.x.swap_remove(index);
            let y = vec.y.swap_remove(index);
            let z = vec.z.swap_remove(index);

            local_vec.x.push(x);
            local_vec.y.push(y);
            local_vec.z.push(z);
        }

        for i in 0..D0 {
            let x_min = T::from_usize(i).unwrap() * self.cutoff;
            let x_max = T::from_usize(i + 1).unwrap() * self.cutoff;
            for j in 0..D1 {
                let y_min = T::from_usize(j).unwrap() * self.cutoff;
                let y_max = T::from_usize(j + 1).unwrap() * self.cutoff;
                for k in 0..D2 {
                    let z_min = T::from_usize(k).unwrap() * self.cutoff;
                    let z_max = T::from_usize(k + 1).unwrap() * self.cutoff;

                    // We know which cell we have and where the particle should be in space.
                    // Iterate over the particles in the cell and check if they belong to the cell.

                    let mut index = 0;
                    while index < self.cells[i][j][k].position.x.len() {
                        if self.cells[i][j][k].position.x[index] < x_min
                            || self.cells[i][j][k].position.x[index] >= x_max
                            || self.cells[i][j][k].position.y[index] < y_min
                            || self.cells[i][j][k].position.y[index] >= y_max
                            || self.cells[i][j][k].position.z[index] < z_min
                            || self.cells[i][j][k].position.z[index] >= z_max
                        {
                            // The particle doesn't belong to this cell
                            // Move it to the local vector
                            move_to_local(&mut self.cells[i][j][k].position, index, &mut positions);
                            move_to_local(
                                &mut self.cells[i][j][k].velocity,
                                index,
                                &mut velocities,
                            );
                            move_to_local(
                                &mut self.cells[i][j][k].old_force,
                                index,
                                &mut old_forces,
                            );
                            move_to_local(
                                &mut self.cells[i][j][k].current_force,
                                index,
                                &mut current_forces,
                            );
                        } else {
                            index += 1;
                        }
                    }
                }
            }
        }

        // Redistribute the particles to the correct cells
        // Note that the data doesn't need to be cleared, as the content is overwritten.
        for i in 0..D0 {
            let x_min = T::from_usize(i).unwrap() * self.cutoff;
            let x_max = T::from_usize(i + 1).unwrap() * self.cutoff;
            for j in 0..D1 {
                let y_min = T::from_usize(j).unwrap() * self.cutoff;
                let y_max = T::from_usize(j + 1).unwrap() * self.cutoff;
                for k in 0..D2 {
                    let z_min = T::from_usize(k).unwrap() * self.cutoff;
                    let z_max = T::from_usize(k + 1).unwrap() * self.cutoff;

                    // Now we need to filter the particles that belong to this cell
                    let mut new_positions = SoAVectorDyn::new(0);
                    let mut new_velocities = SoAVectorDyn::new(0);
                    let mut new_old_forces = SoAVectorDyn::new(0);
                    let mut new_current_forces = SoAVectorDyn::new(0);
                    for l in 0..positions.x.len() {
                        if positions.x[l] >= x_min
                            && positions.x[l] < x_max
                            && positions.y[l] >= y_min
                            && positions.y[l] < y_max
                            && positions.z[l] >= z_min
                            && positions.z[l] < z_max
                        {
                            new_positions.x.push(positions.x[l]);
                            new_positions.y.push(positions.y[l]);
                            new_positions.z.push(positions.z[l]);

                            new_velocities.x.push(velocities.x[l]);
                            new_velocities.y.push(velocities.y[l]);
                            new_velocities.z.push(velocities.z[l]);

                            new_old_forces.x.push(old_forces.x[l]);
                            new_old_forces.y.push(old_forces.y[l]);
                            new_old_forces.z.push(old_forces.z[l]);

                            new_current_forces.x.push(current_forces.x[l]);
                            new_current_forces.y.push(current_forces.y[l]);
                            new_current_forces.z.push(current_forces.z[l]);
                        }
                    }

                    self.cells[i][j][k]
                        .position
                        .x
                        .extend(new_positions.x.iter());
                    self.cells[i][j][k]
                        .position
                        .y
                        .extend(new_positions.y.iter());
                    self.cells[i][j][k]
                        .position
                        .z
                        .extend(new_positions.z.iter());

                    self.cells[i][j][k]
                        .velocity
                        .x
                        .extend(new_velocities.x.iter());
                    self.cells[i][j][k]
                        .velocity
                        .y
                        .extend(new_velocities.y.iter());
                    self.cells[i][j][k]
                        .velocity
                        .z
                        .extend(new_velocities.z.iter());

                    self.cells[i][j][k]
                        .old_force
                        .x
                        .extend(new_old_forces.x.iter());
                    self.cells[i][j][k]
                        .old_force
                        .y
                        .extend(new_old_forces.y.iter());
                    self.cells[i][j][k]
                        .old_force
                        .z
                        .extend(new_old_forces.z.iter());

                    self.cells[i][j][k]
                        .current_force
                        .x
                        .extend(new_current_forces.x.iter());
                    self.cells[i][j][k]
                        .current_force
                        .y
                        .extend(new_current_forces.y.iter());
                    self.cells[i][j][k]
                        .current_force
                        .z
                        .extend(new_current_forces.z.iter());
                }
            }
        }

        // Note that this implementation is not very efficient, as it requires a lot of copying.
        // A more efficient implementation would be to keep track of the indices of the particles in the cells.
    }

    /// Adds particles to the LinkedCell in a grid of given size.
    /// Also redistributes the particles to the correct cells.
    pub fn add_particle_grid(&mut self, grid_size: usize, grid_distance: T) {
        // The different min and max values for the grid
        let x_min = T::zero() + (self.cutoff / T::two());
        let x_max = T::from(grid_size).unwrap() * grid_distance + (self.cutoff / T::two());

        let y_min = T::zero() + (self.cutoff / T::two());
        let y_max = T::from(grid_size).unwrap() * grid_distance + (self.cutoff / T::two());

        let z_min = T::zero() + (self.cutoff / T::two());
        let z_max = T::from(grid_size).unwrap() * grid_distance + (self.cutoff / T::two());

        for i in 0..grid_size {
            for j in 0..grid_size {
                for k in 0..grid_size {
                    let x = x_min
                        + (x_max - x_min) * T::from_usize(i).unwrap()
                            / T::from_usize(grid_size).unwrap();
                    let y = y_min
                        + (y_max - y_min) * T::from_usize(j).unwrap()
                            / T::from_usize(grid_size).unwrap();
                    let z = z_min
                        + (z_max - z_min) * T::from_usize(k).unwrap()
                            / T::from_usize(grid_size).unwrap();

                    let cell_x = (x / self.cutoff).to_usize().unwrap().min(D0 - 1);
                    let cell_y = (y / self.cutoff).to_usize().unwrap().min(D1 - 1);
                    let cell_z = (z / self.cutoff).to_usize().unwrap().min(D2 - 1);

                    self.cells[cell_x][cell_y][cell_z].position.x.push(x);
                    self.cells[cell_x][cell_y][cell_z].position.y.push(y);
                    self.cells[cell_x][cell_y][cell_z].position.z.push(z);

                    self.cells[cell_x][cell_y][cell_z]
                        .velocity
                        .x
                        .push(T::zero());
                    self.cells[cell_x][cell_y][cell_z]
                        .velocity
                        .y
                        .push(T::zero());
                    self.cells[cell_x][cell_y][cell_z]
                        .velocity
                        .z
                        .push(T::zero());

                    self.cells[cell_x][cell_y][cell_z]
                        .old_force
                        .x
                        .push(T::zero());
                    self.cells[cell_x][cell_y][cell_z]
                        .old_force
                        .y
                        .push(T::zero());
                    self.cells[cell_x][cell_y][cell_z]
                        .old_force
                        .z
                        .push(T::zero());

                    self.cells[cell_x][cell_y][cell_z]
                        .current_force
                        .x
                        .push(T::zero());
                    self.cells[cell_x][cell_y][cell_z]
                        .current_force
                        .y
                        .push(T::zero());
                    self.cells[cell_x][cell_y][cell_z]
                        .current_force
                        .z
                        .push(T::zero());
                }
            }
        }

        self.redistribute_particles();
    }

    /// Flushes the forces of all particles in the LinkedCell.
    pub fn flush_forces(&mut self) {
        for i in 0..D0 {
            for j in 0..D1 {
                for k in 0..D2 {
                    self.cells[i][j][k].flush_forces();
                }
            }
        }
    }

    /// Saves the current state of the LinkedCell to a CSV file.
    /// The file will be named 'output_{index}.csv'
    /// in the 'results' directory.
    ///
    /// Note: this function is quite slow, but I couldn't find a way to make it faster.
    /// Even manually writing the data to a file was just as slow.
    pub fn save_to_csv(&self, index: usize)
    where
        T: std::string::ToString,
    {
        let dir = std::path::Path::new("results");
        std::fs::create_dir_all(dir).unwrap();
        let path = dir.join(format!("output_{}.csv", index));

        let mut writer = csv::Writer::from_path(path).unwrap(); // It might look like this writer is slow, but I tried to implement it manually by writing the raw data to a file, and it was just as slow.
                                                                // Format: first line "x,y,z,v" then N lines of data
                                                                // (v is the magnitude of the velocity)
        writer.write_record(["x", "y", "z", "v"]).unwrap();
        for i in 0..D0 {
            for j in 0..D1 {
                for k in 0..D2 {
                    for l in 0..self.cells[i][j][k].position.x.len() {
                        let vx = self.cells[i][j][k].velocity.x[l];
                        let vy = self.cells[i][j][k].velocity.y[l];
                        let vz = self.cells[i][j][k].velocity.z[l];
                        let v = (vx * vx + vy * vy + vz * vz).sqrt();
                        writer
                            .write_record(&[
                                self.cells[i][j][k].position.x[l].to_string(),
                                self.cells[i][j][k].position.y[l].to_string(),
                                self.cells[i][j][k].position.z[l].to_string(),
                                v.to_string(),
                            ])
                            .unwrap();
                    }
                }
            }
        }
    }

    pub fn apply_boundaries(&mut self) {
        for boundary in &self.boundaries {
            for i in 0..D0 {
                for j in 0..D1 {
                    for k in 0..D2 {
                        self.cells[i][j][k].apply_boundary(boundary);
                    }
                }
            }
        }
    }
}
