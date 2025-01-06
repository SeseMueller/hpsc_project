use crate::{force::LjFloat, soavector};

pub struct SoAContainer<T, const N: usize> {
    pub position: soavector::SoAVector<T, N>,
    pub velocity: soavector::SoAVector<T, N>,
    pub old_force: soavector::SoAVector<T, N>,
    pub current_force: soavector::SoAVector<T, N>,
}

impl<T, const N: usize> SoAContainer<T, N> where T: LjFloat + std::string::ToString {
    /// Initialize the SoAContainer with zero forces and velocities 
    /// at random positions within the 10 cube.
    pub fn init_random_pos() -> Self {
        let mut container = SoAContainer {
            position: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
            velocity: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
            old_force: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
            current_force: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
        };
        for i in 0..N {
            container.position.x[i] = T::from_f64(rand::random::<f64>() * 10.0).unwrap();
            container.position.y[i] = T::from_f64(rand::random::<f64>() * 10.0).unwrap();
            container.position.z[i] = T::from_f64(rand::random::<f64>() * 10.0).unwrap();
        }
        container
    }

    /// Distributes the particles in a 3D lattice.
    pub fn init_3d_lattice() -> Self {
        // The size of the lattice such that the particles are evenly distributed
        let lattice_size = N as f64;
        let lattice_size = lattice_size.cbrt().ceil() as usize;

        let mut container = SoAContainer {
            position: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
            velocity: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
            old_force: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
            current_force: soavector::SoAVector {
                x: [T::zero(); N],
                y: [T::zero(); N],
                z: [T::zero(); N],
            },
        };

        let mut i = 0;
        for x in 0..lattice_size {
            for y in 0..lattice_size {
                for z in 0..lattice_size {
                    if i >= N {
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

    /// Save the current state of the container to a CSV file.
    /// It will be saved to a file named `output_{index}.csv`
    /// in the `results` directory.
    pub fn save_to_csv(index: usize, container: &SoAContainer<T, N>) {
        let dir = std::path::Path::new("results");
        std::fs::create_dir_all(dir).unwrap();
        let path = dir.join(format!("output.{}.csv", index));
        let mut wtr = csv::Writer::from_path(path).unwrap();
        // Format: first line "x,y,z,v" then N lines of data
        // (v is the magnitude of the velocity)
        wtr.write_record(["x", "y", "z", "v"]).unwrap();
        for i in 0..N {
            let vx = container.velocity.x[i];
            let vy = container.velocity.y[i];
            let vz = container.velocity.z[i];
            let v = (vx*vx + vy*vy + vz*vz).sqrt();
            wtr.write_record(&[
                container.position.x[i].to_string(),
                container.position.y[i].to_string(),
                container.position.z[i].to_string(),
                v.to_string(),
            ]).unwrap();
        }
    }

    /// Swap the old force with the current force and clear the (now old) force.
    pub fn flush_forces(&mut self) {
        std::mem::swap(&mut self.old_force, &mut self.current_force); // This interacts better with the borrow checker
        self.old_force = soavector::SoAVector {
            x: [T::zero(); N],
            y: [T::zero(); N],
            z: [T::zero(); N],
        };
    }
}
