mod soavector;

mod soacontainer;

use std::io::Write;

use boundary::BoundaryCondition;
use soacontainer::SoAContainer;

mod integrator;

mod force;

mod boundary;

mod linkedcell;

mod soacontainerdyn;

fn main() {
    // run_soa_simulation();
    
    run_linked_cell_simulation();
}

// Controls the number of cells ("containers") in each dimension, the cutoff distance, and the number of particles in each dimension

#[allow(dead_code)]
const TINY_SIM: (usize, f64, usize) = (10, 1.5, 15); // Example: 10^3 containers, each 1.5^3 units in size, with 15^3 particles
#[allow(dead_code)]
const SMALL_SIM: (usize, f64, usize) = (20, 1.5, 30);
#[allow(dead_code)]
const STANDARD_SIM: (usize, f64, usize) = (20, 1.5, 30);
#[allow(dead_code)]
const BIG_SIM: (usize, f64, usize) = (30, 1.5, 40);
#[allow(dead_code)]
const BIG_SLOW_SIM: (usize, f64, usize) = (40, 1.0, 40); // This is possible, but interestingly actually slower. 
#[allow(dead_code)]
const HUGE_SIM: (usize, f64, usize) = (40, 1.5, 60);

const SELECTED_SIM: (usize, f64, usize) = BIG_SIM;

/// Runs a simple simulation using the SoA container.
/// Saves the results to CSV files for visualization.
#[allow(dead_code)]
fn run_soa_simulation() {
    const LOG_TIMESTEPS: usize = 50;
    const ITERATIONS: usize = 5000;

    let mut total_simulated_time = 0.0;
    let mut last_time = std::time::Instant::now();

    let integrator = integrator::VelocityStörmerVerlet::new(0.001);
    let mut cont = SoAContainer::<f64, 1000>::init_3d_lattice();
    force::apply_lj_force_soa(&mut cont); // Calculate initial forces

    let mut log_step = 0;
    for i in 0..ITERATIONS {
        if i % LOG_TIMESTEPS == 0 {
            total_simulated_time += last_time.elapsed().as_secs_f64();

            soacontainer::SoAContainer::save_to_csv(log_step, &cont);
            log_step += 1;

            last_time = std::time::Instant::now();
        }
        integrator.update_position(&mut cont);
        cont.flush_forces();
        force::apply_lj_force_soa(&mut cont);
        integrator.update_velocity(&mut cont);

        print!("Step {}\r", i);
        std::io::stdout().flush().unwrap();
    }
    total_simulated_time += last_time.elapsed().as_secs_f64();
    println!("Successfully ran {} iterations in {} seconds.", ITERATIONS, total_simulated_time);
}



/// Runs a simple simulation using linked cell container.
/// Saves the results to CSV files for visualization.
#[allow(dead_code)]
fn run_linked_cell_simulation() {
    const LOG_TIMESTEPS: usize = 50; // Log every 50 timesteps
    const REDIST_FREQ: usize = 100; // Redistribute particles every 100 timesteps
    const ITERATIONS: usize = 50000; // Number of iterations to run
    const DELTA_T: f64 = 0.00005; // Timestep

    const TOTAL_SIZE: usize = SELECTED_SIM.0; // Total "size" of the simulation: number of cells in each dimension
    const CUTOFF: f64 = SELECTED_SIM.1; // Size of the container cells
    const GRID_SIZE: usize = SELECTED_SIM.2; // Number of particles in each dimension in the grid


    const X_SIZE: usize = TOTAL_SIZE;
    const Y_SIZE: usize = TOTAL_SIZE;
    const Z_SIZE: usize = TOTAL_SIZE;

    // The three boundary conditions for the simulation: don't leave the simulation box
    const X_BOUNDARY: BoundaryCondition<f64> = BoundaryCondition {
        dimension: boundary::Dimension::X,
        min: 0.0,
        max: X_SIZE as f64 * CUTOFF,
    };
    const Y_BOUNDARY: BoundaryCondition<f64> = BoundaryCondition {
        dimension: boundary::Dimension::Y,
        min: 0.0,
        max: Y_SIZE as f64 * CUTOFF,
    };
    const Z_BOUNDARY: BoundaryCondition<f64> = BoundaryCondition {
        dimension: boundary::Dimension::Z,
        min: 0.0,
        max: Z_SIZE as f64 * CUTOFF,
    };

    // Because of the way the grid is set up and the cutoff distance works, the Grid size divided by the total size 
    // should be less than the cutoff distance
    if GRID_SIZE as f64 / TOTAL_SIZE as f64 > CUTOFF {
        eprintln!("The grid size divided by the total size should be less than the cutoff distance. Particles will be lost on the edge; you might simulate fewer than expected!");
    }

    let mut total_simulated_time = 0.0; // records the total time spent simulating
    let mut last_time = std::time::Instant::now();

    let integrator = integrator::VelocityStörmerVerlet::new(DELTA_T); // Simulate with a very small timestep
    let mut cont = linkedcell::LinkedCell::<f64, X_SIZE, Y_SIZE, Z_SIZE>::new(CUTOFF);
    cont.add_particle_grid(GRID_SIZE, 0.80);
    cont.boundaries.push(X_BOUNDARY);
    cont.boundaries.push(Y_BOUNDARY);
    cont.boundaries.push(Z_BOUNDARY);

    force::apply_lj_force_linked_cell(&mut cont); // Calculate initial forces

    let mut log_step = 0;
    for i in 0..ITERATIONS {
        integrator.update_position_linked_par(&mut cont); // Update the position of the particles
        cont.flush_forces_par(); // Clear and swap the forces
        force::apply_lj_force_linked_cell(&mut cont); // Calculate the new forces
        
        integrator.update_velocity_linked_par(&mut cont); // Update the velocity of the particles
        
        cont.apply_boundaries(); // Apply the boundary conditions
        
        if i % REDIST_FREQ == 0 { // If it is time to do so,
            cont.redistribute_particles(); // Redistribute the particles
        }
        
        if i % LOG_TIMESTEPS == 0 { // Log the state of the simulation to a CSV file
            total_simulated_time += last_time.elapsed().as_secs_f64(); // "Stop" the timer

            linkedcell::LinkedCell::save_to_csv(&cont, log_step);
            log_step += 1;

            last_time = std::time::Instant::now(); // "Start" the timer
        }
        print!("Step {}\r", i); // This prints the result to the same line each time
        std::io::stdout().flush().unwrap();
    }
    total_simulated_time += last_time.elapsed().as_secs_f64();
    println!("Successfully ran {} iterations in {} seconds.", ITERATIONS, total_simulated_time);
}
