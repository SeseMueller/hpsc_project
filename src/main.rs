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
    // run_soa_dyn_simulation();
    run_linked_cell_simulation();
}

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

fn run_linked_cell_simulation() {
    const LOG_TIMESTEPS: usize = 50;
    const REDIST_FREQ: usize = 100;
    const ITERATIONS: usize = 5000;

    const X_SIZE: usize = 40;
    const Y_SIZE: usize = 40;
    const Z_SIZE: usize = 40;
    const CUTOFF: f64 = 1.5;

    // The three boundary conditions for the simulation
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


    let mut total_simulated_time = 0.0;
    let mut last_time = std::time::Instant::now();

    let integrator = integrator::VelocityStörmerVerlet::new(0.0005);
    let mut cont = linkedcell::LinkedCell::<f64, X_SIZE, Y_SIZE, Z_SIZE>::new(CUTOFF);
    cont.add_particle_grid(40, 0.9);
    cont.boundaries.push(X_BOUNDARY);
    cont.boundaries.push(Y_BOUNDARY);
    cont.boundaries.push(Z_BOUNDARY);
    force::apply_lj_force_linked_cell(&mut cont); // Calculate initial forces

    let mut log_step = 0;
    for i in 0..ITERATIONS {
        if i % LOG_TIMESTEPS == 0 {
            total_simulated_time += last_time.elapsed().as_secs_f64();

            linkedcell::LinkedCell::save_to_csv(&cont, log_step);
            log_step += 1;

            last_time = std::time::Instant::now();
        }
        integrator.update_position_linked(&mut cont);
        cont.flush_forces();
        force::apply_lj_force_linked_cell(&mut cont);
        // debug_check_linked_cell(&cont);
        integrator.update_velocity_linked(&mut cont);

        cont.apply_boundaries();

        if i % REDIST_FREQ == 0 {
            cont.redistribute_particles();
        }

        print!("Step {}\r", i);
        std::io::stdout().flush().unwrap();
    }
    total_simulated_time += last_time.elapsed().as_secs_f64();
    println!("Successfully ran {} iterations in {} seconds.", ITERATIONS, total_simulated_time);
}

// fn debug_check_linked_cell<T>(cont: &linkedcell::LinkedCell<T, 10, 10, 10>)
// where
//     T: force::LjFloat,
// {
    // // Just for debugging, look at the forces within the linked cell
    // for i in 0..10 {
    //     for j in 0..10 {
    //         for k in 0..10 {
    //             let cell = &cont.cells[i][j][k];
    //             for p in 0..cell.position.x.len() {
    //                 if cell.current_force.x[p] != T::zero() {
    //                     println!("Non-zero force after calculation: {:?}", cell.current_force.x[p]);
    //                 }
    //             }
    //         }
    //     }
    // }
// }