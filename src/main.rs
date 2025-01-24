mod soavector;

mod soacontainer;

use std::io::Write;

use soacontainer::SoAContainer;

mod integrator;

mod force;

mod boundary;

mod linkedcell;

mod soacontainerdyn;

fn main() {
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
