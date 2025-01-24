use crate::soacontainer::SoAContainer;

pub trait LjFloat: Copy + num::Float + num::FromPrimitive {
    fn two() -> Self;
    fn twenty_four() -> Self;
}

impl<T> LjFloat for T
where
    T: num::Float + num::FromPrimitive,
{
    fn two() -> Self {
        T::one() + T::one()
    }
    fn twenty_four() -> Self {
        T::two() * T::two() * T::two() * (T::one() + T::two())
    }
}

/// Applies the Lennard-Jones force to the particles in the container.
pub fn apply_lj_force_soa<T, const N: usize>(container: &mut SoAContainer<T, N>)
where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    for i in 0..N {
        let mut fx: T = T::zero();
        let mut fy: T = T::zero();
        let mut fz: T = T::zero();
        for j in 0..N {
            // Maybe replace with par_iter??? The dependency between the particles makes it hard to parallelize
            if i == j {
                continue;
            }
            let dx = container.position.x[i] - container.position.x[j];
            let dy = container.position.y[i] - container.position.y[j];
            let dz = container.position.z[i] - container.position.z[j];

            // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

            let r2inv = (dx * dx + dy * dy + dz * dz).recip();
            let r6inv = r2inv * r2inv * r2inv;
            let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

            fx += force_scalar * dx;
            fy += force_scalar * dy;
            fz += force_scalar * dz;
        }
        container.current_force.x[i] = fx;
        container.current_force.y[i] = fy;
        container.current_force.z[i] = fz;
    }
}

/// Applies the Lennard-Jones force to the particles of both containers, but not within the same container.
/// If `newton_third_law` is false, the force is only applied to `container_a`.
pub fn apply_lj_force_2soa<T, const N: usize, const M: usize>(
    container_a: &mut SoAContainer<T, N>,
    container_b: &mut SoAContainer<T, M>,
    newton_third_law: bool,
) where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    for i in 0..N {
        let mut fx: T = T::zero();
        let mut fy: T = T::zero();
        let mut fz: T = T::zero();
        for j in 0..M {
            let dx = container_a.position.x[i] - container_b.position.x[j];
            let dy = container_a.position.y[i] - container_b.position.y[j];
            let dz = container_a.position.z[i] - container_b.position.z[j];

            // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

            let r2inv = (dx * dx + dy * dy + dz * dz).recip();
            let r6inv = r2inv * r2inv * r2inv;
            let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

            fx += force_scalar * dx;
            fy += force_scalar * dy;
            fz += force_scalar * dz;
        }
        container_a.current_force.x[i] = fx;
        container_a.current_force.y[i] = fy;
        container_a.current_force.z[i] = fz;

        if newton_third_law {
            container_b.current_force.x[i] = -fx;
            container_b.current_force.y[i] = -fy;
            container_b.current_force.z[i] = -fz;
        }
    }
}

/// Applies the Lennard-Jones force to the particles of a dynamic soa container.
pub fn apply_lj_force_soa_dyn<T>(container: &mut crate::soacontainerdyn::SoAContainerDyn<T>)
where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    let num_particles = container.position.x.len();
    for i in 0..num_particles {
        let mut fx: T = T::zero();
        let mut fy: T = T::zero();
        let mut fz: T = T::zero();
        for j in 0..num_particles {
            if i == j {
                continue;
            }
            let dx = container.position.x[i] - container.position.x[j];
            let dy = container.position.y[i] - container.position.y[j];
            let dz = container.position.z[i] - container.position.z[j];

            // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

            let r2inv = (dx * dx + dy * dy + dz * dz).recip();
            let r6inv = r2inv * r2inv * r2inv;
            let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

            fx += force_scalar * dx;
            fy += force_scalar * dy;
            fz += force_scalar * dz;
        }
        container.current_force.x[i] = fx;
        container.current_force.y[i] = fy;
        container.current_force.z[i] = fz;
    }
}

/// Applies the Lennard-Jones force to the particles of both dynamic soa containers, but not within the same container.
/// If `newton_third_law` is false, the force is only applied to `container_a`.
pub fn apply_lj_force_2soa_dyn<T>(
    container_a: &mut crate::soacontainerdyn::SoAContainerDyn<T>,
    container_b: &mut crate::soacontainerdyn::SoAContainerDyn<T>,
    newton_third_law: bool,
) where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    let num_particles_a = container_a.position.x.len();
    let num_particles_b = container_b.position.x.len();
    for i in 0..num_particles_a {
        let mut fx: T = T::zero();
        let mut fy: T = T::zero();
        let mut fz: T = T::zero();
        for j in 0..num_particles_b {
            let dx = container_a.position.x[i] - container_b.position.x[j];
            let dy = container_a.position.y[i] - container_b.position.y[j];
            let dz = container_a.position.z[i] - container_b.position.z[j];

            // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

            let r2inv = (dx * dx + dy * dy + dz * dz).recip();
            let r6inv = r2inv * r2inv * r2inv;
            let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

            fx += force_scalar * dx;
            fy += force_scalar * dy;
            fz += force_scalar * dz;
        }
        container_a.current_force.x[i] = fx;
        container_a.current_force.y[i] = fy;
        container_a.current_force.z[i] = fz;

        if newton_third_law {
            for j in 0..num_particles_b {
                let dx = container_a.position.x[i] - container_b.position.x[j];
                let dy = container_a.position.y[i] - container_b.position.y[j];
                let dz = container_a.position.z[i] - container_b.position.z[j];

                // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

                let r2inv = (dx * dx + dy * dy + dz * dz).recip();
                let r6inv = r2inv * r2inv * r2inv;
                let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

                container_b.current_force.x[j] = -force_scalar * dx;
                container_b.current_force.y[j] = -force_scalar * dy;
                container_b.current_force.z[j] = -force_scalar * dz;
            }
        }
    }
}

/// Applies the Lennard-Jones force to the particles in the linked cell.
/// TODO!!
pub fn TODO() {
    unimplemented!();
}