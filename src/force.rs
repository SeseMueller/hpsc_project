use crate::{linkedcell::LinkedCell, soacontainer::SoAContainer};

pub trait LjFloat: Copy + num::Float + num::FromPrimitive + std::fmt::Debug {
    fn two() -> Self;
    fn twenty_four() -> Self;
}

impl<T> LjFloat for T
where
    T: num::Float + num::FromPrimitive + std::fmt::Debug,
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
            // TODO: this was wrong, correct it
        }
    }
}

/// Applies the Lennard-Jones force to the particles of a dynamic soa container.
#[inline(never)] // To be able to see the function in the profiler
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
        container.current_force.x[i] += fx;
        container.current_force.y[i] += fy;
        container.current_force.z[i] += fz;
    }
}

/// Applies the Lennard-Jones force to the particles of both dynamic soa containers, but not within the same container.
/// If `newton_third_law` is false, the force is only applied to `container_a`.
#[inline(never)] // To be able to see the function in the profiler
pub fn apply_lj_force_2soa_dyn<T>(
    container_a: &mut crate::soacontainerdyn::SoAContainerDyn<T>,
    container_b: &mut crate::soacontainerdyn::SoAContainerDyn<T>,
    newton_third_law: bool,
) where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::SubAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    let num_particles_a = container_a.position.x.len();
    let num_particles_b = container_b.position.x.len();

    for i in 0..num_particles_a {
        for j in 0..num_particles_b {
            unsafe {
                // yes, I am aware. But the block needs to be unsafe because of the raw pointers
                // and it can't be split up into smaller blocks because the results are needed later

                let dx = *container_a.position.x.get_unchecked(i)
                    - *container_b.position.x.get_unchecked(j);
                let dy = *container_a.position.y.get_unchecked(i)
                    - *container_b.position.y.get_unchecked(j);
                let dz = *container_a.position.z.get_unchecked(i)
                    - *container_b.position.z.get_unchecked(j);

                // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

                let r2inv = (dx * dx + dy * dy + dz * dz).recip();
                let r6inv = r2inv * r2inv * r2inv;
                let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

                // fx += force_scalar * dx;
                // fy += force_scalar * dy;
                // fz += force_scalar * dz;
                *container_a.current_force.x.get_unchecked_mut(i) += force_scalar * dx;
                *container_a.current_force.y.get_unchecked_mut(i) += force_scalar * dy;
                *container_a.current_force.z.get_unchecked_mut(i) += force_scalar * dz;

                if newton_third_law {
                    *container_b.current_force.x.get_unchecked_mut(j) -= force_scalar * dx;
                    *container_b.current_force.y.get_unchecked_mut(j) -= force_scalar * dy;
                    *container_b.current_force.z.get_unchecked_mut(j) -= force_scalar * dz;
                }
            }
        }
    }
}

/// Applies the Lennard-Jones force to the particles of two soa containers by only using the necessary arrays.
#[inline(never)] // To be able to see the function in the profiler
pub fn apply_lj_force_arrays<T>(
    // container_a: &mut crate::soacontainerdyn::SoAContainerDyn<T>,
    // container_b: &mut crate::soacontainerdyn::SoAContainerDyn<T>,
    position_a_x: &[T],
    position_a_y: &[T],
    position_a_z: &[T],
    position_b_x: &[T],
    position_b_y: &[T],
    position_b_z: &[T],
    current_force_a_x: &mut [T],
    current_force_a_y: &mut [T],
    current_force_a_z: &mut [T],
    current_force_b_x: &mut [T],
    current_force_b_y: &mut [T],
    current_force_b_z: &mut [T],
) where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::SubAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    let num_particles_a = position_a_x.len();
    let num_particles_b = position_b_x.len();

    for i in 0..num_particles_a {
        for j in 0..num_particles_b {
            unsafe {
                // yes, I am aware. But the block needs to be unsafe because of the raw pointers
                // and it can't be split up into smaller blocks because the results are needed later

                let dx = *position_a_x.get_unchecked(i) - *position_b_x.get_unchecked(j);
                let dy = *position_a_y.get_unchecked(i) - *position_b_y.get_unchecked(j);
                let dz = *position_a_z.get_unchecked(i) - *position_b_z.get_unchecked(j);

                // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

                let r2inv = (dx * dx + dy * dy + dz * dz).recip();
                let r6inv = r2inv * r2inv * r2inv;
                let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

                *current_force_a_x.get_unchecked_mut(i) += force_scalar * dx;
                *current_force_a_y.get_unchecked_mut(i) += force_scalar * dy;
                *current_force_a_z.get_unchecked_mut(i) += force_scalar * dz;

                *current_force_b_x.get_unchecked_mut(j) -= force_scalar * dx;
                *current_force_b_y.get_unchecked_mut(j) -= force_scalar * dy;
                *current_force_b_z.get_unchecked_mut(j) -= force_scalar * dz;
            }
        }
    }
}

/// Applies the Lennard-Jones force to the particles of both dynamic soa containers, but not within the same container.
/// If `newton_third_law` is false, the force is only applied to `container_a`.
/// Recieves the indexes of the containers to apply the force to.
pub fn apply_lj_force_2soa_dyn_index<T, const D0: usize, const D1: usize, const D2: usize>(
    cell: &mut LinkedCell<T, D0, D1, D2>,
    container_a_index: (usize, usize, usize),
    container_b_index: (usize, usize, usize),
    newton_third_law: bool,
) where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::SubAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    let num_particles_a = cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
        .position
        .x
        .len();
    let num_particles_b = cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
        .position
        .x
        .len();

    // First set the forces to zero, they will be accumulated in the next loop

    for i in 0..num_particles_a {
        for j in 0..num_particles_b {
            let dx = cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
                .position
                .x[i]
                - cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
                    .position
                    .x[j];
            let dy = cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
                .position
                .y[i]
                - cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
                    .position
                    .y[j];
            let dz = cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
                .position
                .z[i]
                - cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
                    .position
                    .z[j];

            // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

            let r2inv = (dx * dx + dy * dy + dz * dz).recip();
            let r6inv = r2inv * r2inv * r2inv;
            let force_scalar = T::twenty_four() * r6inv * ((T::two() * r6inv) - T::one());

            // fx += force_scalar * dx;
            // fy += force_scalar * dy;
            // fz += force_scalar * dz;
            cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
                .current_force
                .x[i] += force_scalar * dx;
            cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
                .current_force
                .y[i] += force_scalar * dy;
            cell.cells[container_a_index.0][container_a_index.1][container_a_index.2]
                .current_force
                .z[i] += force_scalar * dz;

            if newton_third_law {
                cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
                    .current_force
                    .x[j] -= force_scalar * dx;
                cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
                    .current_force
                    .y[j] -= force_scalar * dy;
                cell.cells[container_b_index.0][container_b_index.1][container_b_index.2]
                    .current_force
                    .z[j] -= force_scalar * dz;
            }
        }
    }
}

/// Applies the Lennard-Jones force to the particles in two soa containers of the linked cell.
/// This specialization is needed because the pointer aliasing rule forbids

/// Applies the Lennard-Jones force to the particles in the linked cell.
pub fn apply_lj_force_linked_cell<T, const D0: usize, const D1: usize, const D2: usize>(
    cell: &mut LinkedCell<T, D0, D1, D2>,
) where
    T: LjFloat
        + std::ops::AddAssign
        + std::ops::SubAssign
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Div<Output = T>
        + Copy,
{
    // We need to loop over all cells and apply the force to the particles in the cell
    // AS WELL AS the force to the particles in the neighboring cells

    // We have 26 neighboring cells, for each of them we need to apply the force
    // In order to avoid double counting, we only consider the 13 cells with a "higher index"
    const MINUS1: usize = usize::MAX;
    let stencil = [
        (0, 0, 1),
        (0, 1, 0),
        (0, 1, 1),
        (0, 1, MINUS1),
        (1, 0, 0),
        (1, 0, 1),
        (1, 0, MINUS1),
        (1, 1, 0),
        (1, 1, 1),
        (1, 1, MINUS1),
        (1, MINUS1, 0),
        (1, MINUS1, 1),
        (1, MINUS1, MINUS1),
    ];

    for i in 0..D0 {
        for j in 0..D1 {
            for k in 0..D2 {
                for (delta_0, delta_1, delta_2) in stencil.iter() {
                    let inner_0 = usize::overflowing_add(i, *delta_0).0;
                    let inner_1 = usize::overflowing_add(j, *delta_1).0;
                    let inner_2 = usize::overflowing_add(k, *delta_2).0;

                    if inner_0 >= D0 || inner_1 >= D1 || inner_2 >= D2 {
                        continue; // Skip if we are out of bounds
                    }

                    // Passing the same array twive with mut would violate the pointer aliasing rule,
                    // so we need to use raw pointers

                    // let container_a = &mut cell.cells[i][j][k] as *mut _;
                    // let container_b = &mut cell.cells[inner_0][inner_1][inner_2] as *mut _;
                    // apply_lj_force_2soa_dyn(
                    //     unsafe { &mut *container_a }, // Safety: we know that the pointer is valid and that
                    //     unsafe { &mut *container_b }, // this pointer is pointing somewhere else in the grid
                    //     true,
                    // );

                    // There was some bug there, so I'm trying to do it manually
                    // apply_lj_force_2soa_dyn_index(
                    //     cell,
                    //     (i, j, k),
                    //     (inner_0, inner_1, inner_2),
                    //     true,
                    // );
                    // Nope, both give the same result
                    // Found the error: the single container version overwrote the forces of the other container (= instead of +=)

                    // In an attempt to be faster, I removed unnecessary indirection, but the setup is a bit longer
                    // Basically, the force interactions are an operation between twelve arrays, so we can just pass the arrays
                    // We still need unsafe though.
                    let positions_a_x = &cell.cells[i][j][k].position.x;
                    let positions_a_y = &cell.cells[i][j][k].position.y;
                    let positions_a_z = &cell.cells[i][j][k].position.z;
                    let positions_b_x = &cell.cells[inner_0][inner_1][inner_2].position.x;
                    let positions_b_y = &cell.cells[inner_0][inner_1][inner_2].position.y;
                    let positions_b_z = &cell.cells[inner_0][inner_1][inner_2].position.z;
                    let forces_a_x = &mut *cell.cells[i][j][k].current_force.x as *mut _;
                    let forces_a_y = &mut *cell.cells[i][j][k].current_force.y as *mut _;
                    let forces_a_z = &mut *cell.cells[i][j][k].current_force.z as *mut _;
                    let forces_b_x =
                        &mut *cell.cells[inner_0][inner_1][inner_2].current_force.x as *mut _;
                    let forces_b_y =
                        &mut *cell.cells[inner_0][inner_1][inner_2].current_force.y as *mut _;
                    let forces_b_z =
                        &mut *cell.cells[inner_0][inner_1][inner_2].current_force.z as *mut _;

                    apply_lj_force_arrays(
                        positions_a_x,
                        positions_a_y,
                        positions_a_z,
                        positions_b_x,
                        positions_b_y,
                        positions_b_z,
                        unsafe { &mut *forces_a_x },
                        unsafe { &mut *forces_a_y },
                        unsafe { &mut *forces_a_z },
                        unsafe { &mut *forces_b_x },
                        unsafe { &mut *forces_b_y },
                        unsafe { &mut *forces_b_z },
                    );
                }

                // Now we need to apply the force to the particles in the cell itself
                apply_lj_force_soa_dyn(&mut cell.cells[i][j][k]);
            }
        }
    }
}
