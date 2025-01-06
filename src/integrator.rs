use crate::soacontainer::SoAContainer;

pub struct VelocityStörmerVerlet<T> {
    pub dt: T,
    pub half_dt: T,
    pub half_dt_sq: T,
}

impl<T> VelocityStörmerVerlet<T> where T: num::Float {
    pub fn new(dt: T) -> Self {
        Self {
            dt,
            half_dt: dt / (T::one() + T::one()),
            half_dt_sq: dt * dt / (T::one() + T::one()),
        }
    }
}
impl<T> VelocityStörmerVerlet<T> where T: num::Float + std::ops::AddAssign {

    #[inline(never)] // To be able to see the function in the profiler
    pub fn update_position< const N: usize> (&self, container: &mut SoAContainer<T, N>) {
        for i in 0..N {
            container.position.x[i] += container.velocity.x[i] * self.dt + container.current_force.x[i] * self.half_dt_sq;
            container.position.y[i] += container.velocity.y[i] * self.dt + container.current_force.y[i] * self.half_dt_sq;
            container.position.z[i] += container.velocity.z[i] * self.dt + container.current_force.z[i] * self.half_dt_sq;
        }
    }
    // Note that neither update_position nor update_velocity can be trivially parallelized because of the ownership structure of the SoAContainer.
    #[inline(never)] // To be able to see the function in the profiler
    pub fn update_velocity< const N: usize> (&self, container: &mut SoAContainer<T, N>) {
        for i in 0..N {
            container.velocity.x[i] += (container.current_force.x[i] + container.old_force.x[i]) * self.half_dt;
            container.velocity.y[i] += (container.current_force.y[i] + container.old_force.y[i]) * self.half_dt;
            container.velocity.z[i] += (container.current_force.z[i] + container.old_force.z[i]) * self.half_dt;
        }
    }
}