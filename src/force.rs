use crate::soacontainer::SoAContainer;

pub fn apply_lj_force<T,const N: usize>(container: &mut SoAContainer<T,N>) 
    where T: LjFloat + std::ops::AddAssign + std::ops::Sub<Output = T> + std::ops::Mul<Output = T> + std::ops::Add<Output = T> + std::ops::Div<Output = T> + Copy
{
    for i in 0..N {
        let mut fx: T = T::zero();
        let mut fy: T = T::zero();
        let mut fz: T = T::zero();
        for j in 0..N { // Maybe replace with par_iter??? The dependency between the particles makes it hard to parallelize
            if i == j {
                continue;
            }
            let dx = container.position.x[i] - container.position.x[j];
            let dy = container.position.y[i] - container.position.y[j];
            let dz = container.position.z[i] - container.position.z[j];
            
            // Lennard-Jones potential: 4 eps * (r^-12 - r^-6)

                        let r2inv = (dx*dx + dy*dy + dz*dz).recip();
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

pub trait LjFloat: Copy + num::Float + num::FromPrimitive{
    fn two () -> Self;
    fn twenty_four () -> Self;
}

impl<T> LjFloat for T where T: num::Float + num::FromPrimitive{
    fn two() -> Self {
        T::one() + T::one()
    }
    fn twenty_four() -> Self {
        T::two() * T::two() * T::two() * (T::one() + T::two())
    }
}