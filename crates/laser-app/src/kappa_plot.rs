use crate::plotter::Plotter;
use crate::{LaserApp, timed};
use laser_solver::error::SolverError;

impl LaserApp {
    pub fn kappa_plot(&self) -> Result<Plotter, SolverError> {
        let ((z, kappas), compute_time) = timed(|| {
            (
                self.grid_points.grid(self.fibre.length),
                self.grating.grid(self.grid_points.0),
            )
        });

        // kappas.len() = z.len() - 1, so zip ignores last element of z
        let points = z.into_iter().zip(kappas).map(|(z, k)| [z, k]).collect();

        let mut plot = Plotter::new();
        plot.add_points(points).label("Kappa");
        plot.xlabel("z");
        plot.ylabel("Kappa");
        plot.set_compute_time(compute_time);
        Ok(plot)
    }
}
