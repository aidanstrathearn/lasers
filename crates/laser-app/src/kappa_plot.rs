use crate::plotter::Plotter;
use crate::{dfb::DfbMode, timed};
use laser_solver::error::SolverError;

impl DfbMode {
    pub fn kappa_plot(&mut self) -> Result<Plotter, SolverError> {
        let ((z, kappas), compute_time) = timed(|| {
            (
                self.grid_points.grid(self.fibre.length()),
                self.grating.grid(self.grid_points.0),
            )
        });
        self.compute_time = Some(compute_time);

        // kappas.len() = z.len() - 1, so zip ignores last element of z
        let points = z.into_iter().zip(kappas).map(|(z, k)| [z, k]).collect();

        let mut plot = Plotter::new();
        plot.add_points(points).label("Kappa");
        plot.xlabel("z");
        plot.ylabel("Kappa");
        Ok(plot)
    }
}
