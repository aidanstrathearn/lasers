use crate::plotter::Plotter;
use crate::{dfb::DfbMode, timed};
use laser_solver::error::SolverError;
use laser_solver::grating::sample_grating;
use laser_solver::lase::UniformGrid;

impl DfbMode {
    pub fn kappa_plot(&mut self) -> Result<Plotter, SolverError> {
        let ((z, kappas), compute_time) = timed(|| {
            let grid = UniformGrid::new(self.fibre.geometry.length, self.steps);
            (
                grid.positions().collect::<Vec<_>>(),
                sample_grating(&self.grating, grid.steps()),
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
