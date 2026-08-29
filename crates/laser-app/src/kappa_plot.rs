use crate::plotter::Plotter;
use crate::{dfb::DfbMode, timed};
use laser_solver::error::SolverError;

impl DfbMode {
    pub fn kappa_plot(&mut self) -> Result<Plotter, SolverError> {
        let (points, compute_time) = timed(|| {
            let fibre = self.resolved_fibre();
            fibre
                .grid()
                .positions()
                .zip(fibre.kappas().iter().copied())
                .map(|(z, kappa)| [z, kappa])
                .collect()
        });
        self.compute_time = Some(compute_time);

        let mut plot = Plotter::new();
        plot.add_points(points).label("Kappa");
        plot.xlabel("z");
        plot.ylabel("Kappa");
        Ok(plot)
    }
}
