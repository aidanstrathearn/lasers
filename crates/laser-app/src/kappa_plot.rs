use crate::plotter::Plotter;
use crate::{dfb::DfbMode, timed};
use laser_solver::error::SolverError;
use laser_solver::two_mode::TwoModeSolver;

impl DfbMode {
    pub fn kappa_plot(&mut self) -> Result<Plotter, SolverError> {
        let (points, compute_time) = timed(|| {
            let fibre = self.resolved_fibre();
            let solver = TwoModeSolver::new(&fibre, self.steps);
            solver
                .grid()
                .positions()
                .zip(solver.kappas().iter().copied())
                .map(|(z, kappa)| [z, kappa])
                .collect()
        });
        self.compute_time = Some(compute_time);

        let mut plot = Plotter::new();
        plot.add_points(points).label("Kappa");
        plot.xlabel("Position (m)");
        plot.ylabel("Kappa (m⁻¹)");
        Ok(plot)
    }
}
