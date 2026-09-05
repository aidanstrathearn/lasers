use crate::dfb::LaserParameters;
use laser_solver::two_mode::TwoModeSolver;
use myplotlib::{AppResult, Plotter};

impl LaserParameters {
    pub fn kappa_plot(&mut self) -> AppResult {
        let points = {
            let fibre = self.resolved_fibre();
            let solver = TwoModeSolver::new(&fibre, self.steps);
            solver
                .grid()
                .positions()
                .zip(solver.kappas().iter().copied())
                .map(|(z, kappa)| [z, kappa])
                .collect()
        };

        let mut plot = Plotter::new();
        plot.add_points(points).label("Kappa");
        plot.xlabel("Position (m)");
        plot.ylabel("Kappa (1/m)");
        Ok(plot)
    }
}
