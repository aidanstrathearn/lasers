use crate::dfb::LaserParameters;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::two_mode::TwoModeSolver;
use plot_app::{AppResult, Plotter, Points};

impl LaserParameters {
    pub fn pops_plot(&mut self) -> AppResult {
        let full_profile = true;

        let result = {
            let fibre = self.resolved_fibre();
            let bc = BisectionConfig {
                upper: 2.0 * fibre.pump_flux(self.pump.total).sqrt(),
                ..self.config
            };
            TwoModeSolver::new(&fibre, self.steps).solve_lasing(
                self.pump,
                bc.into(),
                self.picard_config,
                full_profile,
            )
        };
        let result = result?;
        let fibre = self.resolved_fibre();
        let populations = fibre.profile_populations(&result);

        let (ground, excited): (Points, Points) = result
            .z()
            .zip(populations)
            .map(|(z, populations)| ([z, populations.ground], [z, populations.excited]))
            .unzip();

        let mut plot = Plotter::new();
        plot.add_points(ground).label("Ground state");
        plot.add_points(excited).label("Excited state");
        plot.xlabel("Position (m)");
        plot.ylabel("Population fraction");
        Ok(plot)
    }
}
