use crate::plotter::Plotter;
use crate::{LaserParameters, Points, dfb::DfbMode, timed};
use laser_solver::error::SolverError;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::two_mode::TwoModeSolver;

impl DfbMode {
    pub fn pops_plot(&mut self, parameters: &LaserParameters) -> Result<Plotter, SolverError> {
        let full_profile = true;

        let (result, compute_time) = timed(|| {
            let fibre = parameters.resolved_fibre();
            let bc = BisectionConfig {
                upper: 2.0 * fibre.pump_flux(parameters.pump.total).sqrt(),
                ..self.config
            };
            TwoModeSolver::new(&fibre, self.steps).solve_lasing(
                parameters.pump,
                bc.into(),
                self.picard_config,
                full_profile,
            )
        });
        self.compute_time = Some(compute_time);
        let result = result?;
        let fibre = parameters.resolved_fibre();
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
