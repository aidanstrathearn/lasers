use crate::plotter::Plotter;
use crate::{Points, dfb::DfbMode, timed};
use laser_solver::error::SolverError;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::two_mode::TwoModeSolver;

impl DfbMode {
    pub fn pops_plot(&mut self) -> Result<Plotter, SolverError> {
        let full_profile = true;
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };

        let (result, compute_time) = timed(|| {
            let fibre = self.resolved_fibre();
            TwoModeSolver::new(&fibre, self.steps).solve_lasing(
                self.pump,
                bc.into(),
                self.picard_config,
                full_profile,
            )
        });
        self.compute_time = Some(compute_time);
        let result = result?;
        let fibre = self.resolved_fibre();

        let (ground, excited): (Points, Points) = result
            .z
            .into_iter()
            .zip(result.fields)
            .map(|(z, field)| {
                let populations = fibre.populations(field);
                ([z, populations.ground], [z, populations.excited])
            })
            .unzip();

        let mut plot = Plotter::new();
        plot.add_points(ground).label("Ground state");
        plot.add_points(excited).label("Excited state");
        plot.xlabel("z");
        plot.ylabel("Population fraction");
        Ok(plot)
    }
}
