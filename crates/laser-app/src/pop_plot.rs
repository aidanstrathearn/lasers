use crate::plotter::Plotter;
use crate::{Points, dfb::DfbMode, timed};
use laser_solver::error::SolverError;
use laser_solver::maths::rootfind::BisectionConfig;

impl DfbMode {
    pub fn pops_plot(&mut self) -> Result<Plotter, SolverError> {
        let full_profile = true;
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };
        
        
        let (result, compute_time) = timed(|| {
            self.dfb_laser()
                .solve(self.pump, self.dfb_solve_config(bc), full_profile)
        });
        self.compute_time = Some(compute_time);
        let result = result?;
        let fibre = self.resolved_fibre();

        let (ground, excited): (Points, Points) = result
            .z
            .into_iter()
            .zip(result.fields)
            .map(|(z, field)| {
                let (g, e) = fibre.populations(field);
                ([z, g], [z, e])
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
