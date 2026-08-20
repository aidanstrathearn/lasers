use crate::plotter::Plotter;
use crate::{LaserApp, Points};
use eframe::egui::Ui;
use laser_solver::dfb::dfb_solve;
use laser_solver::error::SolverError;
use laser_solver::lase::{Pump, pops};
use laser_solver::picard::PicardConfig;
use laser_solver::rootfind::BisectionConfig;

impl LaserApp {
    pub fn pops_plot(&mut self, ui: &mut Ui) -> Result<(), SolverError> {
        let full_profile = true;
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };

        let picard_config = PicardConfig {
            max_iterations: 5_000,
            relative_tolerance: 1e-6,
            absolute_tolerance: 1e-10,
        };
        let pu = Pump::from_total_and_balance(self.pump.total, self.pump.balance);
        let result = dfb_solve(
            pu,
            self.fibre_params,
            self.grid_points,
            self.grating,
            full_profile,
            bc,
            picard_config,
        )?;

        let (ground, excited): (Points, Points) = result
            .z
            .into_iter()
            .zip(result.fields)
            .map(|(z, field)| {
                let (g, e) = pops(field, self.fibre_params);
                ([z, g], [z, e])
            })
            .unzip();

        let mut plot = Plotter::new();
        plot.add_points(ground).label("Ground state");
        plot.add_points(excited).label("Excited state");
        plot.xlabel("z");
        plot.ylabel("Population fraction");
        plot.show(ui, "population-plot");

        Ok(())
    }
}
