use crate::plotter::Plotter;
use crate::{LaserApp, Points};
use eframe::egui::Ui;
use laser_solver::dfb::dfb_solve;
use laser_solver::lase::{GratingProfile, Pump};
use laser_solver::picard::PicardConfig;
use laser_solver::rootfind::BisectionConfig;
use laser_solver::utils::linspace;

const PI_POSITION_INTERVALS: usize = 40;

impl LaserApp {
    pub fn pi_pos_plot(&mut self, ui: &mut Ui) {
        let pump = Pump::from_total_and_balance(self.pump.total, self.pump.balance);
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };
        let picard_config = PicardConfig {
            max_iterations: 5_000,
            relative_tolerance: 1e-6,
            absolute_tolerance: 1e-10,
        };

        let pi_positions = linspace(0.0, 1.0, PI_POSITION_INTERVALS);
        let mut forward_output: Points = Vec::with_capacity(pi_positions.len());
        let mut backward_output: Points = Vec::with_capacity(pi_positions.len());

        for pi_position in pi_positions {
            let grating = GratingProfile {
                pi_shift_position: pi_position,
                ..self.grating
            };

            if let Ok(profile) = dfb_solve(
                pump,
                self.fibre_params,
                self.grid_points,
                grating,
                false,
                bc,
                picard_config,
            ) {
                let (forward, backward) = profile.output_powers();
                forward_output.push([pi_position, forward]);
                backward_output.push([pi_position, backward]);
            }
        }

        let mut plot = Plotter::new();
        plot.add_points(forward_output).label("Forward");
        plot.add_points(backward_output).label("Backward");
        plot.xlabel("Pi shift position");
        plot.ylabel("Output power");
        plot.xlim(0.0, 1.0);
        plot.show(ui, "pi-position-output-plot");
    }
}
