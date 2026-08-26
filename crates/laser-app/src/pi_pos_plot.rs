use crate::plotter::Plotter;
use crate::{LaserApp, Points, timed};
use laser_solver::dfb::{dfb_solve, GratingProfile};
use laser_solver::error::SolverError;
use laser_solver::rootfind::BisectionConfig;
use laser_solver::utils::linspace;
use std::time::Duration;

const PI_POSITION_INTERVALS: usize = 40;

impl LaserApp {
    pub fn pi_pos_plot(&mut self) -> Result<Plotter, SolverError> {
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };
        let pi_positions = linspace(0.0, 1.0, PI_POSITION_INTERVALS);
        let mut forward_output: Points = Vec::with_capacity(pi_positions.len());
        let mut backward_output: Points = Vec::with_capacity(pi_positions.len());
        let mut compute_time = Duration::ZERO;

        for pi_position in pi_positions {
            let grating = GratingProfile {
                pi_shift_position: pi_position,
                ..self.grating
            };

            let (profile, elapsed) = timed(|| {
                dfb_solve(
                    self.pump,
                    self.fibre_params,
                    self.grid_points,
                    grating,
                    false,
                    bc,
                    self.picard_config,
                )
            });
            compute_time += elapsed;

            if let Ok(profile) = profile {
                let (forward, backward) = profile.output_powers();
                forward_output.push([pi_position, forward]);
                backward_output.push([pi_position, backward]);
            }
        }

        let mut plot = Plotter::new();
        plot.add_points(forward_output).label("Forward");
        plot.add_points(backward_output).label("Backward");
        plot.axvline(self.grating.pi_shift_position)
            .label("Current position");
        plot.xlabel("Pi shift position");
        plot.ylabel("Output power");
        plot.xlim(0.0, 1.0);
        plot.set_compute_time(compute_time);
        Ok(plot)
    }
}
