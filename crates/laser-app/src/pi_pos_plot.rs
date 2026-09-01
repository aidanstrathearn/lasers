use crate::plotter::Plotter;
use crate::{LaserParameters, Points, dfb::DfbMode, timed};
use laser_solver::error::SolverError;
use laser_solver::grating::PiShift;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::maths::utils::linspace;
use laser_solver::two_mode::TwoModeSolver;
use std::time::Duration;

const PI_POSITION_INTERVALS: usize = 40;

impl DfbMode {
    pub fn pi_pos_plot(&mut self, parameters: &LaserParameters) -> Result<Plotter, SolverError> {
        let pi_positions = linspace(0.0, 1.0, PI_POSITION_INTERVALS);
        let mut forward_output: Points = Vec::with_capacity(pi_positions.len());
        let mut backward_output: Points = Vec::with_capacity(pi_positions.len());
        let mut compute_time = Duration::ZERO;

        for pi_position in pi_positions {
            let grating = PiShift {
                pi_shift_position: pi_position,
                ..parameters.fibre.grating
            };

            let (profile, elapsed) = timed(|| {
                let fibre = laser_solver::lase::Fibre {
                    grating,
                    ..parameters.fibre.clone()
                };
                let fibre = fibre.resolve_with_interactions(
                    parameters.pump_mode,
                    parameters.pump_interaction,
                    parameters.sgnl_mode,
                    parameters.signal_interaction,
                );
                let bc = BisectionConfig {
                    upper: 2.0 * fibre.pump_flux(parameters.pump.total).sqrt(),
                    ..self.config
                };
                TwoModeSolver::new(&fibre, self.steps).solve_lasing(
                    parameters.pump,
                    bc.into(),
                    self.picard_config,
                    false,
                )
            });
            compute_time += elapsed;

            if let Ok(profile) = profile {
                let (forward, backward) = profile.output_powers();
                forward_output.push([pi_position, 1_000.0 * forward]);
                backward_output.push([pi_position, 1_000.0 * backward]);
            }
        }
        self.compute_time = Some(compute_time);

        let mut plot = Plotter::new();
        plot.add_points(forward_output).label("Forward");
        plot.add_points(backward_output).label("Backward");
        plot.axvline(parameters.fibre.grating.pi_shift_position)
            .label("Current position");
        plot.xlabel("Pi shift position");
        plot.ylabel("Output power (mW)");
        plot.xlim(0.0, 1.0);
        Ok(plot)
    }
}
