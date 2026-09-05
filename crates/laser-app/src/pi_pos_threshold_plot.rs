use crate::dfb::LaserParameters;
use laser_solver::grating::PiShift;
use laser_solver::lase::Fibre;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::maths::utils::linspace;
use laser_solver::two_mode::TwoModeSolver;
use myplotlib::{AppResult, Plotter, Points};

const PI_POSITION_INTERVALS: usize = 40;

impl LaserParameters {
    pub fn pi_pos_threshold_plot(&mut self) -> AppResult {
        let pi_positions = linspace(0.0, 1.0, PI_POSITION_INTERVALS);
        let mut thresholds: Points = Vec::with_capacity(pi_positions.len());

        for pi_position in pi_positions {
            let grating = PiShift {
                pi_shift_position: pi_position,
                ..self.fibre.grating
            };

            let threshold = {
                let fibre = Fibre {
                    grating,
                    ..self.fibre.clone()
                };
                let fibre = fibre.resolve_with_interactions(
                    self.pump_mode,
                    self.pump_interaction,
                    self.sgnl_mode,
                    self.signal_interaction,
                );
                let threshold_config = self.threshold_range.bisection_config(self.config);
                let pump_boundary_config = BisectionConfig {
                    lower: 1e-8,
                    upper: 1.0,
                    ..self.config
                };

                TwoModeSolver::new(&fibre, self.steps).find_threshold(
                    self.pump.balance,
                    threshold_config.into(),
                    pump_boundary_config.into(),
                )
            };

            if let Ok(threshold) = threshold {
                thresholds.push([pi_position, 1_000.0 * threshold]);
            }
        }
        let mut plot = Plotter::new();
        plot.add_points(thresholds).label("Threshold");
        plot.axvline(self.fibre.grating.pi_shift_position)
            .label("Current position");
        plot.xlabel("Pi shift position");
        plot.ylabel("Threshold pump power (mW)");
        plot.xlim(0.0, 1.0);
        Ok(plot)
    }
}
