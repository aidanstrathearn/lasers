use crate::plotter::Plotter;
use crate::{LaserApp, Points, timed};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
use laser_solver::lase::FieldState;
use laser_solver::propagation::out_field_coupled;
use laser_solver::utils::linspace;

#[derive(Copy, Clone)]
pub struct ResidualRange {
    lower: f64,
    upper: f64,
    num: usize,
}

impl Default for ResidualRange {
    fn default() -> Self {
        Self {
            lower: 1e-6,
            upper: 20.0,
            num: 100,
        }
    }
}

impl LaserApp {
    pub fn residual_plot(&mut self) -> Result<Plotter, SolverError> {
        let inputs = linspace(
            self.residual_range.lower,
            self.residual_range.upper,
            self.residual_range.num,
        );

        let kappas = self.grating.grid(self.grid_points.0);
        let dz = self.grid_points.dz(self.fibre_params.length);
        let trial = |sgnl_b| FieldState {
            sgnl_f: 0.0,
            sgnl_b,
            pump_f: self.pump.forward_amplitude(),
            pump_b: 0.0,
        }; //todo: use picard for backward pump
        let f = |sgnl_b| {
            out_field_coupled(trial(sgnl_b), self.fibre_params, dz, &kappas).sgnl_b / sgnl_b
        };
        let mut compute_time = std::time::Duration::ZERO;
        let residuals = inputs
            .iter()
            .map(|&s| {
                let (residual, elapsed) = timed(|| f(s));
                compute_time += elapsed;
                [s, residual.abs().log10()]
            })
            .collect::<Points>();

        let mut plt = Plotter::new();
        plt.xlabel("Input");
        plt.ylabel("log(|residual|)");
        plt.add_points(residuals);
        plt.set_compute_time(compute_time);
        Ok(plt)
    }
}

pub fn residual_slider_grid(tr: &mut ResidualRange, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("residual").show(ui, |ui| {
        ui.label("high");
        changed |= ui
            .add(egui::Slider::new(&mut tr.upper, 1e-5..=100.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("low");
        changed |= ui
            .add(egui::Slider::new(&mut tr.lower, 1e-6..=tr.upper).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("num");
        changed |= ui
            .add(egui::Slider::new(&mut tr.num, 5..=100).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}
