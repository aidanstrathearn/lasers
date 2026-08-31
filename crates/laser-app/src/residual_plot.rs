use crate::plotter::Plotter;
use crate::{Points, dfb::DfbMode, timed};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
use laser_solver::maths::utils::linspace;
use laser_solver::two_mode::TwoModeSolver;

#[derive(Copy, Clone)]
pub struct ResidualRange {
    lower: f64,
    upper: f64,
    num: usize,
}

impl Default for ResidualRange {
    fn default() -> Self {
        Self {
            lower: 1e-12,
            upper: 400.0,
            num: 100,
        }
    }
}

impl DfbMode {
    pub fn residual_plot(&mut self) -> Result<Plotter, SolverError> {
        let trial_fluxes = linspace(
            self.residual_range.lower,
            self.residual_range.upper,
            self.residual_range.num,
        );

        let fibre = self.resolved_fibre();
        let solver = TwoModeSolver::new(&fibre, self.steps);
        let (residuals, compute_time) =
            timed(|| solver.shooting_residuals(self.pump, &trial_fluxes));
        let residual_points = trial_fluxes
            .into_iter()
            .zip(residuals)
            .map(|(trial_flux, residual)| [trial_flux, residual.abs().log10()])
            .collect::<Points>();
        self.compute_time = Some(compute_time);

        let mut plt = Plotter::new();
        plt.xlabel("Trial backward signal flux (1e25 photons/(m² s))");
        plt.ylabel("log(|residual|)");
        plt.add_points(residual_points);
        Ok(plt)
    }
}

pub fn residual_slider_grid(tr: &mut ResidualRange, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("residual").show(ui, |ui| {
        ui.label("high");
        changed |= ui
            .add(egui::Slider::new(&mut tr.upper, 1e-10..=10_000.0).logarithmic(true))
            .changed();
        ui.end_row();

        ui.label("low");
        changed |= ui
            .add(egui::Slider::new(&mut tr.lower, 1e-12..=tr.upper).logarithmic(true))
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
