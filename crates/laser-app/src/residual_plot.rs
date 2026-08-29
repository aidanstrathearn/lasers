use crate::plotter::Plotter;
use crate::{Points, dfb::DfbMode, timed};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
use laser_solver::lase::{BidirectionalAmplitude, FieldState};
use laser_solver::propagation::out_field_coupled;
use laser_solver::maths::utils::linspace;

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

impl DfbMode {
    pub fn residual_plot(&mut self) -> Result<Plotter, SolverError> {
        let inputs = linspace(
            self.residual_range.lower,
            self.residual_range.upper,
            self.residual_range.num,
        );

        let fibre = self.resolved_fibre();
        let grid = fibre.grid();
        let kappas = fibre.kappas();
        let dz = grid.dz();
        let trial = |sgnl_b| FieldState {
            signal: BidirectionalAmplitude {
                forward: 0.0,
                backward: sgnl_b,
            },
            pump: BidirectionalAmplitude {
                forward: self.pump.forward_amplitude(),
                backward: 0.0,
            },
        }; //todo: use picard for backward pump
        let f = |sgnl_b| {
            out_field_coupled(trial(sgnl_b), |fields| fibre.gain(fields), dz, kappas)
                .signal
                .backward
                / sgnl_b
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
        self.compute_time = Some(compute_time);

        let mut plt = Plotter::new();
        plt.xlabel("Input");
        plt.ylabel("log(|residual|)");
        plt.add_points(residuals);
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
