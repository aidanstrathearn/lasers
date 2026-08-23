use crate::plotter::Plotter;
use crate::{LaserApp, Points};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::dfb::{dfb_pump_scan, out_field};
use laser_solver::lase::FieldState;
use laser_solver::rootfind::BisectionConfig;
use laser_solver::rootfind::RootFindConfig::Newton1d;
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
            upper: 100.0,
            num: 20,
        }
    }
}

impl LaserApp {
    pub fn residual_plot(&mut self, ui: &mut Ui) {
        let inputs = linspace(
            self.residual_range.lower,
            self.residual_range.upper,
            self.residual_range.num,
        );

        let kappas = self.grating.grid(self.grid_points.0);
        let dz = self.grid_points.dz(self.fibre_params.length);
        let trial = |sgnl_b| FieldState {
            sgnl_f: 0.0,
            sgnl_b: sgnl_b,
            pump_f: self.pump.total.sqrt(),
            pump_b: 0.0,
        }; //todo: use picard for backward pump
        let f = |sgnl_b| out_field(trial(sgnl_b), self.fibre_params, dz, &kappas).sgnl_b / sgnl_b;
        let residuals: Vec<[f64; 2]> = inputs.iter().map(|&s| [s, f(s).abs().log10()]).collect();

        let mut plt = Plotter::new();
        plt.xlabel("Input");
        plt.ylabel("log(|residual|)");
        plt.add_points(residuals);
        plt.show(ui, "residual-plot");
    }
}

pub fn residual_slider_grid(tr: &mut ResidualRange, ui: &mut Ui) {
    egui::Grid::new("residual").show(ui, |ui| {
        ui.label("thresh up");
        ui.add(egui::Slider::new(&mut tr.upper, 1e-5..=100.0).step_by(0.01));
        ui.end_row();

        ui.label("thresh low");
        ui.add(egui::Slider::new(&mut tr.lower, 1e-6..=tr.upper).step_by(0.01));
        ui.end_row();

        ui.label("thresh num");
        ui.add(egui::Slider::new(&mut tr.num, 5..=100).step_by(0.01));
        ui.end_row();
    });
}
