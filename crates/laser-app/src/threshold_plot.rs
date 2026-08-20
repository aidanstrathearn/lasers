use crate::plotter::Plotter;
use crate::{LaserApp, Points};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::dfb::dfb_pump_scan;
use laser_solver::picard::PicardConfig;
use laser_solver::rootfind::BisectionConfig;
use laser_solver::utils::linspace;

#[derive(Copy, Clone)]
pub struct ThresholdRange {
    lower: f64,
    upper: f64,
    num: usize,
}

impl Default for ThresholdRange {
    fn default() -> Self {
        Self {
            lower: 1e-6,
            upper: 100.0,
            num: 20,
        }
    }
}

impl LaserApp {
    pub fn threshold_plot(&mut self, ui: &mut Ui) {
        let bc = BisectionConfig {
            upper: self.threshold_range.upper.sqrt(),
            ..self.config
        };

        let picard_config = PicardConfig {
            max_iterations: 5_000,
            relative_tolerance: 1e-6,
            absolute_tolerance: 1e-10,
        };

        let pumps = linspace(
            self.threshold_range.lower,
            self.threshold_range.upper,
            self.threshold_range.num,
        );
        let threshold = dfb_pump_scan(
            &pumps,
            self.pump.balance,
            self.fibre_params,
            self.grid_points,
            self.grating,
            bc,
            picard_config,
        );
        let sgnl_f = threshold.iter().map(|x| x.0.powi(2));
        let sgnl_b = threshold.iter().map(|x| x.1.powi(2));

        let sgnl_f_points: Points = pumps
            .iter()
            .zip(sgnl_f)
            .map(|(&x, y)| [x, y.abs()])
            .collect();
        let sgnl_b_points: Points = pumps
            .iter()
            .zip(sgnl_b)
            .map(|(&x, y)| [x, y.abs()])
            .collect();

        let mut plt = Plotter::new();
        plt.xlabel("Pump");
        plt.ylabel("Signal");
        plt.add_points(sgnl_f_points).label("Forward");
        plt.add_points(sgnl_b_points).label("Backward");
        plt.show(ui, "threshold-plot");
    }
}

pub fn threshold_slider_grid(tr: &mut ThresholdRange, ui: &mut Ui) {
    egui::Grid::new("threshold").show(ui, |ui| {
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
