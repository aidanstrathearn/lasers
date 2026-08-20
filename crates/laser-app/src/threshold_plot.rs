use eframe::egui;
use eframe::egui::Ui;
use egui_plot::{Legend, Line, Plot};
use laser_solver::dfb::dfb_pump_scan;
use laser_solver::rootfind::BisectionConfig;
use laser_solver::utils::linspace;
use crate::{LaserApp, Points};


#[derive(Copy, Clone)]
pub struct ThresholdRange {
    lower: f64,
    upper: f64
}

impl Default for ThresholdRange {
    fn default() -> Self {
        Self{ lower: 1e-6, upper: 10.0}
    }
}

impl LaserApp {
    pub fn threshold_plot(&mut self, ui: &mut Ui) {
        let bc = BisectionConfig {
            upper: self.threshold_range.upper,
            ..self.config
        };

        let pumps = linspace(self.threshold_range.lower, self.threshold_range.upper, 20);
        let threshold = dfb_pump_scan(
            &pumps,
            self.fibre_params,
            self.grid_points,
            self.grating,
            bc,
        );
        let sgnl_f = threshold.iter().map(|x| x.0);
        let sgnl_b = threshold.iter().map(|x| x.1);
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
        Plot::new("threshold")
            .legend(Legend::default())
            .x_axis_label("pump")
            .y_axis_label("signal")
            .show(ui, |plot_ui| {
                plot_ui.line(
                    Line::new("threshold f", sgnl_f_points)
                        .name("Forward Signal")
                        .width(3.0),
                );
                plot_ui.line(
                    Line::new("threshold b", sgnl_b_points)
                        .name("Backward Signal")
                        .width(3.0),
                );
            });

    }
}

pub fn threshold_slider_grid(tr: &mut ThresholdRange, ui: &mut Ui) {
    egui::Grid::new("threshold").show(ui, |ui| {
        ui.label("thresh up");
        ui.add(egui::Slider::new(&mut tr.upper, 1e-5..=15.0).step_by(0.01));
        ui.end_row();

        ui.label("thresh low");
        ui.add(egui::Slider::new(&mut tr.lower, 1e-6..=tr.upper).step_by(0.01));
        ui.end_row();
    });
}