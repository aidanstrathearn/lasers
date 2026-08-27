use crate::plotter::Plotter;
use crate::{Points, dfb::DfbMode, timed};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
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
            upper: 10.0,
            num: 20,
        }
    }
}

impl DfbMode {
    pub fn threshold_plot(&mut self) -> Result<Plotter, SolverError> {
        let bc = BisectionConfig {
            upper: self.threshold_range.upper.sqrt(),
            ..self.config
        };
        // let nc = Newton1dConfig {
        //
        //     iter
        // }

        let pumps = linspace(
            self.threshold_range.lower,
            self.threshold_range.upper,
            self.threshold_range.num,
        );
        let (threshold, compute_time) = timed(|| {
            self.dfb_laser()
                .pump_scan(&pumps, self.pump.balance, self.dfb_solve_config(bc))
        });
        let threshold = threshold?;
        let sgnl_f = threshold
            .iter()
            .map(|output| output.as_ref().map_or(0.0, |output| output.0));
        let sgnl_b = threshold
            .iter()
            .map(|output| output.as_ref().map_or(0.0, |output| output.1));

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
        plt.xlabel("Pump power");
        plt.ylabel("Output power");
        plt.add_points(sgnl_f_points).label("Forward");
        plt.add_points(sgnl_b_points).label("Backward");
        plt.axvline(self.pump.total).label("Current pump");
        plt.xlim(self.threshold_range.lower, self.threshold_range.upper);
        plt.set_compute_time(compute_time);
        Ok(plt)
    }
}

pub fn threshold_slider_grid(tr: &mut ThresholdRange, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("threshold").show(ui, |ui| {
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
