mod profile_plot;
mod threshold_plot;
mod plotter;

use crate::threshold_plot::{ThresholdRange, threshold_slider_grid};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::lase::{FibreParams, GratingProfile, GridPoints, Pump};
use laser_solver::rootfind::BisectionConfig;

#[derive(PartialEq, Default)]
pub enum View {
    #[default]
    Profile,
    Threshold,
}

type Points = Vec<[f64; 2]>;

#[derive(Default)]
pub struct LaserApp {
    view: View,
    pump: PumpParam,
    fibre_params: FibreParams,
    grid_points: GridPoints,
    grating: GratingProfile,
    config: BisectionConfig,
    threshold_range: ThresholdRange,
}

struct PumpParam {
    total: f64,
    balance: f64
}

impl Default for PumpParam {
    fn default() -> Self {
        Self {total: 100.0,  balance: 1.0}
    }
}

impl LaserApp {
    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::default()
    }
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("View: ");
                view_selectors(&mut self.view, ui);
            });

            ui.separator();

            egui::Grid::new("global-params").show(ui, |ui| {
                ui.vertical(|ui| {
                    ui.heading("Fibre");
                    fibre_params_slider_grid(&mut self.fibre_params, ui);
                });
                ui.vertical(|ui| {
                    ui.heading("Bragg");
                    grating_slider_grid(&mut self.grating, ui);
                });
                ui.vertical(|ui| {
                    ui.heading("Solver");
                    bisection_slider_grid(&mut self.config, ui);
                });
                ui.vertical(|ui| {
                    ui.heading("Pump");
                    pump_param_slider_grid(&mut self.pump, ui);
                });
                ui.vertical(|ui| {
                    ui.heading("Threshold");
                    threshold_slider_grid(&mut self.threshold_range, ui);
                });
                ui.end_row();
            });

            ui.separator();

            match self.view {
                View::Threshold => self.threshold_plot(ui),
                View::Profile => self.profile_plot(ui).unwrap_or_else(|error| {
                    ui.colored_label(ui.visuals().error_fg_color, error.to_string());
                }),
            };
        });
    }
}

fn view_selectors(view: &mut View, ui: &mut Ui) {
    ui.horizontal(|ui| {
        ui.selectable_value(view, View::Profile, "Profile");
        ui.selectable_value(view, View::Threshold, "Threshold");
    });
}

fn bisection_slider_grid(config: &mut BisectionConfig, ui: &mut Ui) {
    egui::Grid::new("bisection").show(ui, |ui| {
        ui.label("iters");
        ui.add(egui::Slider::new(&mut config.iteration.max, 10..=2000).step_by(10.0));
        ui.end_row();

        // uh oh - how to do log slider?
        // ui.label("tolerance");
        // ui.add(egui::Slider::new(&mut config.tolerance, 0.0..=1.0).step_by(0.01));
        // ui.end_row();
    });
}
fn grating_slider_grid(grating: &mut GratingProfile, ui: &mut Ui) {
    egui::Grid::new("grating").show(ui, |ui| {
        ui.label("kappa left");
        ui.add(egui::Slider::new(&mut grating.kappa_left, 0.0..=10.0).step_by(0.01));
        ui.end_row();

        ui.label("kappa right");
        ui.add(egui::Slider::new(&mut grating.kappa_right, 0.0..=10.0).step_by(0.01));
        ui.end_row();

        ui.label("pi pos");
        ui.add(egui::Slider::new(&mut grating.pi_shift_position, 0.0..=1.0).step_by(0.01));
        ui.end_row();
    });
}

fn pump_slider_grid(pump: &mut Pump, ui: &mut Ui) {
    egui::Grid::new("pump").show(ui, |ui| {
        ui.label("pump f");
        ui.add(egui::Slider::new(&mut pump.forward, 0.0..=100.0).step_by(0.01));
        ui.end_row();

        ui.label("pump b");
        ui.add(egui::Slider::new(&mut pump.backward, 0.0..=100.0).step_by(0.01));
        ui.end_row();
    });
}

fn pump_param_slider_grid(pump: &mut PumpParam, ui: &mut Ui) {
    egui::Grid::new("pumpp").show(ui, |ui| {
        ui.label("pump-total");
        ui.add(egui::Slider::new(&mut pump.total, 0.0..=1000.0).step_by(0.01));
        ui.end_row();

        ui.label("pump-balance");
        ui.add(egui::Slider::new(&mut pump.balance, -1.0..=1.0).step_by(0.01));
        ui.end_row();
    });
}

fn fibre_params_slider_grid(params: &mut FibreParams, ui: &mut Ui) {
    egui::Grid::new("params").show(ui, |ui| {
        egui::Grid::new("params1").show(ui, |ui| {
            ui.label("pump_em");
            ui.add(egui::Slider::new(&mut params.pump_em, 0.0..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("pump_ab");
            ui.add(egui::Slider::new(&mut params.pump_ab, 0.0..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("sgnl_em");
            ui.add(egui::Slider::new(&mut params.sgnl_em, 0.0..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("sgnl_ab");
            ui.add(egui::Slider::new(&mut params.sgnl_ab, 0.0..=10.0).step_by(0.01));
            ui.end_row();
        });

        egui::Grid::new("params2").show(ui, |ui| {
            ui.label("density");
            ui.add(egui::Slider::new(&mut params.density, 0.1..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("lifetime");
            ui.add(egui::Slider::new(&mut params.lifetime, 0.1..=2.0).step_by(0.01));
            ui.end_row();

            ui.label("length");
            ui.add(egui::Slider::new(&mut params.length, 0.1..=50.0).step_by(0.01));
            ui.end_row();
        });
    });
}
