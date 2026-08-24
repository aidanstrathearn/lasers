mod kappa_plot;
mod pi_pos_plot;
mod plotter;
mod pop_plot;
mod profile_plot;
mod residual_plot;
mod threshold_plot;

use crate::plotter::Plotter;
use crate::residual_plot::{ResidualRange, residual_slider_grid};
use crate::threshold_plot::{ThresholdRange, threshold_slider_grid};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
use laser_solver::lase::{FibreParams, GratingProfile, GridPoints, Pump, PumpParam};
use laser_solver::picard::PicardConfig;
use laser_solver::rootfind::BisectionConfig;
use std::time::Duration;
use web_time::Instant;

#[derive(PartialEq, Default)]
pub enum View {
    #[default]
    Profile,
    Residual,
    Populations,
    Kappa,
    Threshold,
    PiPosition,
}

impl View {
    fn plot_id(&self) -> &'static str {
        match self {
            Self::Threshold => "threshold-plot",
            Self::Profile => "profile-plot",
            Self::Residual => "residual-plot",
            Self::Populations => "population-plot",
            Self::Kappa => "kappa-plot",
            Self::PiPosition => "pi-position-output-plot",
        }
    }

    fn selectors(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;

        ui.horizontal(|ui| {
            changed |= ui
                .selectable_value(self, Self::Profile, "Profile")
                .changed();
            changed |= ui
                .selectable_value(self, Self::Residual, "Residual")
                .changed();
            changed |= ui
                .selectable_value(self, Self::Populations, "Populations")
                .changed();
            changed |= ui.selectable_value(self, Self::Kappa, "Kappa").changed();
            changed |= ui
                .selectable_value(self, Self::Threshold, "Threshold")
                .changed();
            changed |= ui
                .selectable_value(self, Self::PiPosition, "Pi position")
                .changed();
        });

        changed
    }
}

type Points = Vec<[f64; 2]>;

fn timed<T>(compute: impl FnOnce() -> T) -> (T, Duration) {
    let start = Instant::now();
    let result = compute();
    (result, start.elapsed())
}

#[derive(Default)]
pub struct LaserApp {
    view: View,
    pump: PumpParam,
    fibre_params: FibreParams,
    grid_points: GridPoints,
    grating: GratingProfile,
    config: BisectionConfig,
    picard_config: PicardConfig,
    threshold_range: ThresholdRange,
    residual_range: ResidualRange,
    cached_plotter: Option<Result<Plotter, SolverError>>,
}

impl LaserApp {
    fn strong_coupling() -> Self {
        Self {
            pump: PumpParam {
                total: 50.0,
                balance: 1.0,
            },
            fibre_params: FibreParams {
                density: 0.20,
                lifetime: 1.0,
                pump_ab: 0.1,
                pump_em: 0.0,
                sgnl_ab: 0.1,
                sgnl_em: 0.1,
                length: 5.0,
            },
            grid_points: GridPoints::default(),
            grating: GratingProfile {
                kappa_left: 3.0,
                kappa_right: 3.0,
                pi_shift_position: 0.5,
            },
            picard_config: PicardConfig {
                max_iterations: 5_000,
                relative_tolerance: 1e-6,
                absolute_tolerance: 1e-10,
            },
            ..Self::default()
        }
    }

    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::strong_coupling()
    }

    pub fn draw_view_selector(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;

        ui.horizontal(|ui| {
            ui.heading("View: ");
            changed |= self.view.selectors(ui);

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if let Some(Ok(plotter)) = &self.cached_plotter {
                    let milliseconds = plotter.compute_time().as_secs_f64() * 1_000.0;
                    ui.label(format!("Compute: {milliseconds:.3} ms"));
                }
            });
        });

        changed
    }

    pub fn compute_plot(&mut self) -> Result<Plotter, SolverError> {
        match self.view {
            View::Threshold => self.threshold_plot(),
            View::Profile => self.profile_plot(),
            View::Residual => self.residual_plot(),
            View::Populations => self.pops_plot(),
            View::Kappa => self.kappa_plot(),
            View::PiPosition => self.pi_pos_plot(),
        }
    }

    pub fn draw_plot(&self, ui: &mut Ui) {
        match &self.cached_plotter {
            Some(Ok(plotter)) => plotter.show(ui, self.view.plot_id()),
            Some(Err(error)) => {
                ui.colored_label(ui.visuals().error_fg_color, error.to_string());
            }
            None => {}
        }
    }

    pub fn draw_controls(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;

        egui::Grid::new("global-params").show(ui, |ui| {
            ui.vertical(|ui| {
                ui.heading("Fibre");
                changed |= fibre_params_slider_grid(&mut self.fibre_params, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Bragg");
                changed |= grating_slider_grid(&mut self.grating, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Pump");
                changed |= pump_param_slider_grid(&mut self.pump, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Solver");
                changed |= bisection_slider_grid(&mut self.config, ui);
                changed |= gridpoints_slider(&mut self.grid_points, ui);
            });

            match self.view {
                View::Threshold => {
                    ui.vertical(|ui| {
                        ui.heading("Threshold");
                        changed |= threshold_slider_grid(&mut self.threshold_range, ui);
                    });
                }

                View::Residual => {
                    ui.vertical(|ui| {
                        ui.heading("Residual");
                        changed |= residual_slider_grid(&mut self.residual_range, ui);
                    });
                }
                _ => (),
            };

            ui.end_row();
        });

        changed
    }
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let mut changed = false;
        egui::CentralPanel::default().show(ctx, |ui| {
            changed |= self.draw_view_selector(ui);

            ui.separator();

            changed |= self.draw_controls(ui);
            if changed || self.cached_plotter.is_none() {
                self.cached_plotter = Some(self.compute_plot());
                ctx.request_repaint();
            }
            ui.separator();

            self.draw_plot(ui);
        });
    }
}

fn bisection_slider_grid(config: &mut BisectionConfig, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("bisection").show(ui, |ui| {
        ui.label("iters");
        changed |= ui
            .add(egui::Slider::new(&mut config.iteration.max, 10..=2000).step_by(10.0))
            .changed();
        ui.end_row();

        ui.label("tolerance");
        changed |= ui
            .add(
                egui::Slider::new(&mut config.iteration.tol, 1e-9..=1e-2)
                    // 1e-10 causes slider field box to resize
                    .logarithmic(true)
                    .custom_formatter(|value, _| format!("{value:.1e}")),
            )
            .changed();
        ui.end_row();
    });

    changed
}

fn gridpoints_slider(gp: &mut GridPoints, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("grid-points").show(ui, |ui| {
        ui.label("grid points");
        changed |= ui
            .add(egui::Slider::new(&mut gp.0, 10..=1000).step_by(2.0))
            .changed();
        ui.end_row();
    });

    changed
}

fn grating_slider_grid(grating: &mut GratingProfile, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("grating").show(ui, |ui| {
        ui.label("kappa left");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_left, 0.0..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("kappa right");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_right, 0.0..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("pi pos");
        changed |= ui
            .add(egui::Slider::new(&mut grating.pi_shift_position, 0.0..=1.0).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}

fn pump_param_slider_grid(pump: &mut PumpParam, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("pumpp").show(ui, |ui| {
        ui.label("pump-total");
        changed |= ui
            .add(egui::Slider::new(&mut pump.total, 0.0..=100.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("pump-balance");
        changed |= ui
            .add(egui::Slider::new(&mut pump.balance, -1.0..=1.0).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}

fn fibre_params_slider_grid(params: &mut FibreParams, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("params").show(ui, |ui| {
        egui::Grid::new("params1").show(ui, |ui| {
            ui.label("pump_em");
            changed |= ui
                .add(egui::Slider::new(&mut params.pump_em, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("pump_ab");
            changed |= ui
                .add(egui::Slider::new(&mut params.pump_ab, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("sgnl_em");
            changed |= ui
                .add(egui::Slider::new(&mut params.sgnl_em, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("sgnl_ab");
            changed |= ui
                .add(egui::Slider::new(&mut params.sgnl_ab, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();
        });

        egui::Grid::new("params2").show(ui, |ui| {
            ui.label("density");
            changed |= ui
                .add(egui::Slider::new(&mut params.density, 0.1..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("lifetime");
            changed |= ui
                .add(egui::Slider::new(&mut params.lifetime, 0.1..=2.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("length");
            changed |= ui
                .add(egui::Slider::new(&mut params.length, 0.1..=15.0).step_by(0.01))
                .changed();
            ui.end_row();
        });
    });

    changed
}
