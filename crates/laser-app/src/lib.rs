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

#[derive(PartialEq, Default, Copy, Clone)]
pub enum View {
    #[default]
    Profile,
    Residual,
    Populations,
    Kappa,
    Threshold,
    PiPosition,
}

const VIEW_OPTIONS: [(View, &str, egui::Key); 6] = [
    (View::Profile, "[1] Profile", egui::Key::Num1),
    (View::Residual, "[2] Residual", egui::Key::Num2),
    (View::Populations, "[3] Populations", egui::Key::Num3),
    (View::Kappa, "[4] Kappa", egui::Key::Num4),
    (View::Threshold, "[5] Threshold", egui::Key::Num5),
    (View::PiPosition, "[6] Pi position", egui::Key::Num6),
];

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

        if !ui.ctx().wants_keyboard_input() {
            for &(view, _, key) in &VIEW_OPTIONS {
                let shortcut = egui::KeyboardShortcut::new(egui::Modifiers::NONE, key);

                if ui.input_mut(|input| input.consume_shortcut(&shortcut)) {
                    changed |= *self != view;
                    *self = view;
                    break;
                }
            }
        }

        ui.horizontal(|ui| {
            for &(view, label, _) in &VIEW_OPTIONS {
                changed |= ui.selectable_value(self, view, label).changed();
            }
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

    fn clear_physics() -> Self {
        Self {
            pump: PumpParam {
                total: 10.0,
                balance: 1.0,
            },
            fibre_params: FibreParams {
                density: 0.50,
                lifetime: 1.0,
                pump_ab: 1.0,
                pump_em: 0.0,
                sgnl_ab: 0.0,
                sgnl_em: 1.0,
                length: 5.0,
            },
            grid_points: GridPoints::default(),
            grating: GratingProfile {
                kappa_left: 0.6,
                kappa_right: 0.6,
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
    fn reset_params(&mut self) {
        *self = Self {
            view: self.view,
            ..Self::clear_physics()
        };
    }

    fn reset_button(&mut self, ui: &mut Ui) {
        let shortcut = egui::KeyboardShortcut::new(egui::Modifiers::NONE, egui::Key::R);

        let shortcut_pressed = !ui.ctx().wants_keyboard_input()
            && ui.input_mut(|input| input.consume_shortcut(&shortcut));

        let shortcut_text = ui.ctx().format_shortcut(&shortcut);

        let button_clicked = ui
            .add(egui::Button::new("Reset").shortcut_text(shortcut_text))
            .clicked();

        if button_clicked || shortcut_pressed {
            self.reset_params();
        }
    }

    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::clear_physics()
    }

    pub fn draw_view_selector(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;

        ui.horizontal(|ui| {
            ui.heading("View: ");
            changed |= self.view.selectors(ui);

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                self.reset_button(ui);
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
            egui::ScrollArea::both().show(ui, |ui| {
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
        });
    }
}

fn bisection_slider_grid(config: &mut BisectionConfig, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("bisection").show(ui, |ui| {
        ui.label("Iteration no.");
        changed |= ui
            .add(egui::Slider::new(&mut config.iteration.max, 10..=2000).step_by(10.0))
            .changed();
        ui.end_row();

        ui.label("Tolerance");
        changed |= ui
            .add(
                egui::Slider::new(&mut config.iteration.tol, 1e-9..=1e-2)
                    // 1e-10 causes slider field box to resize
                    .logarithmic(true)
                    .custom_formatter(|value, _| format!("{value:.1e}")),
            )
            .changed();
        ui.end_row();

        ui.label("Lower bracket");
        changed |= ui
            .add(
                egui::Slider::new(&mut config.lower, 1e-9..=1e-2)
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
        ui.label("Nz");
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
        ui.label("Kappa left");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_left, 0.0..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Kappa right");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_right, 0.0..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Pi Shift Pos");
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
        ui.label("Total power");
        changed |= ui
            .add(egui::Slider::new(&mut pump.total, 0.0..=100.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Balance");
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
            ui.label("Pump em.");
            changed |= ui
                .add(egui::Slider::new(&mut params.pump_em, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Pump abs.");
            changed |= ui
                .add(egui::Slider::new(&mut params.pump_ab, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Signl em.");
            changed |= ui
                .add(egui::Slider::new(&mut params.sgnl_em, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Signl abs.");
            changed |= ui
                .add(egui::Slider::new(&mut params.sgnl_ab, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();
        });

        egui::Grid::new("params2").show(ui, |ui| {
            ui.label("Dopant density");
            changed |= ui
                .add(egui::Slider::new(&mut params.density, 0.1..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Lifetime");
            changed |= ui
                .add(egui::Slider::new(&mut params.lifetime, 0.1..=2.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Length");
            changed |= ui
                .add(egui::Slider::new(&mut params.length, 0.1..=15.0).step_by(0.01))
                .changed();
            ui.end_row();
        });
    });

    changed
}
