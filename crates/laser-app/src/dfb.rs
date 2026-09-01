use crate::controls::{
    bisection_slider_grid, fibre_params_slider_grid, grating_slider_grid, pump_slider_grid,
    steps_slider,
};
use crate::plotter::Plotter;
use crate::threshold_plot::{ThresholdRange, threshold_slider_grid};
use crate::{LaserParameters, ModeUi};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::BisectionConfig;
use std::time::Duration;

#[derive(PartialEq, Default, Copy, Clone)]
pub(crate) enum DfbView {
    #[default]
    Profile,
    Populations,
    Kappa,
    Threshold,
    PiPosition,
    PiPositionThreshold,
}

const VIEW_OPTIONS: [(DfbView, &str, egui::Key); 6] = [
    (DfbView::Profile, "[1] Profile", egui::Key::Num1),
    (DfbView::Populations, "[2] Populations", egui::Key::Num2),
    (DfbView::Kappa, "[3] Kappa", egui::Key::Num3),
    (DfbView::Threshold, "[4] Threshold", egui::Key::Num4),
    (DfbView::PiPosition, "[5] Pi position", egui::Key::Num5),
    (
        DfbView::PiPositionThreshold,
        "[6] Pi threshold",
        egui::Key::Num6,
    ),
];

impl DfbView {
    fn plot_id(self) -> &'static str {
        match self {
            Self::Threshold => "threshold-plot",
            Self::Profile => "profile-plot",
            Self::Populations => "population-plot",
            Self::Kappa => "kappa-plot",
            Self::PiPosition => "pi-position-output-plot",
            Self::PiPositionThreshold => "pi-position-threshold-plot",
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

pub(crate) struct DfbMode {
    pub(crate) view: DfbView,
    pub(crate) steps: usize,
    pub(crate) config: BisectionConfig,
    pub(crate) picard_config: PicardConfig,
    pub(crate) threshold_range: ThresholdRange,
    cached_plotter: Option<Result<Plotter, SolverError>>,
    pub(crate) compute_time: Option<Duration>,
}

impl DfbMode {
    pub(crate) fn new(parameters: &LaserParameters) -> Self {
        let fibre = parameters.resolved_fibre();
        Self {
            view: DfbView::default(),
            steps: 100,
            config: BisectionConfig::default(),
            picard_config: PicardConfig {
                max_iterations: 5_000,
                relative_tolerance: 1e-6,
                absolute_tolerance: 1e-10,
            },
            threshold_range: ThresholdRange::new_watts(
                fibre.pump_power(1e-6),
                fibre.pump_power(10.0),
                20,
            ),
            cached_plotter: None,
            compute_time: None,
        }
    }
}

impl DfbMode {
    fn compute_plot(&mut self, parameters: &LaserParameters) -> Result<Plotter, SolverError> {
        match self.view {
            DfbView::Threshold => self.threshold_plot(parameters),
            DfbView::Profile => self.profile_plot(parameters),
            DfbView::Populations => self.pops_plot(parameters),
            DfbView::Kappa => self.kappa_plot(parameters),
            DfbView::PiPosition => self.pi_pos_plot(parameters),
            DfbView::PiPositionThreshold => self.pi_pos_threshold_plot(parameters),
        }
    }
}

impl ModeUi for DfbMode {
    fn draw_view_selector(&mut self, ui: &mut Ui) -> bool {
        self.view.selectors(ui)
    }

    fn draw_controls(&mut self, parameters: &mut LaserParameters, ui: &mut Ui) -> bool {
        let mut changed = false;

        egui::Grid::new("global-params").show(ui, |ui| {
            ui.vertical(|ui| {
                ui.heading("Fibre");
                changed |= fibre_params_slider_grid(
                    &mut parameters.fibre,
                    &mut parameters.pump_interaction,
                    &mut parameters.signal_interaction,
                    ui,
                );
            });
            ui.vertical(|ui| {
                ui.heading("Bragg");
                changed |= grating_slider_grid(&mut parameters.fibre.grating, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Pump");
                changed |= pump_slider_grid(&mut parameters.pump, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Solver");
                changed |= bisection_slider_grid(&mut self.config, ui);
                changed |= steps_slider(&mut self.steps, ui);
            });

            match self.view {
                DfbView::Threshold | DfbView::PiPositionThreshold => {
                    ui.vertical(|ui| {
                        ui.heading("Threshold");
                        changed |= threshold_slider_grid(&mut self.threshold_range, ui);
                    });
                }
                _ => (),
            };

            ui.end_row();
        });

        changed
    }

    fn reset(&mut self, parameters: &LaserParameters) {
        *self = Self {
            view: self.view,
            ..Self::new(parameters)
        };
    }

    fn clear_cached_plot(&mut self) {
        self.cached_plotter = None;
    }

    fn has_cached_plot(&self) -> bool {
        self.cached_plotter.is_some()
    }

    fn recompute_plot(&mut self, parameters: &LaserParameters) {
        self.cached_plotter = Some(self.compute_plot(parameters));
    }

    fn compute_time(&self) -> Option<Duration> {
        self.compute_time
    }

    fn draw_plot(&self, ui: &mut Ui) {
        match &self.cached_plotter {
            Some(Ok(plotter)) => plotter.show(ui, self.view.plot_id()),
            Some(Err(error)) => {
                ui.colored_label(ui.visuals().error_fg_color, error.to_string());
            }
            None => {}
        }
    }
}
