use crate::controls::{
    bisection_slider_grid, fibre_params_slider_grid, grating_slider_grid, power_slider_mw,
    pump_slider_grid, steps_slider,
};
use crate::plotter::Plotter;
use crate::{LaserParameters, ModeUi, field_profile_plot, timed};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::error::SolverError;
use laser_solver::lase::Signal;
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::maths::rootfind::Midpoint::Arithmetic;
use laser_solver::two_mode::TwoModeSolver;
use std::time::Duration;

#[derive(PartialEq, Default, Copy, Clone)]
pub(crate) enum AmplifierView {
    #[default]
    Profile,
}

const VIEW_OPTIONS: [(AmplifierView, &str, egui::Key); 1] =
    [(AmplifierView::Profile, "[1] Profile", egui::Key::Num1)];

impl AmplifierView {
    fn plot_id(self) -> &'static str {
        match self {
            Self::Profile => "profile-plot",
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

pub(crate) struct AmplifierMode {
    pub(crate) view: AmplifierView,
    pub(crate) signal: Signal,
    pub(crate) steps: usize,
    pub(crate) config: BisectionConfig,
    cached_plotter: Option<Result<Plotter, SolverError>>,
    pub(crate) compute_time: Option<Duration>,
}

impl AmplifierMode {
    pub(crate) fn new(parameters: &LaserParameters) -> Self {
        let signal_total = parameters.resolved_fibre().signal_power(1.0);
        Self {
            view: AmplifierView::default(),
            signal: Signal {
                total: signal_total,
                ..Signal::default()
            },
            steps: 100,
            config: BisectionConfig::default(),
            cached_plotter: None,
            compute_time: None,
        }
    }
}

pub(crate) fn signal_slider_grid(signal: &mut Signal, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("signal").show(ui, |ui| {
        ui.label("Total power (mW)");
        changed |= power_slider_mw(&mut signal.total, ui);
        ui.end_row();

        ui.label("Balance");
        changed |= ui
            .add(egui::Slider::new(&mut signal.balance, -1.0..=1.0).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}

impl AmplifierMode {
    pub(crate) fn picard_config(&self) -> PicardConfig {
        PicardConfig {
            max_iterations: 5000,
            relative_tolerance: 1e-3, // need higher tolerance than dfb?
            absolute_tolerance: 1e-3,
            ..PicardConfig::default()
        }
    }

    fn compute_plot(&mut self, parameters: &LaserParameters) -> Result<Plotter, SolverError> {
        match self.view {
            AmplifierView::Profile => self.profile_plot(parameters),
        }
    }
}

impl AmplifierMode {
    fn profile_plot(&mut self, parameters: &LaserParameters) -> Result<Plotter, SolverError> {
        let (result, compute_time) = timed(|| {
            let fibre = parameters.resolved_fibre();
            let bc = BisectionConfig {
                upper: fibre.pump_flux(parameters.pump.total).sqrt(),
                lower: 0.0,
                midpoint: Arithmetic,
                ..self.config
            };
            TwoModeSolver::new(&fibre, self.steps).solve_injected(
                parameters.pump,
                self.signal,
                bc.into(),
                self.picard_config(),
            )
        });
        self.compute_time = Some(compute_time);
        Ok(field_profile_plot(&result?))
    }
}

impl ModeUi for AmplifierMode {
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
                ui.heading("Signal");
                changed |= signal_slider_grid(&mut self.signal, ui)
            });
            ui.vertical(|ui| {
                ui.heading("Solver");
                changed |= bisection_slider_grid(&mut self.config, ui);
                changed |= steps_slider(&mut self.steps, ui);
            });

            match self.view {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn injected_mode_uses_shared_pi_shift_grating() {
        let parameters = LaserParameters::default();
        let fibre = parameters.resolved_fibre();
        let solver = TwoModeSolver::new(&fibre, AmplifierMode::new(&parameters).steps);

        assert!(solver.kappas().iter().any(|&kappa| kappa != 0.0));
    }
}
