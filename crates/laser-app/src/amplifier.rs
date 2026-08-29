use crate::controls::{
    bisection_slider_grid, fibre_params_slider_grid, grating_slider_grid, pump_slider_grid,
    steps_slider,
};
use crate::plotter::Plotter;
use crate::residual_plot::{ResidualRange, residual_slider_grid};
use crate::threshold_plot::{ThresholdRange, threshold_slider_grid};
use crate::{ModeUi, field_profile_plot, timed};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::amplifier::{Amplifier, AmplifierSolveConfig};
use laser_solver::dfb::{DfbLaser, DfbSolveConfig, Grating};
use laser_solver::error::SolverError;
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, Pump, Signal, TwoLevelCrossSections, TwoLevelDopant,
};
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::Midpoint::Arithmetic;
use laser_solver::maths::rootfind::{BisectionConfig, Midpoint, RootFindConfig};
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
    pub(crate) pump: Pump,
    pub(crate) fibre: Fibre,
    pub(crate) pump_mode: FieldMode,
    pub(crate) sgnl_mode: FieldMode,
    pub(crate) pump_interaction: TwoLevelCrossSections,
    pub(crate) signal_interaction: TwoLevelCrossSections,
    pub(crate) steps: usize,
    pub(crate) config: BisectionConfig,
    cached_plotter: Option<Result<Plotter, SolverError>>,
    pub(crate) compute_time: Option<Duration>,
}

impl Default for AmplifierMode {
    fn default() -> Self {
        Self {
            view: AmplifierView::default(),
            pump: Pump {
                total: 10.0,
                balance: 1.0,
            },
            signal: Signal::default(),
            fibre: Fibre {
                geometry: FibreGeometry {
                    core_radius: 4e-6,
                    numerical_aperture: 0.1,
                    length: 5.0,
                },
                dopant: TwoLevelDopant {
                    density: 0.50,
                    lifetime: 1.0,
                },
            },
            pump_mode: FieldMode::new(970e-9),
            sgnl_mode: FieldMode::new(1060e-9),
            pump_interaction: TwoLevelCrossSections::new(1.0, 0.0),
            signal_interaction: TwoLevelCrossSections::new(0.0, 1.0),
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
        ui.label("Total power");
        changed |= ui
            .add(egui::Slider::new(&mut signal.total, 0.0..=100.0).step_by(0.01))
            .changed();
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
    pub(crate) fn amplifier(&self) -> Amplifier<'_> {
        Amplifier {
            fibre: self.fibre.resolve_with_interactions(
                self.pump_mode,
                self.pump_interaction,
                self.sgnl_mode,
                self.signal_interaction,
            ),
        }
    }

    pub(crate) fn amplifier_solve_config(
        &self,
        root_find: impl Into<RootFindConfig>,
    ) -> AmplifierSolveConfig {
        AmplifierSolveConfig {
            steps: self.steps,
            root_find: root_find.into(),
            picard: PicardConfig {
                max_iterations: 2000,
                relative_tolerance: 1e-3, // need higher tolerance than dfb?
                absolute_tolerance: 1e-3,
                ..PicardConfig::default()
            },
        }
    }

    fn compute_plot(&mut self) -> Result<Plotter, SolverError> {
        match self.view {
            AmplifierView::Profile => self.profile_plot(),
        }
    }
}

impl AmplifierMode {
    fn profile_plot(&mut self) -> Result<Plotter, SolverError> {
        let full_profile = true;
        let bc = BisectionConfig {
            upper: self.pump.total.sqrt(),
            lower: 0.0,
            midpoint: Arithmetic,
            ..self.config
        };

        let (result, compute_time) = timed(|| {
            self.amplifier().solve(
                self.signal,
                self.pump,
                self.amplifier_solve_config(bc),
                full_profile,
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

    fn draw_controls(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;

        egui::Grid::new("global-params").show(ui, |ui| {
            ui.vertical(|ui| {
                ui.heading("Fibre");
                changed |= fibre_params_slider_grid(
                    &mut self.fibre,
                    &mut self.pump_interaction,
                    &mut self.signal_interaction,
                    ui,
                );
            });
            ui.vertical(|ui| {
                ui.heading("Pump");
                changed |= pump_slider_grid(&mut self.pump, ui);
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

    fn reset(&mut self) {
        *self = Self {
            view: self.view,
            ..Self::default()
        };
    }

    fn has_cached_plot(&self) -> bool {
        self.cached_plotter.is_some()
    }

    fn recompute_plot(&mut self) {
        self.cached_plotter = Some(self.compute_plot());
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
