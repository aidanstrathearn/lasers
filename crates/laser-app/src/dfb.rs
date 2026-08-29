use crate::ModeUi;
use crate::controls::{
    bisection_slider_grid, fibre_params_slider_grid, grating_slider_grid, pump_slider_grid,
    steps_slider,
};
use crate::plotter::Plotter;
use crate::residual_plot::{ResidualRange, residual_slider_grid};
use crate::threshold_plot::{ThresholdRange, threshold_slider_grid};
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::dfb::{DfbLaser, DfbSolveConfig, Grating};
use laser_solver::error::SolverError;
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, Pump, ResolvedFibre, TwoLevelCrossSections, TwoLevelDopant,
};
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::{BisectionConfig, RootFindConfig};
use std::time::Duration;

#[derive(PartialEq, Default, Copy, Clone)]
pub(crate) enum DfbView {
    #[default]
    Profile,
    Residual,
    Populations,
    Kappa,
    Threshold,
    PiPosition,
}

const VIEW_OPTIONS: [(DfbView, &str, egui::Key); 6] = [
    (DfbView::Profile, "[1] Profile", egui::Key::Num1),
    (DfbView::Populations, "[2] Populations", egui::Key::Num2),
    (DfbView::Kappa, "[3] Kappa", egui::Key::Num3),
    (DfbView::Threshold, "[4] Threshold", egui::Key::Num4),
    (DfbView::PiPosition, "[5] Pi position", egui::Key::Num5),
    (DfbView::Residual, "[6] Residual", egui::Key::Num6),
];

impl DfbView {
    fn plot_id(self) -> &'static str {
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

pub(crate) struct DfbMode {
    pub(crate) view: DfbView,
    pub(crate) pump: Pump,
    pub(crate) fibre: Fibre,
    pub(crate) pump_mode: FieldMode,
    pub(crate) sgnl_mode: FieldMode,
    pub(crate) pump_interaction: TwoLevelCrossSections,
    pub(crate) signal_interaction: TwoLevelCrossSections,
    pub(crate) steps: usize,
    pub(crate) grating: Grating,
    pub(crate) config: BisectionConfig,
    pub(crate) picard_config: PicardConfig,
    pub(crate) threshold_range: ThresholdRange,
    pub(crate) residual_range: ResidualRange,
    cached_plotter: Option<Result<Plotter, SolverError>>,
    pub(crate) compute_time: Option<Duration>,
}

impl Default for DfbMode {
    fn default() -> Self {
        Self {
            view: DfbView::default(),
            pump: Pump {
                total: 10.0,
                balance: 1.0,
            },
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
            grating: Grating {
                kappa_left: 0.6,
                kappa_right: 0.6,
                pi_shift_position: 0.5,
            },
            config: BisectionConfig::default(),
            picard_config: PicardConfig {
                max_iterations: 5_000,
                relative_tolerance: 1e-6,
                absolute_tolerance: 1e-10,
            },
            threshold_range: ThresholdRange::default(),
            residual_range: ResidualRange::default(),
            cached_plotter: None,
            compute_time: None,
        }
    }
}

impl DfbMode {
    pub(crate) fn resolved_fibre(&self) -> ResolvedFibre<'_> {
        self.fibre.resolve_with_interactions(
            self.pump_mode,
            self.pump_interaction,
            self.sgnl_mode,
            self.signal_interaction,
        )
    }

    pub(crate) fn dfb_laser(&self) -> DfbLaser<'_> {
        DfbLaser {
            fibre: self.resolved_fibre(),
            grating: self.grating,
        }
    }

    pub(crate) fn dfb_solve_config(&self, root_find: impl Into<RootFindConfig>) -> DfbSolveConfig {
        DfbSolveConfig {
            steps: self.steps,
            root_find: root_find.into(),
            picard: self.picard_config,
        }
    }

    fn compute_plot(&mut self) -> Result<Plotter, SolverError> {
        match self.view {
            DfbView::Threshold => self.threshold_plot(),
            DfbView::Profile => self.profile_plot(),
            DfbView::Residual => self.residual_plot(),
            DfbView::Populations => self.pops_plot(),
            DfbView::Kappa => self.kappa_plot(),
            DfbView::PiPosition => self.pi_pos_plot(),
        }
    }
}

impl ModeUi for DfbMode {
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
                ui.heading("Bragg");
                changed |= grating_slider_grid(&mut self.grating, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Pump");
                changed |= pump_slider_grid(&mut self.pump, ui);
            });
            ui.vertical(|ui| {
                ui.heading("Solver");
                changed |= bisection_slider_grid(&mut self.config, ui);
                changed |= steps_slider(&mut self.steps, ui);
            });

            match self.view {
                DfbView::Threshold => {
                    ui.vertical(|ui| {
                        ui.heading("Threshold");
                        changed |= threshold_slider_grid(&mut self.threshold_range, ui);
                    });
                }
                DfbView::Residual => {
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
