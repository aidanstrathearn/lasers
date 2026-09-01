mod amplifier;
mod controls;
mod dfb;
mod kappa_plot;
mod pi_pos_plot;
mod pi_pos_threshold_plot;
mod plotter;
mod pop_plot;
mod profile_plot;
mod threshold_plot;

use crate::amplifier::AmplifierMode;
use crate::dfb::DfbMode;
use crate::plotter::Plotter;
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::grating::PiShift;
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, FieldProfile, Pump, ResolvedFibre, TwoLevelCrossSections,
    TwoLevelDopant,
};
use std::time::Duration;
use web_time::Instant;

type Points = Vec<[f64; 2]>;

fn timed<T>(compute: impl FnOnce() -> T) -> (T, Duration) {
    let start = Instant::now();
    let result = compute();
    (result, start.elapsed())
}

fn power_points(z: impl Iterator<Item = f64>, power_watts: impl Iterator<Item = f64>) -> Points {
    z.zip(power_watts)
        .map(|(z, power_watts)| [z, 1_000.0 * power_watts])
        .collect()
}

fn field_profile_plot(profile: &FieldProfile) -> Plotter {
    let mut plot = Plotter::new();
    plot.add_points(power_points(profile.z(), profile.signal_forward_power()))
        .label("Forward signal");
    plot.add_points(power_points(profile.z(), profile.signal_backward_power()))
        .label("Backward signal");
    plot.add_points(power_points(profile.z(), profile.pump_forward_power()))
        .label("Forward pump");
    plot.add_points(power_points(profile.z(), profile.pump_backward_power()))
        .label("Backward pump");
    plot.xlabel("Position (m)");
    plot.ylabel("Power (mW)");
    plot
}

trait ModeUi {
    fn draw_view_selector(&mut self, ui: &mut Ui) -> bool;
    fn draw_controls(&mut self, parameters: &mut LaserParameters, ui: &mut Ui) -> bool;
    fn reset(&mut self, parameters: &LaserParameters);
    fn clear_cached_plot(&mut self);
    fn has_cached_plot(&self) -> bool;
    fn recompute_plot(&mut self, parameters: &LaserParameters);
    fn compute_time(&self) -> Option<Duration>;
    fn draw_plot(&self, ui: &mut Ui);
}

pub(crate) struct LaserParameters {
    pub(crate) fibre: Fibre<TwoLevelDopant, PiShift>,
    pub(crate) pump_mode: FieldMode,
    pub(crate) sgnl_mode: FieldMode,
    pub(crate) pump_interaction: TwoLevelCrossSections,
    pub(crate) signal_interaction: TwoLevelCrossSections,
    pub(crate) pump: Pump,
}

impl LaserParameters {
    pub(crate) fn resolved_fibre(&self) -> ResolvedFibre<'_, TwoLevelDopant, PiShift> {
        self.fibre.resolve_with_interactions(
            self.pump_mode,
            self.pump_interaction,
            self.sgnl_mode,
            self.signal_interaction,
        )
    }
}

impl Default for LaserParameters {
    fn default() -> Self {
        let fibre = Fibre {
            geometry: FibreGeometry {
                core_radius: 4e-6,
                numerical_aperture: 0.1,
                length: 5.0,
            },
            dopant: TwoLevelDopant {
                density: 0.50,
                lifetime: 1.0,
            },
            grating: PiShift {
                kappa_left: 0.6,
                kappa_right: 0.6,
                pi_shift_position: 0.5,
            },
        };
        let pump_mode = FieldMode::new(970e-9);
        let sgnl_mode = FieldMode::new(1060e-9);
        let pump_interaction = TwoLevelCrossSections::new(1.0, 0.0);
        let signal_interaction = TwoLevelCrossSections::new(0.0, 1.0);
        let pump_total = fibre
            .resolve_with_interactions(pump_mode, pump_interaction, sgnl_mode, signal_interaction)
            .pump_power(10.0);

        Self {
            fibre,
            pump_mode,
            sgnl_mode,
            pump_interaction,
            signal_interaction,
            pump: Pump {
                total: pump_total,
                balance: 1.0,
            },
        }
    }
}

#[derive(Default, PartialEq, Copy, Clone)]
enum Mode {
    #[default]
    Dfb,
    Amplifier,
}

impl Mode {
    fn selectors(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;
        changed |= ui.selectable_value(self, Self::Dfb, "DFB").changed();
        changed |= ui
            .selectable_value(self, Self::Amplifier, "Amplifier")
            .changed();
        changed
    }
}

pub struct LaserApp {
    selected_mode: Mode,
    parameters: LaserParameters,
    dfb: DfbMode,
    amplifier: AmplifierMode,
}

impl Default for LaserApp {
    fn default() -> Self {
        let parameters = LaserParameters::default();
        Self {
            selected_mode: Mode::default(),
            dfb: DfbMode::new(&parameters),
            amplifier: AmplifierMode::new(&parameters),
            parameters,
        }
    }
}

impl LaserApp {
    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::default()
    }

    fn active_mode_mut(&mut self) -> &mut dyn ModeUi {
        match self.selected_mode {
            Mode::Dfb => &mut self.dfb,
            Mode::Amplifier => &mut self.amplifier,
        }
    }

    fn active_mode_and_parameters_mut(&mut self) -> (&mut LaserParameters, &mut dyn ModeUi) {
        match self.selected_mode {
            Mode::Dfb => (&mut self.parameters, &mut self.dfb),
            Mode::Amplifier => (&mut self.parameters, &mut self.amplifier),
        }
    }

    fn draw_mode_selector(&mut self, ui: &mut Ui) -> bool {
        let mut changed = false;

        ui.horizontal(|ui| {
            ui.heading("Mode: ");
            changed |= self.selected_mode.selectors(ui);
        });

        changed
    }

    fn draw_view_selector(mode: &mut dyn ModeUi, ui: &mut Ui) -> (bool, bool) {
        let mut changed = false;
        let mut reset_requested = false;

        ui.horizontal(|ui| {
            ui.heading("View: ");
            changed |= mode.draw_view_selector(ui);

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                reset_requested = Self::reset_button(ui);
                if let Some(compute_time) = mode.compute_time() {
                    let milliseconds = compute_time.as_secs_f64() * 1_000.0;
                    ui.label(format!("Compute: {milliseconds:.3} ms"));
                }
            });
        });

        (changed, reset_requested)
    }

    fn reset_button(ui: &mut Ui) -> bool {
        let shortcut = egui::KeyboardShortcut::new(egui::Modifiers::NONE, egui::Key::R);

        let shortcut_pressed = !ui.ctx().wants_keyboard_input()
            && ui.input_mut(|input| input.consume_shortcut(&shortcut));

        let shortcut_text = ui.ctx().format_shortcut(&shortcut);

        let button_clicked = ui
            .add(egui::Button::new("Reset").shortcut_text(shortcut_text))
            .clicked();

        button_clicked || shortcut_pressed
    }

    fn reset_active_mode(&mut self) {
        self.parameters = LaserParameters::default();
        match self.selected_mode {
            Mode::Dfb => self.dfb.reset(&self.parameters),
            Mode::Amplifier => self.amplifier.reset(&self.parameters),
        }
        self.dfb.clear_cached_plot();
        self.amplifier.clear_cached_plot();
    }
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let mut changed = false;

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::both().show(ui, |ui| {
                changed |= self.draw_mode_selector(ui);

                ui.separator();

                let (view_changed, reset_requested) = {
                    let mode = self.active_mode_mut();
                    Self::draw_view_selector(mode, ui)
                };
                changed |= view_changed;
                if reset_requested {
                    self.reset_active_mode();
                    changed = true;
                }

                ui.separator();

                let (parameters, mode) = self.active_mode_and_parameters_mut();
                changed |= mode.draw_controls(parameters, ui);
                if changed || !mode.has_cached_plot() {
                    mode.recompute_plot(parameters);
                    ctx.request_repaint();
                }

                ui.separator();

                mode.draw_plot(ui);
            });
        });
    }
}
