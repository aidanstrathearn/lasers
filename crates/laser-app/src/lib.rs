mod controls;
mod dfb;
mod kappa_plot;
mod pi_pos_plot;
mod plotter;
mod pop_plot;
mod profile_plot;
mod residual_plot;
mod threshold_plot;

use crate::dfb::DfbMode;
use crate::plotter::Plotter;
use eframe::egui;
use eframe::egui::Ui;
use laser_solver::lase::FieldProfile;
use std::time::Duration;
use web_time::Instant;

type Points = Vec<[f64; 2]>;

fn timed<T>(compute: impl FnOnce() -> T) -> (T, Duration) {
    let start = Instant::now();
    let result = compute();
    (result, start.elapsed())
}

fn power_points(z: impl Iterator<Item = f64>, amplitude: impl Iterator<Item = f64>) -> Points {
    z.zip(amplitude)
        .map(|(z, amplitude)| [z, amplitude.powi(2)])
        .collect()
}

fn field_profile_plot(profile: &FieldProfile, compute_time: Duration) -> Plotter {
    let mut plot = Plotter::new();
    plot.add_points(power_points(profile.z(), profile.sgnl_f()))
        .label("Forward signal");
    plot.add_points(power_points(profile.z(), profile.sgnl_b()))
        .label("Backward signal");
    plot.add_points(power_points(profile.z(), profile.pump_f()))
        .label("Forward pump");
    plot.add_points(power_points(profile.z(), profile.pump_b()))
        .label("Backward pump");
    plot.xlabel("z");
    plot.ylabel("Power");
    plot.set_compute_time(compute_time);
    plot
}

trait ModeUi {
    fn draw_view_selector(&mut self, ui: &mut Ui) -> bool;
    fn draw_controls(&mut self, ui: &mut Ui) -> bool;
    fn reset(&mut self);
    fn has_cached_plot(&self) -> bool;
    fn recompute_plot(&mut self);
    fn compute_time(&self) -> Option<Duration>;
    fn draw_plot(&self, ui: &mut Ui);
}

#[derive(Default, PartialEq, Copy, Clone)]
enum Mode {
    #[default]
    Dfb,
}

impl Mode {
    fn selectors(&mut self, ui: &mut Ui) -> bool {
        ui.selectable_value(self, Self::Dfb, "DFB").changed()
    }
}

#[derive(Default)]
pub struct LaserApp {
    selected_mode: Mode,
    dfb: DfbMode,
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

    fn draw_view_selector(mode: &mut dyn ModeUi, ui: &mut Ui) -> bool {
        let mut changed = false;

        ui.horizontal(|ui| {
            ui.heading("View: ");
            changed |= mode.draw_view_selector(ui);

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                Self::reset_button(mode, ui);
                if let Some(compute_time) = mode.compute_time() {
                    let milliseconds = compute_time.as_secs_f64() * 1_000.0;
                    ui.label(format!("Compute: {milliseconds:.3} ms"));
                }
            });
        });

        changed
    }

    fn reset_button(mode: &mut dyn ModeUi, ui: &mut Ui) {
        let shortcut = egui::KeyboardShortcut::new(egui::Modifiers::NONE, egui::Key::R);

        let shortcut_pressed = !ui.ctx().wants_keyboard_input()
            && ui.input_mut(|input| input.consume_shortcut(&shortcut));

        let shortcut_text = ui.ctx().format_shortcut(&shortcut);

        let button_clicked = ui
            .add(egui::Button::new("Reset").shortcut_text(shortcut_text))
            .clicked();

        if button_clicked || shortcut_pressed {
            mode.reset();
        }
    }
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let mut changed = false;

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::both().show(ui, |ui| {
                changed |= self.draw_mode_selector(ui);

                ui.separator();

                let mode = self.active_mode_mut();
                changed |= Self::draw_view_selector(mode, ui);

                ui.separator();

                changed |= mode.draw_controls(ui);
                if changed || !mode.has_cached_plot() {
                    mode.recompute_plot();
                    ctx.request_repaint();
                }

                ui.separator();

                mode.draw_plot(ui);
            });
        });
    }
}
