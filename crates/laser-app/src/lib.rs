use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

pub struct LaserApp {
    frequency: f64,
}

impl LaserApp {
    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::default()
    }
}

impl Default for LaserApp {
    fn default() -> Self {
        Self { frequency: 0.05 }
    }
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("eframe + egui_plot");
            ui.horizontal(|ui| {
                ui.add(
                    egui::Slider::new(&mut self.frequency, 0.01..=1.0)
                        .text("frequency")
                        .step_by(0.01),
                );
            });

            let points = PlotPoints::from_iter((0..=200).map(|i| {
                let x = f64::from(i);
                [x, (x * self.frequency).sin()]
            }));

            Plot::new("sine-wave")
                .x_axis_label("x")
                .y_axis_label("sin(x)")
                .show(ui, |plot_ui| {
                    plot_ui.line(Line::new("sine", points));
                });
        });
    }
}
