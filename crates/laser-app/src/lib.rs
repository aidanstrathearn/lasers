use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

#[derive(PartialEq)]
pub enum View {
    Cos,
    Sin,
}
pub struct LaserApp {
    frequency: f64,
    view: View,
}

impl LaserApp {
    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::default()
    }

    pub fn points(&self) -> Vec<[f64; 2]> {
        let x = (0..=200).map(|i| f64::from(i));
        match self.view {
            View::Sin => x.map(move |x| [x, (x * self.frequency).sin()]).collect(),
            View::Cos => x.map(move |x| [x, (x * self.frequency).cos()]).collect(),
        }
    }
}

impl Default for LaserApp {
    fn default() -> Self {
        Self {
            frequency: 0.05,
            view: View::Sin,
        }
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

            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.view, View::Sin, "Sin");
                ui.selectable_value(&mut self.view, View::Cos, "Cos");
            });

            let points = PlotPoints::from(self.points());

            Plot::new("sine-wave")
                .x_axis_label("x")
                .y_axis_label("sin(x)")
                .show(ui, |plot_ui| {
                    plot_ui.line(Line::new("sine", points));
                });
        });
    }
}
