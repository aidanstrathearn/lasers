//! Lightweight native plotting helpers for exploratory solver work.

use eframe::egui;
use egui_plot::{HLine, Legend, Line, LineStyle, Plot, PlotPoints, VLine};

const MATPLOTLIB_COLORS: [egui::Color32; 10] = [
    egui::Color32::from_rgb(31, 119, 180),
    egui::Color32::from_rgb(255, 127, 14),
    egui::Color32::from_rgb(44, 160, 44),
    egui::Color32::from_rgb(214, 39, 40),
    egui::Color32::from_rgb(148, 103, 189),
    egui::Color32::from_rgb(140, 86, 75),
    egui::Color32::from_rgb(227, 119, 194),
    egui::Color32::from_rgb(127, 127, 127),
    egui::Color32::from_rgb(188, 189, 34),
    egui::Color32::from_rgb(23, 190, 207),
];

pub struct PlotLine {
    points: Vec<[f64; 2]>,
    label: Option<String>,
}

impl PlotLine {
    pub fn label(&mut self, label: impl Into<String>) -> &mut Self {
        self.label = Some(label.into());
        self
    }
}

pub struct ReferenceLine {
    value: f64,
    label: Option<String>,
}

#[allow(dead_code)]
impl ReferenceLine {
    pub fn label(&mut self, label: impl Into<String>) -> &mut Self {
        self.label = Some(label.into());
        self
    }
}

#[derive(Default)]
pub struct Plotter {
    series: Vec<PlotLine>,
    horizontal_lines: Vec<ReferenceLine>,
    vertical_lines: Vec<ReferenceLine>,
    x_label: String,
    y_label: String,
    title: String,
}

impl Plotter {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn plot(&mut self, x: &[f64], y: &[f64]) -> &mut PlotLine {
        assert_eq!(
            x.len(),
            y.len(),
            "x and y must contain the same number of values"
        );

        let points = x.iter().zip(y.iter()).map(|(&x, &y)| [x, y]).collect();

        self.series.push(PlotLine {
            points,
            label: None,
        });

        self.series.last_mut().unwrap()
    }

    #[allow(dead_code)]
    pub fn axhline(&mut self, y: f64) -> &mut ReferenceLine {
        self.horizontal_lines.push(ReferenceLine {
            value: y,
            label: None,
        });

        self.horizontal_lines.last_mut().unwrap()
    }

    #[allow(dead_code)]
    pub fn axvline(&mut self, x: f64) -> &mut ReferenceLine {
        self.vertical_lines.push(ReferenceLine {
            value: x,
            label: None,
        });

        self.vertical_lines.last_mut().unwrap()
    }

    pub fn xlabel(&mut self, label: impl Into<String>) {
        self.x_label = label.into();
    }

    pub fn ylabel(&mut self, label: impl Into<String>) {
        self.y_label = label.into();
    }

    pub fn title(&mut self, title: impl Into<String>) {
        self.title = title.into();
    }

    /// Opens a native window and blocks until it is closed.
    pub fn show(self) -> eframe::Result {
        let window_title = if self.title.is_empty() {
            "Plot".to_owned()
        } else {
            self.title.clone()
        };

        let options = eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default()
                .with_inner_size([800.0, 600.0])
                .with_min_inner_size([400.0, 300.0])
                .with_resizable(true),
            ..Default::default()
        };

        eframe::run_native(
            &window_title,
            options,
            Box::new(move |creation_context| {
                creation_context
                    .egui_ctx
                    .set_visuals(egui::Visuals::light());

                Ok(Box::new(PlotApp::new(self)))
            }),
        )
    }
}

struct PlotApp {
    plotter: Plotter,
    reset_view: bool,
}

impl PlotApp {
    fn new(plotter: Plotter) -> Self {
        Self {
            plotter,
            reset_view: false,
        }
    }
}

impl eframe::App for PlotApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            if !self.plotter.title.is_empty() {
                ui.vertical_centered(|ui| {
                    ui.heading(&self.plotter.title);
                });

                ui.add_space(4.0);
            }

            ui.horizontal(|ui| {
                if ui.button("Reset view").clicked() {
                    self.reset_view = true;
                }

                ui.label("Drag to pan • Scroll to zoom • Double-click to reset");
            });

            ui.add_space(4.0);

            let mut plot = Plot::new("main_plot")
                .legend(Legend::default())
                .x_axis_label(&self.plotter.x_label)
                .y_axis_label(&self.plotter.y_label)
                .allow_drag(true)
                .allow_scroll(true)
                .allow_zoom(true)
                .allow_boxed_zoom(true)
                .allow_double_click_reset(true);

            if self.reset_view {
                plot = plot.reset();
                self.reset_view = false;
            }

            plot.show(ui, |plot_ui| {
                for (index, line) in self.plotter.series.iter().enumerate() {
                    let colour = MATPLOTLIB_COLORS[index % MATPLOTLIB_COLORS.len()];
                    let legend_name = line.label.as_deref().unwrap_or_default();
                    let line_name = format!("series_{index}");
                    let points = line.points.clone();

                    let plot_line = Line::new(line_name, PlotPoints::from(points))
                        .name(legend_name)
                        .color(colour);

                    plot_ui.line(plot_line);
                }

                let colour_offset = self.plotter.series.len();

                for (index, line) in self.plotter.horizontal_lines.iter().enumerate() {
                    let colour =
                        MATPLOTLIB_COLORS[(colour_offset + index) % MATPLOTLIB_COLORS.len()];
                    let legend_name = line.label.as_deref().unwrap_or_default();

                    let plot_line = HLine::new(format!("hline_{index}"), line.value)
                        .name(legend_name)
                        .color(colour)
                        .style(LineStyle::dashed_dense());

                    plot_ui.hline(plot_line);
                }

                let colour_offset = colour_offset + self.plotter.horizontal_lines.len();

                for (index, line) in self.plotter.vertical_lines.iter().enumerate() {
                    let colour =
                        MATPLOTLIB_COLORS[(colour_offset + index) % MATPLOTLIB_COLORS.len()];
                    let legend_name = line.label.as_deref().unwrap_or_default();

                    let plot_line = VLine::new(format!("vline_{index}"), line.value)
                        .name(legend_name)
                        .color(colour)
                        .style(LineStyle::dashed_dense());

                    plot_ui.vline(plot_line);
                }
            });
        });
    }
}
