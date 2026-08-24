use crate::Points;
use eframe::egui;
use eframe::egui::Ui;
use egui_plot::{HLine, Legend, Line, LineStyle, Plot, VLine};
use std::hash::Hash;
use std::time::Duration;

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
    points: Points,
    label: Option<String>,
}

pub struct ReferenceLine {
    value: f64,
    label: Option<String>,
}

impl PlotLine {
    pub fn label(&mut self, label: impl Into<String>) -> &mut Self {
        self.label = Some(label.into());
        self
    }
}

impl ReferenceLine {
    pub fn label(&mut self, label: impl Into<String>) -> &mut Self {
        self.label = Some(label.into());
        self
    }
}

#[allow(dead_code)]
#[derive(Default)]
pub struct Plotter {
    series: Vec<PlotLine>,
    horizontal_lines: Vec<ReferenceLine>,
    vertical_lines: Vec<ReferenceLine>,
    x_label: String,
    y_label: String,
    title: String,
    x_limits: Option<(f64, f64)>,
    compute_time: Duration,
}

#[allow(dead_code)]
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

    pub fn add_points(&mut self, points: Points) -> &mut PlotLine {
        self.series.push(PlotLine {
            points,
            label: None,
        });

        self.series.last_mut().unwrap()
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

    pub fn xlim(&mut self, lower: f64, upper: f64) {
        self.x_limits = Some((lower, upper));
    }

    pub fn set_compute_time(&mut self, compute_time: Duration) {
        self.compute_time = compute_time;
    }

    pub fn compute_time(&self) -> Duration {
        self.compute_time
    }

    pub fn axhline(&mut self, y: f64) -> &mut ReferenceLine {
        self.horizontal_lines.push(ReferenceLine {
            value: y,
            label: None,
        });

        self.horizontal_lines.last_mut().unwrap()
    }

    pub fn axvline(&mut self, x: f64) -> &mut ReferenceLine {
        self.vertical_lines.push(ReferenceLine {
            value: x,
            label: None,
        });

        self.vertical_lines.last_mut().unwrap()
    }

    pub fn show(&self, ui: &mut Ui, id: impl Hash) {
        let plot_id = egui::Id::new(id);
        ui.style_mut().text_styles.insert(
            egui::TextStyle::Body,
            egui::FontId::proportional(24.0), // this is for ticks + legend text
        );
        Plot::new(plot_id)
            .legend(Legend::default())
            .x_axis_label(egui::RichText::new(&self.x_label).size(24.0))
            .y_axis_label(egui::RichText::new(&self.y_label).size(24.0))
            .show(ui, |plot_ui| {
                if let Some((lower, upper)) = self.x_limits {
                    plot_ui.set_plot_bounds_x(lower..=upper);
                }
                let colour_offset = self.series.len();

                for (index, line) in self.series.iter().enumerate() {
                    let colour = MATPLOTLIB_COLORS[index % MATPLOTLIB_COLORS.len()];
                    let legend_name = line.label.as_deref().unwrap_or_default();
                    let plot_line = Line::new(legend_name, line.points.clone())
                        .color(colour)
                        .id(plot_id.with(index))
                        .width(3.0);

                    plot_ui.line(plot_line);
                }

                for (index, line) in self.horizontal_lines.iter().enumerate() {
                    let colour =
                        MATPLOTLIB_COLORS[(colour_offset + index) % MATPLOTLIB_COLORS.len()];
                    let legend_name = line.label.as_deref().unwrap_or_default();

                    let plot_line = HLine::new(format!("hline_{index}"), line.value)
                        .name(legend_name)
                        .color(colour)
                        .width(3.0)
                        .style(LineStyle::dashed_dense());

                    plot_ui.hline(plot_line);
                }

                let colour_offset = colour_offset + self.horizontal_lines.len();

                for (index, line) in self.vertical_lines.iter().enumerate() {
                    let colour =
                        MATPLOTLIB_COLORS[(colour_offset + index) % MATPLOTLIB_COLORS.len()];
                    let legend_name = line.label.as_deref().unwrap_or_default();

                    let plot_line = VLine::new(format!("vline_{index}"), line.value)
                        .name(legend_name)
                        .color(colour)
                        .width(3.0)
                        .style(LineStyle::dashed_dense());

                    plot_ui.vline(plot_line);
                }
            });
    }
}
