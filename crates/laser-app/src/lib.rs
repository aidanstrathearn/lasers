use eframe::egui;
use eframe::egui::Ui;
use egui_plot::{Line, Plot, PlotPoints};
use laser_solver::lase::{FibreParams, GratingProfile, GridPoints, Pump, dfb_solve};
use laser_solver::rootfind::Newton1dConfig;

#[derive(PartialEq)]
pub enum View {
    Cos,
    Sin,
}
pub struct LaserApp {
    frequency: f64,
    view: View,
    fibre_params: FibreParams,
    grating: GratingProfile,
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

    pub fn plot(&self, ui: &mut Ui) {
        let points = PlotPoints::from(self.points());

        Plot::new("sine-wave")
            .x_axis_label("x")
            .y_axis_label("sin(x)")
            .show(ui, |plot_ui| {
                plot_ui.line(Line::new("sine", points));
            });
    }

    pub fn plot_fields(&self, ui: &mut Ui) {
        let pu = Pump {
            forward: 100.0,
            backward: 0.0,
        };
        let gp = GridPoints(500);
        let full_profile = true;
        let nc = Newton1dConfig {
            tolerance: 1e-8f64,
            max_iters: 100usize,
            initial: pu.forward,
            dx: 1e-6,
        };
        let result = dfb_solve(pu, self.fibre_params, gp, self.grating, full_profile, nc).unwrap();
        let sgnl_f = result.plotpoints("sgnl_f");
        let sgnl_b = result.plotpoints("sgnl_b");
        let pump_f = result.plotpoints("pump_f");
        let pump_b = result.plotpoints("pump_b");
        Plot::new("field-profile")
            .x_axis_label("z")
            .y_axis_label("fields")
            .show(ui, |plot_ui| {
                plot_ui.line(Line::new("field", sgnl_f));
                plot_ui.line(Line::new("field", sgnl_b));
                plot_ui.line(Line::new("field", pump_f));
                plot_ui.line(Line::new("field", pump_b));
            });
    }
}

impl Default for LaserApp {
    fn default() -> Self {
        let fp = FibreParams {
            density: 1.0,
            lifetime: 1.0,
            pump_ab: 0.01,
            pump_em: 0.0,
            sgnl_ab: 0.0,
            sgnl_em: 1.0,
            length: 10.0,
        };

        let kp = GratingProfile {
            kappa_max: 1.0,
            pi_shift_position: 0.45,
        };

        Self {
            frequency: 0.05,
            view: View::Sin,
            fibre_params: fp,
            grating: kp,
        }
    }
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("eframe + egui_plot");
            ui.horizontal(|ui| {
                ui.label("pi pos");
                ui.add(
                    egui::Slider::new(&mut self.grating.pi_shift_position, 0.01..=1.0)
                        .step_by(0.01),
                );
                ui.separator();

                ui.label("density");
                ui.add(
                    egui::Slider::new(&mut self.fibre_params.density, 0.01..=10.0).step_by(0.01),
                );
                ui.separator();

                ui.label("length");
                ui.add(egui::Slider::new(&mut self.fibre_params.length, 0.5..=20.0).step_by(0.01));
                ui.separator();
            });
            ui.vertical(|ui| {
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut self.view, View::Sin, "Sin");
                    ui.selectable_value(&mut self.view, View::Cos, "Cos");
                });
                self.plot_fields(ui);
            })
        });
    }
}
