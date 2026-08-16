use eframe::egui;
use eframe::egui::{Response, Ui};
use egui_plot::{Legend, Line, Plot, PlotPoints};
use laser_solver::lase::{FibreParams, GratingProfile, GridPoints, Pump, dfb_solve};
use laser_solver::rootfind::{Newton1dConfig, RootFindError};
use std::sync::mpsc::Receiver;
use std::sync::{Arc, mpsc};
use std::thread;

pub fn f64_slider(handle: &mut f64, name: &str, min: f64, max: f64, step: f64, ui: &mut Ui) -> Response {
    ui.horizontal(|ui| {
        ui.label(name);
        ui.add(egui::Slider::new(handle, min..=max).step_by(step))
    })
    .inner
}

#[derive(PartialEq, Default)]
pub enum View {
    Cos,
    #[default]
    Sin,
}

type Points = Vec<[f64; 2]>;

#[derive(Default)]
pub struct ProfilePlot {
    pump: Pump,
    fibre_params: FibreParams,
    grid_points: GridPoints,
    grating: GratingProfile,
    pending: Option<Receiver<[Points; 4]>>,
    result: Option<[Points; 4]>
}

impl ProfilePlot {
    #[cfg(not(target_arch = "wasm32"))]
    fn start_compute(&mut self, ctx: egui::Context) {
        let full_profile = true;
        let nc = Newton1dConfig {
            tolerance: 1e-8f64,
            max_iters: 100usize,
            initial: self.pump.forward,
            dx: 1e-6,
        };
        let pump = self.pump;
        let fibre_params = self.fibre_params;
        let grid_points = self.grid_points;
        let grating = self.grating;
        let compute_fn = move || {
            let result = dfb_solve(pump, fibre_params, grid_points, grating, full_profile, nc)?;
            Ok([result.plotpoints("sgnl_f"),
                result.plotpoints("sgnl_b"),
                result.plotpoints("pump_f"),
                result.plotpoints("pump_f")])
        };

        let (tx, rx) = mpsc::channel();

        self.pending = Some(rx);

        thread::spawn(move || {
            //thread::sleep(Duration::from_millis(100));
            let points: Result<[Points;4], RootFindError> = compute_fn();
            let _ = tx.send(points.unwrap());
            ctx.request_repaint();
        });
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn fetch_result(&mut self) {
        // &self.pending: &Option<..>, and rust treats LHS as &Some(ref rx) -- generally derefs through pattern matching
        // e.g. (a, b) = &my_tuple derefs the outer tuple and makes a,b: &...
        // using self.pending.as_ref(): Option<&..> would mean nothing implicit happens
        if let Some(rx) = &self.pending {
            if let Ok(points) = rx.try_recv() {
                self.result = Some(points);
                self.pending = None;
            }
        }
    }
}

#[derive(Default)]
pub struct LaserApp {
    view: View,
    pump: Pump,
    fibre_params: FibreParams,
    grid_points: GridPoints,
    grating: GratingProfile,
    pending: Option<Receiver<Points>>,
    result: Option<Points>,
}

impl LaserApp {
    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::default()
    }

    pub fn compute_and_plot(&mut self, ui: &mut Ui) -> Result<(), RootFindError> {
        let full_profile = true;
        let nc = Newton1dConfig {
            tolerance: 1e-8f64,
            max_iters: 100usize,
            initial: self.pump.forward,
            dx: 1e-6,
        };
        let result = dfb_solve(
            self.pump,
            self.fibre_params,
            self.grid_points,
            self.grating,
            full_profile,
            nc,
        )?;
        let sgnl_f = result.plotpoints("sgnl_f");
        let sgnl_b = result.plotpoints("sgnl_b");
        let pump_f = result.plotpoints("pump_f");
        let pump_b = result.plotpoints("pump_b");
        Plot::new("field-profile")
            .legend(Legend::default())
            .x_axis_label("z")
            .y_axis_label("fields")
            .show(ui, |plot_ui| {
                plot_ui.line(Line::new("field", sgnl_f).name("Forward Signal").width(3.0));
                plot_ui.line(
                    Line::new("field", sgnl_b)
                        .name("Backward Signal")
                        .width(3.0),
                );
                plot_ui.line(Line::new("field", pump_f).name("Forward Pump").width(3.0));
                plot_ui.line(Line::new("field", pump_b).name("Backward Pump").width(3.0));
            });
        Ok(())
    }

    pub fn param_sliders(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label("pi pos");
            ui.add(
                egui::Slider::new(&mut self.grating.pi_shift_position, 0.01..=1.0).step_by(0.01),
            );
            ui.separator();

            ui.label("density");
            ui.add(egui::Slider::new(&mut self.fibre_params.density, 0.01..=10.0).step_by(0.01));
            ui.separator();

            ui.label("length");
            ui.add(egui::Slider::new(&mut self.fibre_params.length, 0.5..=20.0).step_by(0.01));
            ui.separator();
        });
        ui.horizontal(|ui| {
            ui.label("grid points");
            ui.add(egui::Slider::new(&mut self.grid_points.0, 10..=1000).step_by(0.01));
            ui.separator();

            ui.label("pump");
            ui.add(egui::Slider::new(&mut self.pump.forward, 0.1..=200.0).step_by(0.01));
            ui.separator();

            f64_slider(&mut self.fibre_params.pump_em, "pump_em", 0.0, 10.0, 0.01, ui);
            f64_slider(&mut self.fibre_params.lifetime, "lifetime", 0.01, 10.0, 0.01, ui);
        });
    }

    pub fn view_selectors(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.view, View::Sin, "Sin");
            ui.selectable_value(&mut self.view, View::Cos, "Cos");
        });
    }


}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Lasers");
            self.param_sliders(ui);
            self.view_selectors(ui);
            self.compute_and_plot(ui).unwrap_or_else(|error| {
                ui.colored_label(ui.visuals().error_fg_color, error.to_string());
            });
        });
    }
}
