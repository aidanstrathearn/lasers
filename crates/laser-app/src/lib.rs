use eframe::egui;
use eframe::egui::Ui;
use egui_plot::{Legend, Line, Plot};
use laser_solver::lase::{
    FibreParams, GratingProfile, GridPoints, Pump, dfb_solve, dfb_threshold_curve_with_zeros,
    linspace,
};
use laser_solver::rootfind::{BisectionConfig, Midpoint, Newton1dConfig, RootFindError};
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread;

// pub fn f64_slider(
//     handle: &mut f64,
//     name: &str,
//     min: f64,
//     max: f64,
//     step: f64,
//     ui: &mut Ui,
// ) -> Response {
//     ui.horizontal(|ui| {
//         ui.label(name);
//         ui.add(egui::Slider::new(handle, min..=max).step_by(step))
//     })
//     .inner
// }

#[derive(PartialEq, Default)]
pub enum View {
    Profile,
    #[default]
    Threshold,
}

type Points = Vec<[f64; 2]>;

#[allow(dead_code)]
#[derive(Default)]
pub struct ProfilePlot {
    pump: Pump,
    fibre_params: FibreParams,
    grid_points: GridPoints,
    grating: GratingProfile,
    pending: Option<Receiver<[Points; 4]>>,
    result: Option<[Points; 4]>,
}

impl ProfilePlot {
    #[allow(dead_code)]
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
            Ok([
                result.plotpoints("sgnl_f"),
                result.plotpoints("sgnl_b"),
                result.plotpoints("pump_f"),
                result.plotpoints("pump_f"),
            ])
        };

        let (tx, rx) = mpsc::channel();

        self.pending = Some(rx);

        thread::spawn(move || {
            //thread::sleep(Duration::from_millis(100));
            let points: Result<[Points; 4], RootFindError> = compute_fn();
            let _ = tx.send(points.unwrap());
            ctx.request_repaint();
        });
    }
    #[allow(dead_code)]
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
}

impl LaserApp {
    pub fn new(creation_context: &eframe::CreationContext<'_>) -> Self {
        creation_context
            .egui_ctx
            .set_visuals(egui::Visuals::light());
        Self::default()
    }

    pub fn threshold_plot(&mut self, ui: &mut Ui) {
        let bc = BisectionConfig {
            tolerance: 1e-8f64,
            max_iters: 100usize,
            upper: self.pump.forward,
            lower: 1e-6,
            midpoint: Midpoint::Geometric,
        };

        let pumps = linspace(0.1, 20.0, 20);
        let threshold = dfb_threshold_curve_with_zeros(
            &pumps,
            self.fibre_params,
            self.grid_points,
            self.grating,
            bc,
        );
        let sgnl_f = threshold.iter().map(|x| x.0);
        let sgnl_b = threshold.iter().map(|x| x.1);
        let sgnl_f_points: Points = pumps
            .iter()
            .zip(sgnl_f)
            .map(|(&x, y)| [x, y.abs()])
            .collect();
        let sgnl_b_points: Points = pumps
            .iter()
            .zip(sgnl_b)
            .map(|(&x, y)| [x, y.abs()])
            .collect();
        Plot::new("threshold")
            .legend(Legend::default())
            .x_axis_label("pump")
            .y_axis_label("signal")
            .show(ui, |plot_ui| {
                plot_ui.line(
                    Line::new("field", sgnl_f_points)
                        .name("Forward Signal")
                        .width(3.0),
                );
                plot_ui.line(
                    Line::new("field", sgnl_b_points)
                        .name("Backward Signal")
                        .width(3.0),
                );
            });
    }

    pub fn profile_plot(&mut self, ui: &mut Ui) -> Result<(), RootFindError> {
        let full_profile = true;
        let nc = Newton1dConfig {
            tolerance: 1e-8f64,
            max_iters: 200usize,
            initial: self.pump.forward,
            dx: 1e-6,
        };
        let bc = BisectionConfig {
            tolerance: 1e-9f64,
            max_iters: 200usize,
            upper: self.pump.forward,
            lower: 1e-10,
            midpoint: Midpoint::Geometric,
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

    pub fn view_selectors(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.view, View::Threshold, "Threshold");
            ui.selectable_value(&mut self.view, View::Profile, "Profile");
        });
    }
}

fn grating_slider_grid(grating: &mut GratingProfile, ui: &mut Ui) {
    egui::Grid::new("grating").show(ui, |ui| {
        ui.label("kappa");
        ui.add(egui::Slider::new(&mut grating.kappa_max, 0.0..=10.0).step_by(0.01));
        ui.end_row();

        ui.label("pi pos");
        ui.add(egui::Slider::new(&mut grating.pi_shift_position, 0.0..=1.0).step_by(0.01));
        ui.end_row();
    });
}

fn fibre_params_slider_grid(params: &mut FibreParams, ui: &mut Ui) {
    egui::Grid::new("params").show(ui, |ui| {
        egui::Grid::new("params1").show(ui, |ui| {
            ui.label("pump_em");
            ui.add(egui::Slider::new(&mut params.pump_em, 0.0..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("pump_ab");
            ui.add(egui::Slider::new(&mut params.pump_ab, 0.0..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("sgnl_em");
            ui.add(egui::Slider::new(&mut params.sgnl_em, 0.0..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("sgnl_ab");
            ui.add(egui::Slider::new(&mut params.sgnl_ab, 0.0..=10.0).step_by(0.01));
            ui.end_row();
        });

        egui::Grid::new("params2").show(ui, |ui| {
            ui.label("density");
            ui.add(egui::Slider::new(&mut params.density, 0.1..=10.0).step_by(0.01));
            ui.end_row();

            ui.label("lifetime");
            ui.add(egui::Slider::new(&mut params.lifetime, 0.1..=2.0).step_by(0.01));
            ui.end_row();

            ui.label("length");
            ui.add(egui::Slider::new(&mut params.length, 0.1..=50.0).step_by(0.01));
            ui.end_row();
        });
    });
}

impl eframe::App for LaserApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("Lasers");
                ui.separator();
                self.view_selectors(ui);
            });

            egui::Grid::new("all-params").show(ui, |ui| {
                fibre_params_slider_grid(&mut self.fibre_params, ui);
                grating_slider_grid(&mut self.grating, ui);
                ui.end_row();
            });

            match self.view {
                View::Threshold => self.threshold_plot(ui),
                View::Profile => self.profile_plot(ui).unwrap_or_else(|error| {
                    ui.colored_label(ui.visuals().error_fg_color, error.to_string());
                }),
            };
        });
    }
}
