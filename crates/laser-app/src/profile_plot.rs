use crate::{LaserApp, Points};
use eframe::egui;
use eframe::egui::Ui;
use egui_plot::{Legend, Line, Plot};
use laser_solver::dfb::{dfb_solve, dfb_solve_shooting};
use laser_solver::error::SolverError;
use laser_solver::lase::{FibreParams, FieldProfile, GratingProfile, GridPoints, Pump};
use laser_solver::picard::PicardConfig;
use laser_solver::rootfind::{BisectionConfig, Newton1dConfig};
use laser_solver::utils::IterationConfig;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread;

trait FieldProfileExt {
    fn plotpoints(&self, field: &str) -> Vec<[f64; 2]>;
}

impl FieldProfileExt for FieldProfile {
    fn plotpoints(&self, field: &str) -> Vec<[f64; 2]> {
        let z = self.z();
        match field {
            "sgnl_b" => {
                let field = self.sgnl_b();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }

            "sgnl_f" => {
                let field = self.sgnl_f();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }

            "pump_f" => {
                let field = self.pump_f();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }

            "pump_b" => {
                let field = self.pump_b();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }
            _ => panic!(),
        }
    }
}

impl LaserApp {
    pub fn profile_plot(&mut self, ui: &mut Ui) -> Result<(), SolverError> {
        let full_profile = true;
        // let nc = Newton1dConfig {
        //     tolerance: 1e-8f64,
        //     max_iters: 200usize,
        //     initial: self.pump.forward,
        //     dx: 1e-6,
        // };
        let bc = BisectionConfig {
            upper: (self.pump.forward.powi(2) + self.pump.backward.powi(2)).sqrt(),
            ..self.config
        };

        let picard_config = PicardConfig {
            max_iterations: 5_000,
            relative_tolerance: 1e-6,
            absolute_tolerance: 1e-10,
        };
        let result = dfb_solve(
            self.pump,
            self.fibre_params,
            self.grid_points,
            self.grating,
            full_profile,
            bc,
            picard_config,
        )?;
        let sgnl_f = result.plotpoints("sgnl_f");
        let sgnl_b = result.plotpoints("sgnl_b");
        let pump_f = result.plotpoints("pump_f");
        let pump_b = result.plotpoints("pump_b");
        Plot::new("field-profile")
            .legend(Legend::default())
            .x_axis_label(egui::RichText::new("z").size(24.0))
            .y_axis_label(egui::RichText::new("fields").size(24.0))
            .show(ui, |plot_ui| {
                plot_ui.line(
                    Line::new("profile sf", sgnl_f)
                        .name("Forward Signal")
                        .width(3.0),
                );
                plot_ui.line(
                    Line::new("profile sb", sgnl_b)
                        .name("Backward Signal")
                        .width(3.0),
                );
                plot_ui.line(
                    Line::new("profile pf", pump_f)
                        .name("Forward Pump")
                        .width(3.0),
                );
                plot_ui.line(
                    Line::new("profile pb", pump_b)
                        .name("Backward Pump")
                        .width(3.0),
                );
            });
        Ok(())
    }
}

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
            iteration: IterationConfig::default(),
            initial: self.pump.forward,
            dx: 1e-6,
        };
        let pump = self.pump;
        let fibre_params = self.fibre_params;
        let grid_points = self.grid_points;
        let grating = self.grating;
        let compute_fn = move || {
            let result =
                dfb_solve_shooting(pump, fibre_params, grid_points, grating, full_profile, nc)?;
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
            let points: Result<[Points; 4], SolverError> = compute_fn();
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
