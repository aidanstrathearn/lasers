use crate::plotter::Plotter;
use crate::{LaserApp, Points, timed};
use eframe::egui;
use laser_solver::dfb::{dfb_solve, dfb_solve_shooting, Grating};
use laser_solver::error::SolverError;
use laser_solver::lase::{Fibre, FieldProfile, GridPoints, Pump};
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
                z.zip(field).map(|(x, y)| [x, y.powi(2)]).collect()
            }

            "sgnl_f" => {
                let field = self.sgnl_f();
                z.zip(field).map(|(x, y)| [x, y.powi(2)]).collect()
            }

            "pump_f" => {
                let field = self.pump_f();
                z.zip(field).map(|(x, y)| [x, y.powi(2)]).collect()
            }

            "pump_b" => {
                let field = self.pump_b();
                z.zip(field).map(|(x, y)| [x, y.powi(2)]).collect()
            }
            _ => panic!(),
        }
    }
}

impl LaserApp {
    pub fn profile_plot(&mut self) -> Result<Plotter, SolverError> {
        let full_profile = true;
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };

        let (result, compute_time) = timed(|| {
            dfb_solve(
                self.pump,
                self.fibre_params,
                self.grid_points,
                self.grating,
                full_profile,
                bc,
                self.picard_config,
            )
        });
        let result = result?;
        let sgnl_f = result.plotpoints("sgnl_f");
        let sgnl_b = result.plotpoints("sgnl_b");
        let pump_f = result.plotpoints("pump_f");
        let pump_b = result.plotpoints("pump_b");

        let mut plt = Plotter::new();
        plt.add_points(sgnl_f).label("Forward signal");
        plt.add_points(sgnl_b).label("Backward signal");
        plt.add_points(pump_f).label("Forward pump");
        plt.add_points(pump_b).label("Backward pump");
        plt.xlabel("z");
        plt.ylabel("Power");
        plt.set_compute_time(compute_time);
        Ok(plt)
    }
}

#[allow(dead_code)]
#[derive(Default)]
pub struct ProfilePlot {
    pump: Pump,
    fibre_params: Fibre,
    grid_points: GridPoints,
    grating: Grating,
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
            initial: self.pump.forward_amplitude(),
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
        if let Some(rx) = &self.pending &&
            let Ok(points) = rx.try_recv() {
                self.result = Some(points);
                self.pending = None;
            }
        
    }
}
