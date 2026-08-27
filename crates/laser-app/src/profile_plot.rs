use crate::plotter::Plotter;
use crate::{Points, dfb::DfbMode, field_profile_plot, power_points, timed};
use eframe::egui;
use laser_solver::dfb::{DfbLaser, DfbSolveConfig, Grating};
use laser_solver::error::SolverError;
use laser_solver::lase::{Fibre, GridPoints, Pump};
use laser_solver::rootfind::{BisectionConfig, Newton1dConfig};
use laser_solver::utils::IterationConfig;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread;

impl DfbMode {
    pub fn profile_plot(&mut self) -> Result<Plotter, SolverError> {
        let full_profile = true;
        let bc = BisectionConfig {
            upper: 2.0 * self.pump.total.sqrt(),
            ..self.config
        };

        let (result, compute_time) = timed(|| {
            self.dfb_laser()
                .solve(self.pump, self.dfb_solve_config(bc), full_profile)
        });
        self.compute_time = Some(compute_time);
        Ok(field_profile_plot(&result?))
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
            let laser = DfbLaser {
                fibre: fibre_params,
                grating,
            };
            let result = laser.solve_shooting(
                pump,
                DfbSolveConfig {
                    grid_points,
                    root_find: nc.into(),
                    picard: Default::default(),
                },
                full_profile,
            )?;
            Ok([
                power_points(result.z(), result.sgnl_f()),
                power_points(result.z(), result.sgnl_b()),
                power_points(result.z(), result.pump_f()),
                power_points(result.z(), result.pump_f()),
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
        if let Some(rx) = &self.pending
            && let Ok(points) = rx.try_recv()
        {
            self.result = Some(points);
            self.pending = None;
        }
    }
}
