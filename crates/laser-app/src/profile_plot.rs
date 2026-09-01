use crate::plotter::Plotter;
use crate::{LaserParameters, dfb::DfbMode, field_profile_plot, timed};
use laser_solver::error::SolverError;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::two_mode::TwoModeSolver;

impl DfbMode {
    pub fn profile_plot(&mut self, parameters: &LaserParameters) -> Result<Plotter, SolverError> {
        let full_profile = true;

        let (result, compute_time) = timed(|| {
            let fibre = parameters.resolved_fibre();
            let bc = BisectionConfig {
                upper: 2.0 * fibre.pump_flux(parameters.pump.total).sqrt(),
                ..self.config
            };
            TwoModeSolver::new(&fibre, self.steps).solve_lasing(
                parameters.pump,
                bc.into(),
                self.picard_config,
                full_profile,
            )
        });
        self.compute_time = Some(compute_time);
        Ok(field_profile_plot(&result?))
    }
}

// #[allow(dead_code)]
// pub struct ProfilePlot {
//     pump: Pump,
//     fibre_params: Fibre<TwoLevelDopant, PiShift>,
//     pump_mode: FieldMode,
//     sgnl_mode: FieldMode,
//     pump_interaction: TwoLevelCrossSections,
//     signal_interaction: TwoLevelCrossSections,
//     steps: usize,
//     pending: Option<Receiver<[Points; 4]>>,
//     result: Option<[Points; 4]>,
// }
//
// impl Default for ProfilePlot {
//     fn default() -> Self {
//         Self {
//             pump: Pump::default(),
//             fibre_params: Fibre::default(),
//             pump_mode: FieldMode::default(),
//             sgnl_mode: FieldMode::default(),
//             pump_interaction: TwoLevelCrossSections::new(0.01, 0.0),
//             signal_interaction: TwoLevelCrossSections::new(0.0, 1.0),
//             steps: 100,
//             pending: None,
//             result: None,
//         }
//     }
// }
//
// impl ProfilePlot {
//     #[allow(dead_code)]
//     #[cfg(not(target_arch = "wasm32"))]
//     fn start_compute(&mut self, ctx: egui::Context) {
//         let full_profile = true;
//         let nc = Newton1dConfig {
//             iteration: IterationConfig::default(),
//             initial: self.pump.forward_amplitude(),
//             dx: 1e-6,
//         };
//         let pump = self.pump;
//         let fibre_params = self.fibre_params.clone();
//         let pump_mode = self.pump_mode;
//         let sgnl_mode = self.sgnl_mode;
//         let pump_interaction = self.pump_interaction;
//         let signal_interaction = self.signal_interaction;
//         let steps = self.steps;
//         let compute_fn = move || {
//             let fibre = fibre_params.resolve_with_interactions(
//                 pump_mode,
//                 pump_interaction,
//                 sgnl_mode,
//                 signal_interaction,
//             );
//             let laser = TwoModeSolver::new(&fibre, steps);
//             let result = laser.solve_lasing(
//                 pump,
//                 nc.into(),
//                 Default::default(),
//                 full_profile,
//             )?;
//             Ok([
//                 power_points(result.z(), result.sgnl_f()),
//                 power_points(result.z(), result.sgnl_b()),
//                 power_points(result.z(), result.pump_f()),
//                 power_points(result.z(), result.pump_f()),
//             ])
//         };
//
//         let (tx, rx) = mpsc::channel();
//
//         self.pending = Some(rx);
//
//         thread::spawn(move || {
//             //thread::sleep(Duration::from_millis(100));
//             let points: Result<[Points; 4], SolverError> = compute_fn();
//             let _ = tx.send(points.unwrap());
//             ctx.request_repaint();
//         });
//     }
//     #[allow(dead_code)]
//     #[cfg(not(target_arch = "wasm32"))]
//     fn fetch_result(&mut self) {
//         // &self.pending: &Option<..>, and rust treats LHS as &Some(ref rx) -- generally derefs through pattern matching
//         // e.g. (a, b) = &my_tuple derefs the outer tuple and makes a,b: &...
//         // using self.pending.as_ref(): Option<&..> would mean nothing implicit happens
//         if let Some(rx) = &self.pending
//             && let Ok(points) = rx.try_recv()
//         {
//             self.result = Some(points);
//             self.pending = None;
//         }
//     }
// }
