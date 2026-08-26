pub mod picard;
mod shooting;

use self::picard::{PicardConfig, PicardSolver};
use crate::error::SolverError;
use crate::lase::{
    Fibre, FieldProfile, FieldState, GridPoints, Pump, PumpScan,
    find_threshold_and_slope as scan_for_threshold, gain, pump_scan as scan_pump_totals,
};
use crate::rootfind::RootFindConfig;
use crate::utils::IterationConfig;

#[derive(Copy, Clone)]
pub struct DfbSolveConfig {
    pub grid_points: GridPoints,
    pub root_find: RootFindConfig,
    pub picard: PicardConfig,
}

pub struct DfbLaser {
    pub fibre: Fibre,
    pub grating: Grating,
}

impl DfbLaser {
    pub fn solve(
        &self,
        pump: Pump,
        solve_config: DfbSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let use_picard = pump.backward_amplitude() > 0.0;
        if use_picard {
            self.solve_picard(pump, solve_config, full_profile)
        } else {
            self.solve_shooting(pump, solve_config, full_profile)
        }
    }

    pub fn find_threshold_and_slope(
        &self,
        pump_start: Pump,
        pump_step: f64,
        iteration: IterationConfig,
        solve_config: DfbSolveConfig,
    ) -> Result<(f64, f64, f64), SolverError> {
        pump_start.amplitudes();
        let use_picard = pump_start.balance != 1.0;
        if use_picard {
            let mut solver = PicardSolver::new(pump_start, self.fibre, solve_config.grid_points);
            let f = |total| {
                self.output_power_picard(
                    Pump {
                        total,
                        ..pump_start
                    },
                    solve_config,
                    &mut solver,
                )
            };
            scan_for_threshold(pump_start.total, pump_step, iteration, f)
        } else {
            let f = |total| {
                self.output_power_shooting(
                    Pump {
                        total,
                        ..pump_start
                    },
                    solve_config,
                )
            };
            scan_for_threshold(pump_start.total, pump_step, iteration, f)
        }
    }

    pub fn pump_scan(
        &self,
        pump_totals: &[f64],
        balance: f64,
        solve_config: DfbSolveConfig,
    ) -> Result<PumpScan, SolverError> {
        let Some(&pump_start) = pump_totals.first() else {
            return Ok(Vec::new());
        };

        let use_picard = balance != 1.0;
        if use_picard {
            let mut solver = PicardSolver::new(
                Pump {
                    total: pump_start,
                    balance,
                },
                self.fibre,
                solve_config.grid_points,
            );
            scan_pump_totals(pump_totals, |total| {
                self.output_power_picard(Pump { total, balance }, solve_config, &mut solver)
            })
        } else {
            scan_pump_totals(pump_totals, |total| {
                self.output_power_shooting(Pump { total, balance }, solve_config)
            })
        }
    }
}

#[derive(Copy, Clone)]
pub struct Grating {
    pub kappa_left: f64,
    pub kappa_right: f64,
    pub pi_shift_position: f64,
}

impl Default for Grating {
    fn default() -> Self {
        Self {
            kappa_left: 1.0,
            kappa_right: 1.0,
            pi_shift_position: 0.45,
        }
    }
}

impl Grating {
    pub fn grid(self, n: usize) -> Vec<f64> {
        (0..n)
            .map(|j| {
                let z = j as f64 / n as f64;
                if z < self.pi_shift_position {
                    self.kappa_left
                } else {
                    -self.kappa_right
                }
            })
            .collect()
    }
}

impl FieldState {
    pub fn propagate(self, fp: Fibre, kappa: f64, dz: f64) -> Self {
        self.general_step(self, fp, kappa, dz)
    }

    pub fn general_step(self, other: Self, fp: Fibre, kappa: f64, dz: f64) -> Self {
        let (gp, gs) = gain(other, fp);
        let (a, b, c, d) = transfer(gs, kappa, dz);
        let expg = (0.5 * gp * dz).exp();

        FieldState {
            sgnl_f: a * self.sgnl_f + b * self.sgnl_b,
            sgnl_b: c * self.sgnl_f + d * self.sgnl_b,
            pump_f: self.pump_f * expg,
            pump_b: self.pump_b / expg,
        }
    }
}

pub fn transfer(gain: f64, kappa: f64, dz: f64) -> (f64, f64, f64, f64) {
    let g_dz = 0.5 * gain * dz;
    let k_dz = kappa * dz;
    let x = (g_dz * g_dz + k_dz * k_dz).sqrt();

    let cosh = x.cosh();
    let sinch = if x > 1e-30 { x.sinh() / x } else { 1.0_f64 };
    (
        cosh + g_dz * sinch,
        k_dz * sinch,
        k_dz * sinch,
        cosh - g_dz * sinch,
    )
}

pub fn solve_profile(fs: FieldState, fp: Fibre, dz: f64, kappas: &[f64]) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(kappas.len() + 1);
    result.push(current);
    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
        result.push(current);
    }
    result
}

pub fn out_field(fs: FieldState, fp: Fibre, dz: f64, kappas: &[f64]) -> FieldState {
    let mut current = fs;
    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
    }
    current
}
