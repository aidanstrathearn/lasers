pub mod picard;
mod shooting;

use self::picard::{PicardConfig, PicardSolver, dfb_output_power_picard, dfb_solve_picard};
use crate::error::SolverError;
use crate::lase::{
    Fibre, FieldProfile, FieldState, GridPoints, Pump, PumpScan, find_threshold_and_slope, gain,
    pump_scan,
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
            dfb_solve_picard(
                pump,
                self.fibre,
                solve_config.grid_points,
                self.grating,
                full_profile,
                solve_config.root_find,
                solve_config.picard,
            )
        } else {
            self.solve_shooting(pump, solve_config, full_profile)
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

pub fn dfb_find_threshold_and_slope(
    pump_start: Pump,
    pump_step: f64,
    ip: IterationConfig,
    fp: Fibre,
    gp: GridPoints,
    kp: Grating,
    config: impl Into<RootFindConfig> + Copy,
    picard_config: PicardConfig,
) -> Result<(f64, f64, f64), SolverError> {
    pump_start.amplitudes();
    let use_picard = pump_start.balance != 1.0;
    if use_picard {
        let mut solver = PicardSolver::new(pump_start, fp, gp);
        let f = |total| {
            dfb_output_power_picard(
                Pump {
                    total,
                    ..pump_start
                },
                fp,
                gp,
                kp,
                config,
                &mut solver,
                picard_config,
            )
        };
        find_threshold_and_slope(pump_start.total, pump_step, ip, f)
    } else {
        let laser = DfbLaser {
            fibre: fp,
            grating: kp,
        };
        let solve_config = DfbSolveConfig {
            grid_points: gp,
            root_find: config.into(),
            picard: picard_config,
        };
        let f = |total| {
            laser.output_power_shooting(
                Pump {
                    total,
                    ..pump_start
                },
                solve_config,
            )
        };
        find_threshold_and_slope(pump_start.total, pump_step, ip, f)
    }
}

pub fn dfb_pump_scan(
    pumps: &[f64],
    balance: f64,
    fp: Fibre,
    gp: GridPoints,
    kp: Grating,
    config: impl Into<RootFindConfig> + Copy,
    picard_config: PicardConfig,
) -> Result<PumpScan, SolverError> {
    let Some(&pump_start) = pumps.first() else {
        return Ok(Vec::new());
    };

    let use_picard = balance != 1.0;
    if use_picard {
        let mut solver = PicardSolver::new(
            Pump {
                total: pump_start,
                balance,
            },
            fp,
            gp,
        );
        pump_scan(pumps, |total| {
            dfb_output_power_picard(
                Pump { total, balance },
                fp,
                gp,
                kp,
                config,
                &mut solver,
                picard_config,
            )
        })
    } else {
        let laser = DfbLaser {
            fibre: fp,
            grating: kp,
        };
        let solve_config = DfbSolveConfig {
            grid_points: gp,
            root_find: config.into(),
            picard: picard_config,
        };
        pump_scan(pumps, |total| {
            laser.output_power_shooting(Pump { total, balance }, solve_config)
        })
    }
}

pub fn dfb_solve(
    pump: Pump,
    fp: Fibre,
    gp: GridPoints,
    kp: Grating,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    picard_config: PicardConfig,
) -> Result<FieldProfile, SolverError> {
    let use_picard = pump.backward_amplitude() > 0.0;
    if use_picard {
        dfb_solve_picard(pump, fp, gp, kp, full_profile, config, picard_config)
    } else {
        DfbLaser {
            fibre: fp,
            grating: kp,
        }
        .solve_shooting(
            pump,
            DfbSolveConfig {
                grid_points: gp,
                root_find: config.into(),
                picard: picard_config,
            },
            full_profile,
        )
    }
}
