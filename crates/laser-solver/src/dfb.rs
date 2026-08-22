use crate::error::SolverError;
use crate::lase::{FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump, gain};
use crate::picard::{PicardConfig, dfb_pump_scan_picard, dfb_solve_picard};
use crate::rootfind::{RootFindConfig, rootfind_1d};
use crate::utils::{IterationConfig, relative_diff};

impl FieldState {
    pub fn propagate(self, fp: FibreParams, kappa: f64, dz: f64) -> Self {
        self.general_step(self, fp, kappa, dz)
    }

    pub fn general_step(self, other: Self, fp: FibreParams, kappa: f64, dz: f64) -> Self {
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

pub fn solve_profile(fs: FieldState, fp: FibreParams, dz: f64, kappas: &[f64]) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(kappas.len() + 1);
    result.push(current);
    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
        result.push(current);
    }
    result
}

pub fn out_field(fs: FieldState, fp: FibreParams, dz: f64, kappas: &[f64]) -> FieldState {
    let mut current = fs;
    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
    }
    current
}

pub fn dfb_solve_shooting(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
) -> Result<FieldProfile, SolverError> {
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let trial = |sgnl_b| FieldState {
        sgnl_f: 0.0,
        sgnl_b: sgnl_b,
        pump_f: pu.forward,
        pump_b: 0.0, // shooting method requires pump.backward = 0
    };
    let f = |sgnl_b| out_field(trial(sgnl_b), fp, dz, &kappas).sgnl_b;
    let sgnl_b = rootfind_1d(f, config)?;

    if full_profile {
        let z = gp.grid(fp.length);
        let fields = solve_profile(trial(sgnl_b), fp, dz, &kappas);
        Ok(FieldProfile::new(z, fields))
    } else {
        let z = vec![0.0_f64, fp.length];
        let out_left = trial(sgnl_b);
        let fields = vec![out_left, out_field(out_left, fp, dz, &kappas)];
        Ok(FieldProfile::new(z, fields))
    }
}

pub fn dfb_output_power_shooting(
    pump_power: f64,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> (f64, f64, bool) {
    let pu = Pump {
        forward: pump_power.sqrt(),
        backward: 0.0,
    };

    dfb_solve_shooting(pu, fp, gp, kp, false, config).map_or((0.0, 0.0, false), |result| {
        (
            result.sgnl_f().last().unwrap().powi(2),
            result.sgnl_b().next().unwrap().powi(2),
            true,
        )
    })
}

pub fn dfb_pump_scan_shooting(
    pumps: &[f64],
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Vec<(f64, f64, bool)> {
    pumps
        .iter()
        .map(|&pmp| dfb_output_power_shooting(pmp, fp, gp, kp, config))
        .collect()
}

pub fn find_threshold_and_slope(
    pump_start: f64,
    pump_step: f64,
    ip: IterationConfig,
    mut output_power: impl FnMut(f64) -> (f64, f64, bool),
) -> Result<(f64, f64, f64), SolverError> {
    let mut current_pump = pump_start;
    let mut total_diff = -1.0;
    let mut sf = 0.0;
    let mut sb = 0.0;
    for _ in 0..ip.max {
        let (new_sf, new_sb, success) = output_power(current_pump);
        if !success {
            current_pump += pump_step;
            continue;
        }

        let new_total_diff = (new_sf + new_sb) - (sb + sf);

        if relative_diff(new_total_diff, total_diff) < ip.tol && new_total_diff > 0.0 {
            let slope_f = (new_sf - sf) / pump_step;
            let slope_b = (new_sb - sb) / pump_step;
            let threshold = current_pump - (new_sf + new_sb) / (slope_b + slope_f);
            return Ok((slope_f, slope_b, threshold));
        } else {
            current_pump += pump_step;
            total_diff = new_total_diff;
            sb = new_sb;
            sf = new_sf;
        }
    }
    Err(SolverError::ThresholdNotFound)
}

pub fn dfb_find_threshold_and_slope_shooting(
    pump_start: f64,
    pump_step: f64,
    ip: IterationConfig,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Result<(f64, f64, f64), SolverError> {

    let f  = |p| dfb_output_power_shooting(p, fp, gp, kp, config);

    find_threshold_and_slope(pump_start, pump_step, ip, f)
}

pub fn dfb_pump_scan(
    pumps: &[f64],
    balance: f64,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
    picard_config: PicardConfig,
) -> Vec<(f64, f64, bool)> {
    let use_picard = balance != 1.0;
    if use_picard {
        dfb_pump_scan_picard(pumps, balance, fp, gp, kp, config, picard_config)
    } else {
        dfb_pump_scan_shooting(pumps, fp, gp, kp, config)
    }
}

pub fn dfb_solve(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    picard_config: PicardConfig,
) -> Result<FieldProfile, SolverError> {
    let use_picard = pu.backward.abs() > 0.0;
    if use_picard {
        dfb_solve_picard(pu, fp, gp, kp, full_profile, config, picard_config)
    } else {
        dfb_solve_shooting(pu, fp, gp, kp, full_profile, config)
    }
}
