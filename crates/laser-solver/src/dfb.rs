use crate::error::SolverError;
use crate::lase::{
    FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, OutputPower, Pump, PumpScan,
    find_threshold_and_slope, gain, pump_scan,
};
use crate::picard::{PicardConfig, PicardDfbSolver, dfb_output_power_picard, dfb_solve_picard};
use crate::rootfind::{RootFindConfig, rootfind_1d};
use crate::utils::IterationConfig;

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
    pump: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
) -> Result<FieldProfile, SolverError> {
    let (pump_forward, pump_backward) = pump.amplitudes();
    assert_eq!(
        pump_backward, 0.0,
        "shooting solver requires a forward-only pump"
    );
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let trial = |sgnl_b| FieldState {
        sgnl_f: 0.0,
        sgnl_b,
        pump_f: pump_forward,
        pump_b: 0.0, // shooting method requires zero backward pump amplitude
    };
    let f = |sgnl_b| out_field(trial(sgnl_b), fp, dz, &kappas).sgnl_b / sgnl_b;
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
    pump: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Result<OutputPower, SolverError> {
    let profile = dfb_solve_shooting(pump, fp, gp, kp, false, config)?;
    Ok(profile.output_powers())
}

pub fn dfb_find_threshold_and_slope(
    pump_start: Pump,
    pump_step: f64,
    ip: IterationConfig,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
    picard_config: PicardConfig,
) -> Result<(f64, f64, f64), SolverError> {
    pump_start.amplitudes();
    let use_picard = pump_start.balance != 1.0;
    if use_picard {
        let mut solver = PicardDfbSolver::new(pump_start, fp, gp);
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
        let f = |total| {
            dfb_output_power_shooting(
                Pump {
                    total,
                    ..pump_start
                },
                fp,
                gp,
                kp,
                config,
            )
        };
        find_threshold_and_slope(pump_start.total, pump_step, ip, f)
    }
}

pub fn dfb_pump_scan(
    pumps: &[f64],
    balance: f64,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
    picard_config: PicardConfig,
) -> Result<PumpScan, SolverError> {
    let Some(&pump_start) = pumps.first() else {
        return Ok(Vec::new());
    };

    let use_picard = balance != 1.0;
    if use_picard {
        let mut solver = PicardDfbSolver::new(
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
        pump_scan(pumps, |total| {
            dfb_output_power_shooting(Pump { total, balance }, fp, gp, kp, config)
        })
    }
}

pub fn dfb_solve(
    pump: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    picard_config: PicardConfig,
) -> Result<FieldProfile, SolverError> {
    let use_picard = pump.backward_amplitude() > 0.0;
    if use_picard {
        dfb_solve_picard(pump, fp, gp, kp, full_profile, config, picard_config)
    } else {
        dfb_solve_shooting(pump, fp, gp, kp, full_profile, config)
    }
}
