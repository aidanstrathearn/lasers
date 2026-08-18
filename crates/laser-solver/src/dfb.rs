use crate::lase::{FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump, gain};
use crate::rootfind::{RootFindConfig, RootFindError, rootfind_1d};

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

pub fn dfb_solve(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
) -> Result<FieldProfile, RootFindError> {
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let trial = |sgnl_b| FieldState {
        sgnl_f: 0.0,
        sgnl_b: sgnl_b,
        pump_f: pu.forward,
        pump_b: 0.0, // bidirectional pump not implemented yet
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

pub fn dfb_threshold_curve(
    pumps: Vec<f64>,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Vec<Result<f64, RootFindError>> {
    pumps
        .iter()
        .map(|&pmp| {
            let pu = Pump {
                forward: pmp,
                backward: 0.0,
            };
            let result = dfb_solve(pu, fp, gp, kp, false, config)?;
            Ok(result.sgnl_f().last().unwrap())
        })
        .collect()
}

pub fn dfb_threshold_curve_with_zeros(
    pumps: &[f64],
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Vec<(f64, f64)> {
    let full_profile = false;
    pumps
        .iter()
        .map(|&pmp| {
            let pu = Pump {
                forward: pmp,
                backward: 0.0,
            };
            let result = dfb_solve(pu, fp, gp, kp, full_profile, config)?;
            Ok((
                result.sgnl_f().last().unwrap(),
                result.sgnl_b().next().unwrap(),
            ))
        })
        .map(|result: Result<(f64, f64), RootFindError>| result.unwrap_or((0.0, 0.0)))
        .collect()

    // pumps
    //     .iter()
    //     .map(|&pmp| {
    //         let pu = Pump {
    //             forward: pmp,
    //             backward: 0.0,
    //         };
    //
    //         dfb_solve(pu, fp, gp, kp, full_profile, config).map_or((0.0, 0.0), |result| {
    //             (
    //                 result.sgnl_f().last().unwrap(),
    //                 result.sgnl_b().next().unwrap(),
    //             )
    //         })
    //     })
    //     .collect()
}
