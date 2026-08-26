use crate::lase::{Fibre, FieldState, gain};

impl FieldState {
    pub fn coupled_step_shooting(self, fp: Fibre, kappa: f64, dz: f64) -> Self {
        self.coupled_step_general(self, fp, kappa, dz)
    }

    pub fn coupled_step_general(self, other: Self, fp: Fibre, kappa: f64, dz: f64) -> Self {
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

    pub fn uncoupled_step_shooting(self, fp: Fibre, dz: f64) -> Self {
        self.uncoupled_step_general(self, fp, dz)
    }

    pub fn uncoupled_step_general(self, other: Self, fibre: Fibre, dz: f64) -> Self {
        let (pump_gain, signal_gain) = gain(other, fibre);

        let pump_factor = (0.5 * pump_gain * dz).exp();
        let signal_factor = (0.5 * signal_gain * dz).exp();

        Self {
            sgnl_f: self.sgnl_f * signal_factor,
            sgnl_b: self.sgnl_b / signal_factor,
            pump_f: self.pump_f * pump_factor,
            pump_b: self.pump_b / pump_factor,
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

pub fn solve_profile_uncoupled(
    fs: FieldState,
    fp: Fibre,
    dz: f64,
    nsteps: usize,
) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(nsteps + 1);
    result.push(current);
    for _ in 0..nsteps {
        current = current.uncoupled_step_shooting(fp, dz);
        result.push(current);
    }
    result
}

pub fn out_field_uncoupled(fs: FieldState, fp: Fibre, dz: f64, nsteps: usize) -> FieldState {
    let mut current = fs;
    for _ in 0..nsteps {
        current = current.uncoupled_step_shooting(fp, dz);
    }
    current
}

pub fn solve_profile_coupled(fs: FieldState, fp: Fibre, dz: f64, kappas: &[f64]) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(kappas.len() + 1);
    result.push(current);
    for &kappa in kappas {
        current = current.coupled_step_shooting(fp, kappa, dz);
        result.push(current);
    }
    result
}

pub fn out_field_coupled(fs: FieldState, fp: Fibre, dz: f64, kappas: &[f64]) -> FieldState {
    let mut current = fs;
    for &kappa in kappas {
        current = current.coupled_step_shooting(fp, kappa, dz);
    }
    current
}