use crate::lase::{FieldState, Gain};

impl FieldState {
    pub fn coupled_step(self, gain: Gain, kappa: f64, dz: f64) -> Self {
        let (a, b, c, d) = transfer(gain.signal, kappa, dz);
        let expg = (0.5 * gain.pump * dz).exp();

        FieldState {
            sgnl_f: a * self.sgnl_f + b * self.sgnl_b,
            sgnl_b: c * self.sgnl_f + d * self.sgnl_b,
            pump_f: self.pump_f * expg,
            pump_b: self.pump_b / expg,
        }
    }

    pub fn uncoupled_step(self, gain: Gain, dz: f64) -> Self {
        let pump_factor = (0.5 * gain.pump * dz).exp();
        let signal_factor = (0.5 * gain.signal * dz).exp();

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
    gain: impl Fn(FieldState) -> Gain,
    dz: f64,
    nsteps: usize,
) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(nsteps + 1);
    result.push(current);
    for _ in 0..nsteps {
        current = current.uncoupled_step(gain(current), dz);
        result.push(current);
    }
    result
}

pub fn out_field_uncoupled(
    fs: FieldState,
    gain: impl Fn(FieldState) -> Gain,
    dz: f64,
    nsteps: usize,
) -> FieldState {
    let mut current = fs;
    for _ in 0..nsteps {
        current = current.uncoupled_step(gain(current), dz);
    }
    current
}

pub fn solve_profile_coupled(
    fs: FieldState,
    gain: impl Fn(FieldState) -> Gain,
    dz: f64,
    kappas: &[f64],
) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(kappas.len() + 1);
    result.push(current);
    for &kappa in kappas {
        current = current.coupled_step(gain(current), kappa, dz);
        result.push(current);
    }
    result
}

pub fn out_field_coupled(
    fs: FieldState,
    gain: impl Fn(FieldState) -> Gain,
    dz: f64,
    kappas: &[f64],
) -> FieldState {
    let mut current = fs;
    for &kappa in kappas {
        current = current.coupled_step(gain(current), kappa, dz);
    }
    current
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn uncoupled_step_applies_supplied_gain() {
        let fields = FieldState {
            sgnl_f: 2.0,
            sgnl_b: 3.0,
            pump_f: 5.0,
            pump_b: 7.0,
        };
        let gain = Gain {
            pump: 4.0,
            signal: 2.0,
        };

        let stepped = fields.uncoupled_step(gain, 0.5);

        assert_eq!(stepped.sgnl_f, 2.0 * 0.5_f64.exp());
        assert_eq!(stepped.sgnl_b, 3.0 / 0.5_f64.exp());
        assert_eq!(stepped.pump_f, 5.0 * 1.0_f64.exp());
        assert_eq!(stepped.pump_b, 7.0 / 1.0_f64.exp());
    }

    #[test]
    fn profile_gain_is_evaluated_from_each_current_field() {
        let evaluated_at = RefCell::new(Vec::new());
        let initial = FieldState {
            pump_f: 2.0,
            ..FieldState::default()
        };

        let profile = solve_profile_uncoupled(
            initial,
            |fields| {
                evaluated_at.borrow_mut().push(fields.pump_f);
                Gain {
                    pump: 2.0,
                    signal: 0.0,
                }
            },
            1.0,
            2,
        );

        assert_eq!(
            evaluated_at.into_inner(),
            vec![profile[0].pump_f, profile[1].pump_f]
        );
    }
}
