use crate::fibre::BidirectionalAmplitude;
use crate::two_mode::{FieldState, Gain};

impl BidirectionalAmplitude {
    pub fn coupled_step(self, gain: f64, kappa: f64, dz: f64) -> Self {
        let (a, b, c, d) = transfer(gain, kappa, dz);
        Self {
            forward: a * self.forward + b * self.backward,
            backward: c * self.forward + d * self.backward,
        }
    }

    pub fn uncoupled_step(self, gain: f64, dz: f64) -> Self {
        let factor = (0.5 * gain * dz).exp();
        Self {
            forward: self.forward * factor,
            backward: self.backward / factor,
        }
    }
}

impl FieldState {
    pub fn coupled_step(self, gain: Gain, kappa: f64, dz: f64) -> Self {
        Self {
            signal: self.signal.coupled_step(gain.signal, kappa, dz),
            pump: self.pump.uncoupled_step(gain.pump, dz),
        }
    }

    pub fn uncoupled_step(self, gain: Gain, dz: f64) -> Self {
        Self {
            signal: self.signal.uncoupled_step(gain.signal, dz),
            pump: self.pump.uncoupled_step(gain.pump, dz),
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
            signal: BidirectionalAmplitude {
                forward: 2.0,
                backward: 3.0,
            },
            pump: BidirectionalAmplitude {
                forward: 5.0,
                backward: 7.0,
            },
        };
        let gain = Gain {
            pump: 4.0,
            signal: 2.0,
        };

        let stepped = fields.uncoupled_step(gain, 0.5);

        assert_eq!(stepped.signal.forward, 2.0 * 0.5_f64.exp());
        assert_eq!(stepped.signal.backward, 3.0 / 0.5_f64.exp());
        assert_eq!(stepped.pump.forward, 5.0 * 1.0_f64.exp());
        assert_eq!(stepped.pump.backward, 7.0 / 1.0_f64.exp());
    }

    #[test]
    fn profile_gain_is_evaluated_from_each_current_field() {
        let evaluated_at = RefCell::new(Vec::new());
        let initial = FieldState {
            pump: BidirectionalAmplitude {
                forward: 2.0,
                ..BidirectionalAmplitude::default()
            },
            ..FieldState::default()
        };

        let profile = solve_profile_uncoupled(
            initial,
            |fields| {
                evaluated_at.borrow_mut().push(fields.pump.forward);
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
            vec![profile[0].pump.forward, profile[1].pump.forward]
        );
    }
}
