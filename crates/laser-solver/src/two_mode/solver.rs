use crate::error::SolverError;
use crate::fibre::{BidirectionalAmplitude, transfer};
use crate::grating::GratingModel;
use crate::lase::{DopantModel, FieldState};
use crate::maths::picard::PicardConfig;
use crate::two_mode::{FieldProfile, Pump, ResolvedFibre, Signal};

pub struct TwoModeSolver<'a, D: DopantModel, G: GratingModel> {
    fibre: &'a ResolvedFibre<'a, D, G>,
}

impl<D: DopantModel, G: GratingModel> TwoModeSolver<'_, D, G> {
    pub fn solve_injected(
        &self,
        pump: Pump,
        signal: Signal,
        picard_config: PicardConfig,
    ) -> Result<FieldProfile, SolverError> {
        let use_shooting = pump.balance == 1.0 && signal.balance == 1.0;
        let solution = if use_shooting {
            vec![FieldState::default(); self.fibre.grid.points()]
        } else {
            vec![FieldState::default(); self.fibre.grid.points()]
        };

        Result::Ok::<FieldProfile, SolverError>(FieldProfile::new(
            self.fibre.grid.positions().collect(),
            solution,
        ))
    }

    pub fn solve_lasing() -> Result<FieldProfile, SolverError> {
        Result::Ok::<FieldProfile, SolverError>(FieldProfile::new(
            vec![0.0],
            vec![FieldState::default()],
        ))
    }

    fn injected_left_boundary(
        &self,
        pump: Pump,
        signal: Signal,
        profile: &[FieldState],
    ) -> FieldState {
        let fibre = self.fibre;
        let grid = fibre.grid();
        assert_eq!(profile.len(), grid.points());

        let dz = grid.dz();
        let mut pump_optical_depth = 0.0;
        let mut signal_transfer_bottom_row = (0.0, 1.0); // calculate bottom row of transfer matrix

        for (&field, &kappa) in profile[..grid.steps()].iter().zip(fibre.kappas()).rev() {
            let gain = fibre.gain(field);
            pump_optical_depth += 0.5 * gain.pump * dz;
            signal_transfer_bottom_row = update_signal_transfer_bottom_row(
                signal_transfer_bottom_row,
                gain.signal,
                kappa,
                dz,
            );
        }

        let (signal_forward, signal_backward_right) = signal.amplitudes();
        let (pump_forward, pump_backward_right) = pump.amplitudes();
        let (t21, t22) = signal_transfer_bottom_row;

        FieldState {
            signal: BidirectionalAmplitude {
                forward: signal_forward,
                backward: (signal_backward_right - t21 * signal_forward) / t22,
            },
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: pump_backward_right * pump_optical_depth.exp(),
            },
        }
    }
}

fn update_signal_transfer_bottom_row(
    (t21, t22): (f64, f64),
    gain: f64,
    kappa: f64,
    dz: f64,
) -> (f64, f64) {
    if kappa == 0.0 {
        let forward_factor = (0.5 * gain * dz).exp();
        (t21 * forward_factor, t22 / forward_factor)
    } else {
        let (a, b, c, d) = transfer(gain, kappa, dz);
        (t21 * a + t22 * c, t21 * b + t22 * d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dopant::{TwoLevelCrossSections, TwoLevelDopant};
    use crate::fibre::{Fibre, FibreGeometry, FieldMode};
    use crate::grating::{NoGrating, PiShift};

    fn zero_gain_fibre<G: GratingModel>(grating: G) -> Fibre<TwoLevelDopant, G> {
        Fibre {
            geometry: FibreGeometry {
                core_radius: 4e-6,
                numerical_aperture: 0.1,
                length: 1.0,
            },
            dopant: TwoLevelDopant {
                density: 0.0,
                lifetime: 1.0,
            },
            grating,
        }
    }

    fn resolve_zero_gain<G: GratingModel>(
        fibre: &Fibre<TwoLevelDopant, G>,
        steps: usize,
    ) -> ResolvedFibre<'_, TwoLevelDopant, G> {
        fibre.resolve_with_interactions(
            FieldMode::new(970e-9),
            TwoLevelCrossSections::new(1.0, 0.0),
            FieldMode::new(1060e-9),
            TwoLevelCrossSections::new(0.0, 1.0),
            steps,
        )
    }

    fn assert_near(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
    }

    #[test]
    fn injected_boundary_reduces_to_inputs_without_gain_or_grating() {
        let fibre = zero_gain_fibre(NoGrating);
        let fibre = resolve_zero_gain(&fibre, 4);
        let solver = TwoModeSolver { fibre: &fibre };
        let signal = Signal {
            total: 5.0,
            balance: 0.2,
        };
        let pump = Pump {
            total: 7.0,
            balance: -0.4,
        };
        let profile = vec![FieldState::default(); fibre.grid().points()];

        let boundary = solver.injected_left_boundary(pump, signal, &profile);
        let (signal_forward, signal_backward) = signal.amplitudes();
        let (pump_forward, pump_backward) = pump.amplitudes();

        assert_eq!(boundary.signal.forward, signal_forward);
        assert_eq!(boundary.signal.backward, signal_backward);
        assert_eq!(boundary.pump.forward, pump_forward);
        assert_eq!(boundary.pump.backward, pump_backward);
    }

    #[test]
    fn injected_boundary_accounts_for_grating_transfer() {
        let fibre = zero_gain_fibre(PiShift {
            kappa_left: 0.7,
            kappa_right: 0.3,
            pi_shift_position: 0.5,
        });
        let fibre = resolve_zero_gain(&fibre, 4);
        let solver = TwoModeSolver { fibre: &fibre };
        let signal = Signal {
            total: 5.0,
            balance: 0.2,
        };
        let pump = Pump {
            total: 7.0,
            balance: -0.4,
        };
        let profile = vec![FieldState::default(); fibre.grid().points()];

        let boundary = solver.injected_left_boundary(pump, signal, &profile);
        let mut propagated_signal = boundary.signal;
        for &kappa in fibre.kappas() {
            propagated_signal = propagated_signal.coupled_step(0.0, kappa, fibre.grid().dz());
        }

        assert_near(propagated_signal.backward, signal.backward_amplitude());
        assert_near(boundary.pump.backward, pump.backward_amplitude());
    }
}
