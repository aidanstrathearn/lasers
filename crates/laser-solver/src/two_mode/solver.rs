use crate::error::SolverError;
use crate::fibre::{BidirectionalAmplitude, transfer};
use crate::grating::GratingModel;
use crate::lase::{DopantModel, FieldState, profile_convergence_error};
use crate::maths::picard::{PicardConfig, PicardSolver};
use crate::maths::rootfind::RootFindConfig;
use crate::two_mode::propagation::solve_profile;
use crate::two_mode::{FieldProfile, Pump, ResolvedFibre, Signal};

pub struct TwoModeSolver<'a, D: DopantModel, G: GratingModel> {
    fibre: &'a ResolvedFibre<'a, D, G>,
}

impl<'a, D: DopantModel, G: GratingModel> TwoModeSolver<'a, D, G> {
    pub fn new(fibre: &'a ResolvedFibre<'a, D, G>) -> Self {
        Self { fibre }
    }

    pub fn solve_injected(
        &self,
        pump: Pump,
        signal: Signal,
        root_find_config: RootFindConfig,
        picard_config: PicardConfig,
    ) -> Result<FieldProfile, SolverError> {
        // Reserved for injected solves that require shooting/root-finding.
        let _ = root_find_config;
        let dz = self.fibre.grid.dz();
        let kappas = self.fibre.kappas();
        let (signal_forward, signal_backward) = signal.amplitudes();
        let (pump_forward, pump_backward) = pump.amplitudes();
        let left_boundary = FieldState {
            signal: BidirectionalAmplitude {
                forward: signal_forward,
                backward: signal_backward,
            },
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: pump_backward,
            },
        };
        let use_shooting = pump_backward == 0.0
            && signal_backward == 0.0
            && kappas.iter().all(|&kappa| kappa == 0.0); // reflections affect boundary

        let solution = if use_shooting {
            solve_profile(left_boundary, |fields| self.fibre.gain(fields), dz, kappas)
        } else {
            let mut solver = PicardSolver::filled(self.fibre.grid.points(), left_boundary);
            let set_boundary =
                |current: &[FieldState]| self.injected_left_boundary(pump, signal, current);
            let step = |new_previous: &FieldState, old_current: &FieldState, i| {
                new_previous.step_if(self.fibre.gain(*old_current), kappas[i], dz)
            };

            let error = |current: &[FieldState], previous: &[FieldState]| {
                profile_convergence_error(
                    current,
                    previous,
                    picard_config.absolute_tolerance,
                    picard_config.relative_tolerance,
                )
            };

            solver.solve(picard_config.max_iterations, set_boundary, step, error)?;

            solver.profile().to_vec()
        };

        Ok(FieldProfile::new(
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

        // todo: need to catch when t22=0
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
    use crate::maths::rootfind::BisectionConfig;

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

    #[test]
    fn forward_injected_grating_solve_satisfies_both_signal_boundaries() {
        let fibre = zero_gain_fibre(PiShift {
            kappa_left: 0.7,
            kappa_right: 0.3,
            pi_shift_position: 0.5,
        });
        let fibre = resolve_zero_gain(&fibre, 4);
        let solver = TwoModeSolver::new(&fibre);
        let signal = Signal {
            total: 5.0,
            balance: 1.0,
        };
        let pump = Pump {
            total: 7.0,
            balance: 1.0,
        };

        let profile = solver
            .solve_injected(
                pump,
                signal,
                BisectionConfig::default().into(),
                PicardConfig::default(),
            )
            .expect("forward-injected grating solve should converge");
        let left = profile.fields.first().unwrap();
        let right = profile.fields.last().unwrap();

        assert_near(left.signal.forward, signal.forward_amplitude());
        assert_near(right.signal.backward, signal.backward_amplitude());
        assert_near(left.pump.forward, pump.forward_amplitude());
        assert_near(right.pump.backward, pump.backward_amplitude());
        assert!(
            left.signal.backward.abs() > 1e-6,
            "grating solve should produce a reflected signal"
        );
    }

    #[test]
    fn forward_injected_no_grating_solve_preserves_zero_gain_fields() {
        let fibre = zero_gain_fibre(NoGrating);
        let fibre = resolve_zero_gain(&fibre, 4);
        let solver = TwoModeSolver::new(&fibre);
        let signal = Signal {
            total: 5.0,
            balance: 1.0,
        };
        let pump = Pump {
            total: 7.0,
            balance: 1.0,
        };

        let profile = solver
            .solve_injected(
                pump,
                signal,
                BisectionConfig::default().into(),
                PicardConfig::default(),
            )
            .expect("forward-injected no-grating solve should succeed");
        let expected = FieldState {
            signal: BidirectionalAmplitude {
                forward: signal.forward_amplitude(),
                backward: 0.0,
            },
            pump: BidirectionalAmplitude {
                forward: pump.forward_amplitude(),
                backward: 0.0,
            },
        };

        assert_eq!(profile.fields.len(), fibre.grid().points());
        for field in profile.fields {
            assert_near(field.signal.forward, expected.signal.forward);
            assert_near(field.signal.backward, expected.signal.backward);
            assert_near(field.pump.forward, expected.pump.forward);
            assert_near(field.pump.backward, expected.pump.backward);
        }
    }

    #[test]
    fn bidirectional_no_grating_solve_satisfies_all_injected_boundaries() {
        let fibre = zero_gain_fibre(NoGrating);
        let fibre = resolve_zero_gain(&fibre, 4);
        let solver = TwoModeSolver::new(&fibre);
        let signal = Signal {
            total: 2.0,
            balance: 0.0,
        };
        let pump = Pump {
            total: 2.0,
            balance: 0.0,
        };

        let profile = solver
            .solve_injected(
                pump,
                signal,
                BisectionConfig::default().into(),
                PicardConfig::default(),
            )
            .expect("bidirectional no-grating solve should converge");
        let left = profile.fields.first().unwrap();
        let right = profile.fields.last().unwrap();

        assert_near(left.signal.forward, signal.forward_amplitude());
        assert_near(right.signal.backward, signal.backward_amplitude());
        assert_near(left.pump.forward, pump.forward_amplitude());
        assert_near(right.pump.backward, pump.backward_amplitude());
    }
}
