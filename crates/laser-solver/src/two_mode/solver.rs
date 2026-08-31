use crate::error::SolverError;
use crate::fibre::{BidirectionalAmplitude, transfer};
use crate::grating::{GratingModel, sample_grating};
use crate::lase::{
    DopantModel, FieldState, PumpScan, UniformGrid, profile_convergence_error,
    pump_scan as scan_pump_totals,
};
use crate::maths::picard::{PicardConfig, PicardError, PicardSolver};
use crate::maths::rootfind::{RootFindConfig, rootfind_1d, try_rootfind_1d};
use crate::two_mode::propagation::{out_field_coupled, solve_profile, solve_profile_coupled};
use crate::two_mode::{FieldProfile, OutputPower, Pump, ResolvedFibre, Signal};

pub struct TwoModeSolver<'a, D: DopantModel, G: GratingModel> {
    fibre: &'a ResolvedFibre<'a, D, G>,
    grid: UniformGrid,
    kappas: Box<[f64]>,
}

impl<'a, D: DopantModel, G: GratingModel> TwoModeSolver<'a, D, G> {
    pub fn new(fibre: &'a ResolvedFibre<'a, D, G>, steps: usize) -> Self {
        let grid = UniformGrid::new(fibre.length(), steps);
        let kappas = sample_grating(fibre.grating(), grid.steps()).into_boxed_slice();
        Self {
            fibre,
            grid,
            kappas,
        }
    }

    pub fn grid(&self) -> UniformGrid {
        self.grid
    }

    pub fn kappas(&self) -> &[f64] {
        &self.kappas
    }

    fn field_profile(&self, z: Vec<f64>, fields: Vec<FieldState>) -> FieldProfile {
        let (pump_flux_per_watt, signal_flux_per_watt) = self.fibre.flux_per_watt();
        FieldProfile::new(
            z,
            fields,
            pump_flux_per_watt,
            signal_flux_per_watt,
        )
    }

    fn solve_picard(
        &self,
        solver: &mut PicardSolver<FieldState>,
        config: PicardConfig,
        set_boundary: impl FnMut(&[FieldState]) -> FieldState,
    ) -> Result<(), PicardError> {
        assert_eq!(self.kappas.len() + 1, solver.profile().len());
        let dz = self.grid.dz();
        let step = |new_previous: &FieldState, old_current: &FieldState, i| {
            new_previous.step_if(self.fibre.gain(*old_current), self.kappas[i], dz)
        };
        let error = |current: &[FieldState], previous: &[FieldState]| {
            profile_convergence_error(
                current,
                previous,
                config.absolute_tolerance,
                config.relative_tolerance,
            )
        };

        solver
            .solve(config.max_iterations, set_boundary, step, error)
            .map(|_| ())
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
        let dz = self.grid.dz();
        let kappas = self.kappas();
        let (signal_forward, signal_backward) = self.fibre.signal_flux_amplitudes(signal);
        let (pump_forward, pump_backward) = self.fibre.pump_flux_amplitudes(pump);
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
            let mut solver = PicardSolver::filled(self.grid.points(), left_boundary);
            let set_boundary =
                |current: &[FieldState]| self.injected_left_boundary(pump, signal, current);
            self.solve_picard(&mut solver, picard_config, set_boundary)?;

            solver.profile().to_vec()
        };

        Ok(self.field_profile(
            self.grid.positions().collect(),
            solution,
        ))
    }

    pub fn solve_lasing(
        &self,
        pump: Pump,
        root_find_config: RootFindConfig,
        picard_config: PicardConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let profile = if self.fibre.pump_flux_amplitudes(pump).1 > 0.0 {
            self.solve_lasing_picard(pump, root_find_config, picard_config)?
        } else {
            self.solve_lasing_shooting(pump, root_find_config)?
        };

        Ok(if full_profile {
            profile
        } else {
            profile.into_endpoints()
        })
    }

    pub fn pump_scan(
        &self,
        pump_totals: &[f64],
        balance: f64,
        root_find_config: RootFindConfig,
        picard_config: PicardConfig,
    ) -> Result<PumpScan, SolverError> {
        scan_pump_totals(pump_totals, |total| {
            self.lasing_output(
                Pump { total, balance },
                root_find_config,
                picard_config,
            )
        })
    }

    fn lasing_output(
        &self,
        pump: Pump,
        root_find_config: RootFindConfig,
        picard_config: PicardConfig,
    ) -> Result<OutputPower, SolverError> {
        let profile = self.solve_lasing(pump, root_find_config, picard_config, false)?;
        Ok(profile.output_powers())
    }

    fn solve_lasing_shooting(
        &self,
        pump: Pump,
        root_find_config: RootFindConfig,
    ) -> Result<FieldProfile, SolverError> {
        let (pump_forward, pump_backward) = self.fibre.pump_flux_amplitudes(pump);
        assert_eq!(
            pump_backward, 0.0,
            "shooting solver requires a forward-only pump"
        );
        let grid = self.grid();
        let kappas = self.kappas();
        let dz = grid.dz();
        let trial = |sgnl_b| FieldState {
            signal: BidirectionalAmplitude {
                forward: 0.0,
                backward: sgnl_b,
            },
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: 0.0,
            },
        };
        let residual = |sgnl_b| {
            out_field_coupled(
                trial(sgnl_b),
                |fields| self.fibre.gain(fields),
                dz,
                kappas,
            )
            .signal
            .backward
                / sgnl_b
        };
        let sgnl_b = rootfind_1d(residual, root_find_config)?;

        Ok(self.field_profile(
            grid.positions().collect(),
            solve_profile_coupled(
                trial(sgnl_b),
                |fields| self.fibre.gain(fields),
                dz,
                kappas,
            ),
        ))
    }

    fn solve_lasing_picard(
        &self,
        pump: Pump,
        root_find_config: RootFindConfig,
        picard_config: PicardConfig,
    ) -> Result<FieldProfile, SolverError> {
        let (pump_forward, pump_backward) = self.fibre.pump_flux_amplitudes(pump);
        let initial = FieldState {
            signal: BidirectionalAmplitude::default(),
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: pump_backward,
            },
        };
        let mut solver = PicardSolver::filled(self.grid.points(), initial);
        let residual = |sgnl_b| -> Result<f64, SolverError> {
            let set_boundary = |current: &[FieldState]| FieldState {
                signal: BidirectionalAmplitude {
                    forward: 0.0,
                    backward: sgnl_b,
                },
                pump: BidirectionalAmplitude {
                    forward: pump_forward,
                    backward: self.lasing_pump_backward(pump_backward, current),
                },
            };
            self.solve_picard(&mut solver, picard_config, set_boundary)?;
            Ok(solver.profile().last().unwrap().signal.backward / sgnl_b)
        };
        let _sgnl_b = try_rootfind_1d(residual, root_find_config)?;

        Ok(self.field_profile(
            self.grid.positions().collect(),
            solver.profile().to_vec(),
        ))
    }

    fn lasing_pump_backward(&self, pump_backward: f64, profile: &[FieldState]) -> f64 {
        let dz = self.grid.dz();
        let optical_depth: f64 = profile[..profile.len() - 1]
            .iter()
            .map(|&field| 0.5 * self.fibre.gain(field).pump * dz)
            .sum();
        pump_backward * optical_depth.exp()
    }

    fn injected_left_boundary(
        &self,
        pump: Pump,
        signal: Signal,
        profile: &[FieldState],
    ) -> FieldState {
        let fibre = self.fibre;
        let grid = self.grid;
        assert_eq!(profile.len(), grid.points());

        let dz = grid.dz();
        let mut pump_optical_depth = 0.0;
        let mut signal_transfer_bottom_row = (0.0, 1.0); // calculate bottom row of transfer matrix

        for (&field, &kappa) in profile[..grid.steps()].iter().zip(self.kappas()).rev() {
            let gain = fibre.gain(field);
            pump_optical_depth += 0.5 * gain.pump * dz;
            signal_transfer_bottom_row = update_signal_transfer_bottom_row(
                signal_transfer_bottom_row,
                gain.signal,
                kappa,
                dz,
            );
        }

        let (signal_forward, signal_backward_right) =
            self.fibre.signal_flux_amplitudes(signal);
        let (pump_forward, pump_backward_right) = self.fibre.pump_flux_amplitudes(pump);
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
    use crate::grating::{NoGrating, PiShift, sample_grating};
    use crate::maths::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
    use crate::maths::utils::IterationConfig;
    use crate::two_mode::profile_max_diff;

    const LASING_STEPS: usize = 500;
    const LASING_PUMP_AMPLITUDE: f64 = 100.0;
    const LASING_PUMP_FLUX: f64 = LASING_PUMP_AMPLITUDE * LASING_PUMP_AMPLITUDE;
    const LASING_ITERATION: IterationConfig = IterationConfig {
        max: 500,
        tol: 1e-10,
    };
    const LASING_PICARD: PicardConfig = PicardConfig {
        max_iterations: 500,
        relative_tolerance: 1e-10,
        absolute_tolerance: 1e-12,
    };
    const LASING_BISECTION: BisectionConfig = BisectionConfig {
        iteration: LASING_ITERATION,
        upper: LASING_PUMP_AMPLITUDE,
        lower: 1e-8,
        midpoint: Midpoint::Geometric,
    };
    const LASING_NEWTON: Newton1dConfig = Newton1dConfig {
        iteration: LASING_ITERATION,
        initial: LASING_PUMP_AMPLITUDE,
        dx: 1e-6,
    };

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
    ) -> ResolvedFibre<'_, TwoLevelDopant, G> {
        fibre.resolve_with_interactions(
            FieldMode::new(970e-9),
            TwoLevelCrossSections::new(1.0, 0.0),
            FieldMode::new(1060e-9),
            TwoLevelCrossSections::new(0.0, 1.0),
        )
    }

    fn active_lasing_fibre() -> Fibre<TwoLevelDopant, PiShift> {
        Fibre {
            geometry: FibreGeometry {
                core_radius: 4e-6,
                numerical_aperture: 0.1,
                length: 10.0,
            },
            dopant: TwoLevelDopant {
                density: 1.0,
                lifetime: 1.0,
            },
            grating: PiShift {
                kappa_left: 1.0,
                kappa_right: 1.0,
                pi_shift_position: 0.45,
            },
        }
    }

    fn resolve_active_lasing(
        fibre: &Fibre<TwoLevelDopant, PiShift>,
    ) -> ResolvedFibre<'_, TwoLevelDopant, PiShift> {
        fibre.resolve_with_interactions(
            FieldMode::new(970e-9),
            TwoLevelCrossSections::new(1.0, 0.0),
            FieldMode::new(1060e-9),
            TwoLevelCrossSections::new(0.0, 1.0),
        )
    }

    fn assert_near(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
    }

    fn pump_for_flux<D: DopantModel, G: GratingModel>(
        fibre: &ResolvedFibre<'_, D, G>,
        total_flux: f64,
        balance: f64,
    ) -> Pump {
        Pump {
            total: fibre.pump_power(total_flux),
            balance,
        }
    }

    fn assert_profiles_identical(left: &FieldProfile, right: &FieldProfile) {
        assert_eq!(left.z, right.z);
        assert_eq!(left.fields.len(), right.fields.len());
        assert_eq!(profile_max_diff(&left.fields, &right.fields), 0.0);
    }

    fn assert_lasing_algorithms_agree(root_find_config: RootFindConfig) {
        let fibre = active_lasing_fibre();
        let fibre = resolve_active_lasing(&fibre);
        let solver = TwoModeSolver::new(&fibre, LASING_STEPS);
        let pump = pump_for_flux(&fibre, LASING_PUMP_FLUX, 1.0);
        let shooting = solver
            .solve_lasing_shooting(pump, root_find_config)
            .expect("shooting lasing solve failed");
        let picard = solver
            .solve_lasing_picard(pump, root_find_config, LASING_PICARD)
            .expect("Picard lasing solve failed");

        assert_profiles_identical(&shooting, &picard);
    }

    #[test]
    fn reduced_lasing_profile_preserves_full_profile_endpoints() {
        let fibre = active_lasing_fibre();
        let fibre = resolve_active_lasing(&fibre);
        let solver = TwoModeSolver::new(&fibre, LASING_STEPS);
        let pump = pump_for_flux(&fibre, LASING_PUMP_FLUX, 1.0);
        let full = solver
            .solve_lasing(
                pump,
                RootFindConfig::Newton1d(LASING_NEWTON),
                LASING_PICARD,
                true,
            )
            .expect("full lasing solve failed");
        let reduced = solver
            .solve_lasing(
                pump,
                RootFindConfig::Newton1d(LASING_NEWTON),
                LASING_PICARD,
                false,
            )
            .expect("reduced lasing solve failed");
        let expected_fields = [full.fields[0], *full.fields.last().unwrap()];

        assert_eq!(reduced.z, [full.z[0], *full.z.last().unwrap()]);
        assert_eq!(reduced.len(), 2);
        assert_eq!(profile_max_diff(&reduced.fields, &expected_fields), 0.0);
        assert_eq!(reduced.output_powers(), full.output_powers());
    }

    #[test]
    fn direct_and_picard_lasing_profiles_agree() {
        let fibre = active_lasing_fibre();
        let fibre = resolve_active_lasing(&fibre);
        let solver = TwoModeSolver::new(&fibre, LASING_STEPS);
        let sgnl_b = 1.0;
        let pump = pump_for_flux(&fibre, LASING_PUMP_FLUX, 1.0);
        let pump_forward = fibre.pump_flux_amplitudes(pump).0;
        let boundary = FieldState {
            signal: BidirectionalAmplitude {
                forward: 0.0,
                backward: sgnl_b,
            },
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: 0.0,
            },
        };
        let direct = solver.field_profile(
            solver.grid.positions().collect(),
            solve_profile_coupled(
                boundary,
                |fields| solver.fibre.gain(fields),
                solver.grid.dz(),
                solver.kappas(),
            ),
        );
        let mut picard_solver = PicardSolver::filled(solver.grid.points(), boundary);
        solver
            .solve_picard(
                &mut picard_solver,
                LASING_PICARD,
                |_: &[FieldState]| boundary,
            )
            .expect("Picard profile solve failed");
        let picard = solver.field_profile(
            direct.z.clone(),
            picard_solver.profile().to_vec(),
        );

        assert_profiles_identical(&direct, &picard);
    }

    #[test]
    fn shooting_and_picard_lasing_solvers_agree_with_newton() {
        assert_lasing_algorithms_agree(RootFindConfig::Newton1d(LASING_NEWTON));
    }

    #[test]
    fn shooting_and_picard_lasing_solvers_agree_with_bisection() {
        assert_lasing_algorithms_agree(RootFindConfig::Bisection(LASING_BISECTION));
    }

    #[test]
    fn no_grating_solver_owns_fixed_zero_profile() {
        let fibre = zero_gain_fibre(NoGrating);
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);

        assert_eq!(solver.grid().steps(), 4);
        assert_eq!(solver.grid().points(), 5);
        assert_eq!(solver.grid().position(4), fibre.length());
        assert_eq!(solver.kappas(), &[0.0; 4]);
    }

    #[test]
    fn pi_shift_solver_caches_left_edge_samples() {
        let grating = PiShift {
            kappa_left: 2.0,
            kappa_right: 3.0,
            pi_shift_position: 0.5,
        };
        let fibre = zero_gain_fibre(grating);
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);

        assert_eq!(solver.kappas().len(), solver.grid().steps());
        assert_eq!(solver.kappas(), sample_grating(&grating, 4));
    }

    #[test]
    fn injected_boundary_reduces_to_inputs_without_gain_or_grating() {
        let fibre = zero_gain_fibre(NoGrating);
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);
        let signal = Signal {
            total: 5.0,
            balance: 0.2,
        };
        let pump = Pump {
            total: 7.0,
            balance: -0.4,
        };
        let profile = vec![FieldState::default(); solver.grid().points()];

        let boundary = solver.injected_left_boundary(pump, signal, &profile);
        let (signal_forward, signal_backward) = fibre.signal_flux_amplitudes(signal);
        let (pump_forward, pump_backward) = fibre.pump_flux_amplitudes(pump);

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
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);
        let signal = Signal {
            total: 5.0,
            balance: 0.2,
        };
        let pump = Pump {
            total: 7.0,
            balance: -0.4,
        };
        let profile = vec![FieldState::default(); solver.grid().points()];

        let boundary = solver.injected_left_boundary(pump, signal, &profile);
        let mut propagated_signal = boundary.signal;
        for &kappa in solver.kappas() {
            propagated_signal = propagated_signal.coupled_step(0.0, kappa, solver.grid().dz());
        }

        assert_near(
            propagated_signal.backward,
            fibre.signal_flux_amplitudes(signal).1,
        );
        assert_near(
            boundary.pump.backward,
            fibre.pump_flux_amplitudes(pump).1,
        );
    }

    #[test]
    fn forward_injected_grating_solve_satisfies_both_signal_boundaries() {
        let fibre = zero_gain_fibre(PiShift {
            kappa_left: 0.7,
            kappa_right: 0.3,
            pi_shift_position: 0.5,
        });
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);
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
        let (signal_forward, signal_backward) = fibre.signal_flux_amplitudes(signal);
        let (pump_forward, pump_backward) = fibre.pump_flux_amplitudes(pump);

        assert_near(left.signal.forward, signal_forward);
        assert_near(right.signal.backward, signal_backward);
        assert_near(left.pump.forward, pump_forward);
        assert_near(right.pump.backward, pump_backward);
        assert!(
            left.signal.backward.abs() > 1e-6,
            "grating solve should produce a reflected signal"
        );
    }

    #[test]
    fn forward_injected_no_grating_solve_preserves_zero_gain_fields() {
        let fibre = zero_gain_fibre(NoGrating);
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);
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
                forward: fibre.signal_flux_amplitudes(signal).0,
                backward: 0.0,
            },
            pump: BidirectionalAmplitude {
                forward: fibre.pump_flux_amplitudes(pump).0,
                backward: 0.0,
            },
        };

        assert_eq!(profile.fields.len(), solver.grid().points());
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
        let fibre = resolve_zero_gain(&fibre);
        let solver = TwoModeSolver::new(&fibre, 4);
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
        let (signal_forward, signal_backward) = fibre.signal_flux_amplitudes(signal);
        let (pump_forward, pump_backward) = fibre.pump_flux_amplitudes(pump);

        assert_near(left.signal.forward, signal_forward);
        assert_near(right.signal.backward, signal_backward);
        assert_near(left.pump.forward, pump_forward);
        assert_near(right.pump.backward, pump_backward);
    }
}
