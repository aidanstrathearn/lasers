use crate::error::SolverError;
use crate::fibre::BidirectionalAmplitude;
use crate::grating::GratingModel;
use crate::lase::{DopantModel, FieldState};
use crate::maths::rootfind::{RootFindConfig, RootFindError, rootfind_1d, try_rootfind_1d};
use crate::two_mode::propagation::out_field_coupled;
use crate::two_mode::{Pump, TwoModeSolver};

impl<D: DopantModel, G: GratingModel> TwoModeSolver<'_, D, G> {
    /// Finds the pump power in watts where the small-signal lasing residual is zero.
    ///
    /// The threshold root variable is total pump power in watts. For bidirectional pumping, the
    /// pump-boundary root variable is the ratio between the backward pump amplitudes at the left
    /// and right boundaries.
    pub fn find_threshold(
        &self,
        pump_balance: f64,
        threshold_root_find_config: RootFindConfig,
        pump_boundary_root_find_config: RootFindConfig,
    ) -> Result<f64, SolverError> {
        let (_, pump_backward) = self.fibre().pump_flux_amplitudes(Pump {
            total: 1.0,
            balance: pump_balance,
        });

        if pump_backward == 0.0 {
            rootfind_1d(
                |total| {
                    self.threshold_residual_shooting(Pump {
                        total,
                        balance: pump_balance,
                    })
                },
                threshold_root_find_config,
            )
            .map_err(SolverError::from)
        } else {
            // this one works using nested rootfinds
            // inner one finds correct pump boundary
            // outer one finds threshold
            try_rootfind_1d(
                |total| {
                    self.threshold_residual_bidirectional(
                        Pump {
                            total,
                            balance: pump_balance,
                        },
                        pump_boundary_root_find_config,
                    )
                },
                threshold_root_find_config,
            )
            .map_err(SolverError::from)
        }
    }

    fn threshold_residual_shooting(&self, pump: Pump) -> f64 {
        let (pump_forward, pump_backward) = self.fibre().pump_flux_amplitudes(pump);
        assert_eq!(
            pump_backward, 0.0,
            "shooting threshold requires a forward-only pump"
        );
        self.threshold_residual_from_left_pump(pump_forward, 0.0)
    }

    fn threshold_residual_bidirectional(
        &self,
        pump: Pump,
        pump_boundary_root_find_config: RootFindConfig,
    ) -> Result<f64, RootFindError> {
        let (pump_forward, target_pump_backward) = self.fibre().pump_flux_amplitudes(pump);
        assert!(
            target_pump_backward > 0.0,
            "bidirectional threshold requires a backward pump"
        );
        let pump_boundary_residual = |left_ratio| {
            let initial = FieldState {
                signal: BidirectionalAmplitude::default(),
                pump: BidirectionalAmplitude {
                    forward: pump_forward,
                    backward: left_ratio * target_pump_backward,
                },
            };
            let output = out_field_coupled(
                initial,
                |fields| self.fibre().gain(fields),
                self.grid().dz(),
                self.kappas(),
            );
            output.pump.backward / target_pump_backward - 1.0
        };
        let left_ratio = rootfind_1d(pump_boundary_residual, pump_boundary_root_find_config)?;

        Ok(self.threshold_residual_from_left_pump(pump_forward, left_ratio * target_pump_backward))
    }

    fn threshold_residual_from_left_pump(&self, pump_forward: f64, pump_backward: f64) -> f64 {
        let initial = FieldState {
            signal: BidirectionalAmplitude {
                forward: 0.0,
                backward: 1.0,
            },
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: pump_backward,
            },
        };

        out_field_coupled(
            initial,
            |fields| {
                self.fibre().gain(FieldState {
                    signal: BidirectionalAmplitude::default(),
                    ..fields
                })
            },
            self.grid().dz(),
            self.kappas(),
        )
        .signal
        .backward
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dopant::{TwoLevelCrossSections, TwoLevelDopant};
    use crate::fibre::{Fibre, FibreGeometry, FieldMode};
    use crate::grating::PiShift;
    use crate::maths::rootfind::{BisectionConfig, Midpoint};
    use crate::maths::utils::IterationConfig;
    use crate::two_mode::ResolvedFibre;

    const STEPS: usize = 500;
    const PUMP_AMPLITUDE: f64 = 100.0;
    const PUMP_FLUX: f64 = PUMP_AMPLITUDE * PUMP_AMPLITUDE;
    const ITERATION: IterationConfig = IterationConfig {
        max: 500,
        tol: 1e-10,
    };

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

    fn threshold_configs(
        fibre: &ResolvedFibre<'_, TwoLevelDopant, PiShift>,
    ) -> (f64, f64, BisectionConfig, BisectionConfig) {
        let lower = fibre.pump_power(1e-6);
        let upper = fibre.pump_power(PUMP_FLUX);
        let threshold = BisectionConfig {
            iteration: ITERATION,
            lower,
            upper,
            midpoint: Midpoint::Geometric,
        };
        let pump_boundary = BisectionConfig {
            iteration: ITERATION,
            lower: 1e-8,
            upper: 1.0,
            midpoint: Midpoint::Geometric,
        };
        (lower, upper, threshold, pump_boundary)
    }

    #[test]
    fn shooting_threshold_zeros_small_signal_residual() {
        let fibre = active_lasing_fibre();
        let fibre = resolve_active_lasing(&fibre);
        let solver = TwoModeSolver::new(&fibre, STEPS);
        let (lower, upper, threshold_config, pump_boundary_config) = threshold_configs(&fibre);

        let threshold = solver
            .find_threshold(1.0, threshold_config.into(), pump_boundary_config.into())
            .expect("shooting threshold solve failed");

        assert!(threshold > lower && threshold < upper);
        assert!(
            solver
                .threshold_residual_shooting(Pump {
                    total: threshold,
                    balance: 1.0,
                })
                .abs()
                < ITERATION.tol
        );
    }

    #[test]
    fn bidirectional_threshold_zeros_nested_shooting_residual() {
        let fibre = active_lasing_fibre();
        let fibre = resolve_active_lasing(&fibre);
        let solver = TwoModeSolver::new(&fibre, STEPS);
        let (lower, upper, threshold_config, pump_boundary_config) = threshold_configs(&fibre);
        let balance = 0.95;

        let threshold = solver
            .find_threshold(
                balance,
                threshold_config.into(),
                pump_boundary_config.into(),
            )
            .expect("bidirectional threshold solve failed");
        let residual = solver
            .threshold_residual_bidirectional(
                Pump {
                    total: threshold,
                    balance,
                },
                pump_boundary_config.into(),
            )
            .expect("backward pump boundary solve failed");

        assert!(threshold > lower && threshold < upper);
        assert!(residual.abs() < ITERATION.tol);
    }
}
