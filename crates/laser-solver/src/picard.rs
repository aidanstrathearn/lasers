use crate::lase::FieldState;
use std::fmt;

#[derive(Debug)]
pub enum PicardError {
    DidNotConverge,
}

impl fmt::Display for PicardError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DidNotConverge => {
                write!(formatter, "picard iteration did not converge")
            }
        }
    }
}

impl std::error::Error for PicardError {}

#[derive(Copy, Clone, Debug)]
pub struct PicardConfig {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
}

impl Default for PicardConfig {
    fn default() -> Self {
        Self {
            max_iterations: 500,
            relative_tolerance: 1e-8,
            absolute_tolerance: 1e-12,
        }
    }
}

pub struct PicardSolver {
    current: Vec<FieldState>,
    new: Vec<FieldState>,
}

impl PicardSolver {
    pub fn from_initial(initial: Vec<FieldState>) -> Self {
        let n = initial.len();
        Self {
            current: initial,
            new: vec![FieldState::default(); n],
        }
    }

    pub fn profile(&self) -> &[FieldState] {
        &self.current
    }

    pub fn solve<SetBoundary, Step>(
        &mut self,
        config: PicardConfig,
        mut set_boundary: SetBoundary,
        mut step: Step,
    ) -> Result<&[FieldState], PicardError>
    where
        SetBoundary: FnMut(&[FieldState]) -> FieldState,
        Step: FnMut(FieldState, FieldState, usize) -> FieldState,
    {
        for _ in 0..config.max_iterations {
            self.new[0] = set_boundary(&self.current);
            for i in 0..self.new.len() - 1 {
                self.new[i + 1] = step(self.new[i], self.current[i], i);
            }

            std::mem::swap(&mut self.current, &mut self.new);
            let error = profile_convergence_error(&self.current, &self.new, config);
            if error <= 1.0 {
                return Ok(&self.current);
            }
        }
        Err(PicardError::DidNotConverge)
    }
}

pub fn profile_convergence_error(
    current: &[FieldState],
    new: &[FieldState],
    config: PicardConfig,
) -> f64 {
    assert_eq!(current.len(), new.len());
    let mut max_dif_s = 0.0_f64;
    let mut max_dif_p = 0.0_f64;
    let mut max_mag_s = 0.0_f64;
    let mut max_mag_p = 0.0_f64;

    for (&current, &new) in current.iter().zip(new) {
        let current_powers = current.field_powers();
        let new_powers = new.field_powers();
        if !current_powers[0].is_finite()
            || !new_powers[0].is_finite()
            || !current_powers[1].is_finite()
            || !new_powers[1].is_finite()
        {
            return f64::INFINITY;
        }
        max_dif_s = max_dif_s.max((current_powers[0] - new_powers[0]).abs().sqrt());
        max_dif_p = max_dif_p.max((current_powers[1] - new_powers[1]).abs().sqrt());
        max_mag_s = max_mag_s.max(current_powers[0].max(new_powers[0]).sqrt());
        max_mag_p = max_mag_p.max(current_powers[1].max(new_powers[1]).sqrt());
    }
    (max_dif_p / (config.absolute_tolerance + config.relative_tolerance * max_mag_p))
        .max(max_dif_s / (config.absolute_tolerance + config.relative_tolerance * max_mag_s))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn convergence_check_config() -> PicardConfig {
        PicardConfig {
            max_iterations: 1,
            relative_tolerance: 1e-6,
            absolute_tolerance: 1e-10,
        }
    }

    #[test]
    fn step_receives_new_previous_old_current_and_interval_index() {
        let initial = vec![
            FieldState {
                pump_f: 1.0,
                ..FieldState::default()
            },
            FieldState {
                pump_f: 2.0,
                ..FieldState::default()
            },
            FieldState {
                pump_f: 3.0,
                ..FieldState::default()
            },
        ];
        let mut solver = PicardSolver::from_initial(initial);
        let mut steps = Vec::new();
        let config = PicardConfig {
            max_iterations: 1,
            relative_tolerance: 0.0,
            absolute_tolerance: 100.0,
        };

        solver
            .solve(
                config,
                |current| {
                    assert_eq!(
                        current.iter().map(|field| field.pump_f).collect::<Vec<_>>(),
                        vec![1.0, 2.0, 3.0]
                    );
                    FieldState {
                        pump_f: 10.0,
                        ..FieldState::default()
                    }
                },
                |new_previous, old_current, i| {
                    steps.push((new_previous.pump_f, old_current.pump_f, i));
                    FieldState {
                        pump_f: new_previous.pump_f + 1.0,
                        ..FieldState::default()
                    }
                },
            )
            .expect("single Picard iteration should satisfy the loose tolerance");

        assert_eq!(steps, vec![(10.0, 1.0, 0), (11.0, 2.0, 1)]);
        assert_eq!(
            solver
                .profile()
                .iter()
                .map(|field| field.pump_f)
                .collect::<Vec<_>>(),
            vec![10.0, 11.0, 12.0]
        );
    }

    #[test]
    fn identical_profiles_have_zero_error() {
        let profile = vec![FieldState {
            pump_f: 100.0,
            pump_b: 10.0,
            sgnl_f: 1.0,
            sgnl_b: -1.0,
        }];

        let error = profile_convergence_error(&profile, &profile, convergence_check_config());

        assert_eq!(error, 0.0);
    }

    #[test]
    fn tiny_zero_crossing_converges() {
        let current = vec![FieldState {
            sgnl_b: 1e-12,
            ..FieldState::default()
        }];
        let new = vec![FieldState {
            sgnl_b: -1e-12,
            ..FieldState::default()
        }];

        let error = profile_convergence_error(&current, &new, convergence_check_config());

        assert!(error <= 1.0, "tiny zero crossing error was {error:e}");
    }

    #[test]
    fn localized_error_is_not_hidden() {
        let current = vec![
            FieldState {
                pump_f: 1.0,
                ..FieldState::default()
            };
            100
        ];
        let mut new = current.clone();
        new[50].pump_f = 1.01;

        let error = profile_convergence_error(&current, &new, convergence_check_config());

        assert!(error > 1.0, "localized profile error was {error:e}");
    }

    #[test]
    fn fields_are_scaled_independently() {
        let current = vec![FieldState {
            pump_f: 1e6,
            sgnl_b: 1e-6,
            ..FieldState::default()
        }];
        let new = vec![FieldState {
            pump_f: 1e6,
            sgnl_b: 2e-6,
            ..FieldState::default()
        }];

        let error = profile_convergence_error(&current, &new, convergence_check_config());

        assert!(error > 1.0, "signal error was hidden by pump scale");
    }

    #[test]
    fn non_finite_values_fail_convergence() {
        let current = vec![FieldState::default()];
        let new = vec![FieldState {
            pump_b: f64::NAN,
            ..FieldState::default()
        }];

        let error = profile_convergence_error(&current, &new, convergence_check_config());

        assert!(error.is_infinite());
    }
}
