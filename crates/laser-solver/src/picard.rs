use crate::error::SolverError;
use crate::lase::{FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump, gain};
use crate::rootfind::{RootFindConfig, try_rootfind_1d};
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

pub struct PicardDfbSolver {
    initial: Vec<FieldState>,
    current: Vec<FieldState>,
    new: Vec<FieldState>,
}

impl PicardDfbSolver {
    pub fn from_initial(initial: Vec<FieldState>) -> Self {
        let n = initial.len();
        Self {
            initial,
            current: vec![FieldState::default(); n],
            new: vec![FieldState::default(); n],
        }
    }

    pub fn new(pump: Pump, fp: FibreParams, gp: GridPoints) -> Self {
        Self::from_initial(initial_profile(pump, fp, gp).fields)
    }

    pub fn solve_profile_picard(
        &mut self,
        sgnl_b: f64,
        pump: Pump,
        fp: FibreParams,
        config: PicardConfig,
        kappas: &[f64],
        dz: f64,
    ) -> Result<&[FieldState], PicardError> {
        assert_eq!(kappas.len() + 1, self.initial.len());

        self.current.copy_from_slice(&self.initial);
        let boundary = FieldState {
            sgnl_f: 0.0,
            sgnl_b: sgnl_b,
            pump_f: pump.forward,
            pump_b: 0.0,
        };
        for _ in 0..config.max_iterations {
            self.new[0] = FieldState {
                pump_b: find_pump_b(pump, &self.current, fp, dz),
                ..boundary
            };
            for j in 1..self.new.len() {
                self.new[j] =
                    self.new[j - 1].general_step(self.current[j - 1], fp, kappas[j - 1], dz);
            }
            let error = profile_convergence_error(&self.current, &self.new, config);
            //current = new.clone();
            std::mem::swap(&mut self.current, &mut self.new);
            if error <= 1.0 {
                self.initial.copy_from_slice(&self.current);
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

pub fn initial_profile(pump: Pump, fp: FibreParams, gp: GridPoints) -> FieldProfile {
    let g = -fp.pump_ab * fp.density;
    let zs = gp.grid(fp.length);
    let end_factor = (0.5 * g * fp.length).exp();

    let fields = zs
        .iter()
        .map(|z| {
            let f = (0.5 * g * z).exp(); // &f64 * f64 -> f64 apparently, so no need to deref
            let b = end_factor / f;

            FieldState {
                sgnl_f: 0.0,
                sgnl_b: 0.0,
                pump_f: f * pump.forward,
                pump_b: b * pump.backward,
            }
        })
        .collect();
    FieldProfile::new(zs, fields)
}
pub fn find_pump_b(pump: Pump, profile: &Vec<FieldState>, fp: FibreParams, dz: f64) -> f64 {
    let expg: f64 = profile[..profile.len() - 1]
        .iter()
        .map(|&field| {
            let (g, _) = gain(field, fp);
            0.5 * g * dz
        })
        .sum::<f64>() // dont know why it couldnt infer f64 here
        .exp();
    pump.backward * expg
}

pub fn dfb_solve_from_picard_solver(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    solver: &mut PicardDfbSolver,
    picard_config: PicardConfig,
) -> Result<FieldProfile, SolverError> {
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let f = |sgnl_b| -> Result<f64, SolverError> {
        let fields = solver.solve_profile_picard(sgnl_b, pu, fp, picard_config, &kappas, dz)?;
        Ok(fields.last().unwrap().sgnl_b)
    };
    // try_rootfind_1d muts the solver which leaves the lasing solution in the 'current' buffer
    let _sgnl_b = try_rootfind_1d(f, config)?;
    if full_profile {
        Ok(FieldProfile::new(
            gp.grid(fp.length),
            solver.current.clone(),
        ))
    } else {
        Ok(FieldProfile::new(
            vec![0.0_f64, fp.length],
            vec![solver.current[0], solver.current.last().copied().unwrap()],
        ))
    }
}

pub fn dfb_solve_picard(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    picard_config: PicardConfig,
) -> Result<FieldProfile, SolverError> {
    let mut solver = PicardDfbSolver::new(pu, fp, gp);
    dfb_solve_from_picard_solver(
        pu,
        fp,
        gp,
        kp,
        full_profile,
        config,
        &mut solver,
        picard_config,
    )
}

pub fn dfb_pump_scan_picard(
    pumps: &[f64],
    balance: f64,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
    picard_config: PicardConfig,
) -> Vec<(f64, f64, bool)> {
    let full_profile = false;
    let mut solver = PicardDfbSolver::new(Pump::from_total_and_balance(pumps[0], balance), fp, gp);
    pumps
        .iter()
        .map(|&pmp| {
            let pu = Pump::from_total_and_balance(pmp, balance);

            let result = dfb_solve_from_picard_solver(
                pu,
                fp,
                gp,
                kp,
                full_profile,
                config,
                &mut solver,
                picard_config,
            )
            .map_or((0.0, 0.0, false), |result| {
                (
                    result.sgnl_f().last().unwrap(),
                    result.sgnl_b().next().unwrap(),
                    true,
                )
            });
            result
        })
        .collect()
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
