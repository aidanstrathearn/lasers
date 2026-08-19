use crate::dfb::{out_field, solve_profile};
use crate::error::SolverError;
use crate::lase::{gain, profile_avg_diff, profile_max_diff, FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump};
use crate::rootfind::{RootFindConfig, rootfind_1d, try_rootfind_1d};
use crate::utils::IterationConfig;
use std::fmt;

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

pub struct PicardDfbSolver {
    initial: Vec<FieldState>,
    current: Vec<FieldState>,
    new: Vec<FieldState>,
}

impl PicardDfbSolver {
    pub fn init(initial: Vec<FieldState>) -> Self {
        let n = initial.len();
        Self {
            initial,
            current: vec![FieldState::default(); n],
            new: vec![FieldState::default(); n],
        }
    }

    pub fn solve_profile_picard(
        &mut self,
        sgnl_b: f64,
        pump: Pump,
        fp: FibreParams,
        ic: IterationConfig,
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
        for _ in 0..ic.max {
            self.new[0] = FieldState {
                pump_b: find_pump_b(pump, &self.current, fp, dz),
                ..boundary
            };
            for j in 1..self.new.len() {
                self.new[j] =
                    self.new[j - 1].general_step(self.current[j - 1], fp, kappas[j - 1], dz);
            }
            let diff = profile_max_diff(&self.current, &self.new);
            //current = new.clone();
            std::mem::swap(&mut self.current, &mut self.new);
            if diff < ic.tol {
                self.initial.copy_from_slice(&self.current);
                return Ok(&self.current);
            }
        }
        Err(PicardError::DidNotConverge)
    }
}

pub fn solve_profile_picard(
    sgnl_b: f64,
    initial: Vec<FieldState>,
    pump: Pump,
    fp: FibreParams,
    ic: IterationConfig,
    kappas: &[f64],
    dz: f64,
) -> Result<Vec<FieldState>, PicardError> {
    assert_eq!(kappas.len() + 1, initial.len());

    let mut current = initial;
    let mut new = current.clone();
    let boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b: sgnl_b,
        pump_f: pump.forward,
        pump_b: 0.0,
    };
    for _ in 0..ic.max {
        new[0] = FieldState {
            pump_b: find_pump_b(pump, &current, fp, dz),
            ..boundary
        };
        for j in 1..new.len() {
            new[j] = new[j - 1].general_step(current[j - 1], fp, kappas[j - 1], dz);
        }
        let diff = profile_avg_diff(&current, &new);
        //current = new.clone();
        std::mem::swap(&mut current, &mut new);
        if diff < ic.tol {
            return Ok(current);
        }
    }
    Err(PicardError::DidNotConverge)
}

pub fn dfb_solve_picard(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    ic: IterationConfig,
) -> Result<FieldProfile, SolverError> {
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let initial = initial_profile(pu, fp, gp);
    let z = initial.z().collect();
    let initial_fields = initial.fields;

    let mut f = |sgnl_b| -> Result<f64, SolverError> {
        let fields = solve_profile_picard(sgnl_b, initial_fields.clone(), pu, fp, ic, &kappas, dz)?;
        Ok(fields.last().unwrap().sgnl_b)
    };
    let sgnl_b = try_rootfind_1d(&mut f, config)?;
    let fields = solve_profile_picard(sgnl_b, initial_fields, pu, fp, ic, &kappas, dz)?;
    if full_profile {
        Ok(FieldProfile::new(z, fields))
    } else {
        let z = vec![0.0_f64, fp.length];
        let fields = vec![fields[0], fields.last().copied().unwrap()];
        Ok(FieldProfile::new(z, fields))
    }
}

pub fn dfb_solve_picard_buffers(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
    ic: IterationConfig,
) -> Result<FieldProfile, SolverError> {
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let initial = initial_profile(pu, fp, gp);
    let z = initial.z().collect();
    let initial_fields = initial.fields;
    let mut solver = PicardDfbSolver::init(initial_fields.clone());

    let f = |sgnl_b| -> Result<f64, SolverError> {
        let fields = solver.solve_profile_picard(sgnl_b, pu, fp, ic, &kappas, dz)?;
        Ok(fields.last().unwrap().sgnl_b)
    };
    let sgnl_b = try_rootfind_1d(f, config)?;
    let fields = solve_profile_picard(sgnl_b, initial_fields, pu, fp, ic, &kappas, dz)?;
    if full_profile {
        Ok(FieldProfile::new(z, fields))
    } else {
        let z = vec![0.0_f64, fp.length];
        let fields = vec![fields[0], fields.last().copied().unwrap()];
        Ok(FieldProfile::new(z, fields))
    }
}
