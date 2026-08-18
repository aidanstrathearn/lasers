use crate::lase::{
    FibreParams, FieldProfile, FieldState, GridPoints, Pump, gain, profile_max_diff,
};
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
pub fn find_pump_b(pump: Pump, profile: &FieldProfile, fp: FibreParams, dz: f64) -> f64 {
    let expg: f64 = profile.fields[..profile.fields.len() - 1]
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

pub fn solve_profile_picard(
    sgnl_b: f64,
    pump: Pump,
    fp: FibreParams,
    gp: GridPoints,
    ic: IterationConfig,
    kappas: &[f64],
) -> Result<FieldProfile, PicardError> {
    let dz = gp.dz(fp.length);
    let mut current = initial_profile(pump, fp, gp);
    assert_eq!(kappas.len() + 1, current.fields.len());
    let mut new_fields = current.fields.clone();
    let boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b: sgnl_b,
        pump_f: pump.forward,
        pump_b: 0.0,
    };
    for _ in 0..ic.max {
        new_fields[0] = FieldState {
            pump_b: find_pump_b(pump, &current, fp, dz),
            ..boundary
        };
        for j in 1..new_fields.len() {
            new_fields[j] =
                new_fields[j - 1].general_step(current.fields[j - 1], fp, kappas[j - 1], dz);
        }
        let diff = profile_max_diff(&current.fields, &new_fields);
        current.fields = new_fields.clone();
        if diff < ic.tol {
            return Ok(current);
        }
    }
    Err(PicardError::DidNotConverge)
}
