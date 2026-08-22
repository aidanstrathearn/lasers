use crate::error::SolverError;
use crate::utils::{linspace, relative_diff, IterationConfig};

#[derive(Copy, Clone)]
pub struct FibreParams {
    pub density: f64,
    pub lifetime: f64,
    pub pump_ab: f64,
    pub pump_em: f64,
    pub sgnl_ab: f64,
    pub sgnl_em: f64,
    pub length: f64,
}

impl Default for FibreParams {
    fn default() -> Self {
        Self {
            density: 1.0,
            lifetime: 1.0,
            pump_ab: 0.01,
            pump_em: 0.0,
            sgnl_ab: 0.0,
            sgnl_em: 1.0,
            length: 5.0,
        }
    }
}

#[allow(dead_code)]
#[derive(Copy, Clone)]
pub struct Reflectivities {
    left: f64,
    right: f64,
}

#[derive(Copy, Clone)]
pub struct GridPoints(pub usize);

impl Default for GridPoints {
    fn default() -> Self {
        Self(100)
    }
}

impl GridPoints {
    pub fn grid(self, length: f64) -> Vec<f64> {
        linspace(0.0, length, self.0)
    }

    pub fn dz(self, length: f64) -> f64 {
        length / self.0 as f64
    }
}

#[derive(Copy, Clone)]
pub struct GratingProfile {
    pub kappa_left: f64,
    pub kappa_right: f64,
    pub pi_shift_position: f64,
}

impl Default for GratingProfile {
    fn default() -> Self {
        Self {
            kappa_left: 1.0,
            kappa_right: 1.0,
            pi_shift_position: 0.45,
        }
    }
}

impl GratingProfile {
    pub fn grid(self, n: usize) -> Vec<f64> {
        (0..n)
            .map(|j| {
                let z = j as f64 / n as f64;
                if z < self.pi_shift_position {
                    self.kappa_left
                } else {
                    -self.kappa_right
                }
            })
            .collect()
    }
}

#[derive(Copy, Clone, Debug, Default)]
pub struct FieldState {
    pub sgnl_f: f64,
    pub sgnl_b: f64,
    pub pump_f: f64,
    pub pump_b: f64,
}

impl FieldState {
    pub(crate) fn field_powers(self) -> [f64; 2] {
        [
            self.sgnl_f * self.sgnl_f + self.sgnl_b * self.sgnl_b,
            self.pump_f * self.pump_f + self.pump_b * self.pump_b,
        ]
    }
}

pub fn field_max_diff(f1: FieldState, f2: FieldState) -> f64 {
    let diffs = [
        relative_diff(f1.pump_f, f2.pump_f),
        relative_diff(f1.pump_b, f2.pump_b),
        relative_diff(f1.sgnl_f, f2.sgnl_f),
        relative_diff(f1.sgnl_b, f2.sgnl_b),
    ];
    diffs.into_iter().fold(0.0, f64::max)
}
pub fn profile_max_diff(p1: &Vec<FieldState>, p2: &Vec<FieldState>) -> f64 {
    p1.iter()
        .zip(p2.iter())
        .map(|(&f1, &f2)| field_max_diff(f1, f2))
        .reduce(f64::max)
        .unwrap_or(f64::NAN)
}

pub fn field_avg_diff(f1: FieldState, f2: FieldState) -> f64 {
    let diffs = [
        relative_diff(f1.pump_f, f2.pump_f),
        relative_diff(f1.pump_b, f2.pump_b),
        relative_diff(f1.sgnl_f, f2.sgnl_f),
        relative_diff(f1.sgnl_b, f2.sgnl_b),
    ];
    diffs.into_iter().sum::<f64>() / 4.0
}
pub fn profile_avg_diff(p1: &Vec<FieldState>, p2: &Vec<FieldState>) -> f64 {
    p1.iter()
        .zip(p2.iter())
        .map(|(&f1, &f2)| field_max_diff(f1, f2))
        .sum::<f64>()
        / p1.len() as f64
}

#[derive(Clone)]
pub struct FieldProfile {
    pub z: Vec<f64>,
    pub fields: Vec<FieldState>,
}

impl FieldProfile {
    pub fn new(z: Vec<f64>, fields: Vec<FieldState>) -> Self {
        assert_eq!(z.len(), fields.len());
        Self { z, fields }
    }

    pub fn sgnl_f(&self) -> impl Iterator<Item = f64> + '_ {
        // borrows from self so '_ lifetime needs to match self
        // but looks like rust can infer this so doesnt need to be explicit
        self.fields.iter().map(|x| x.sgnl_f)
    }

    pub fn sgnl_b(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.sgnl_b)
    }

    pub fn pump_f(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.pump_f)
    }

    pub fn pump_b(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.pump_b)
    }

    pub fn z(&self) -> impl Iterator<Item = f64> {
        //manual deref needed here and not for self.fields because struct field access e.g. x.sgnl_b derefs implicitly.
        //could also do self.z.iter().copied()
        self.z.iter().map(|&z| z)
    }
}

pub fn pops(fs: FieldState, fp: FibreParams) -> (f64, f64) {
    let pump_flux = fs.pump_f * fs.pump_f + fs.pump_b * fs.pump_b;
    let sgnl_flux = fs.sgnl_f * fs.sgnl_f + fs.sgnl_b * fs.sgnl_b;
    let gamma_up = pump_flux * fp.pump_ab + sgnl_flux * fp.sgnl_ab;
    let gamma_dn = pump_flux * fp.pump_em + sgnl_flux * fp.sgnl_em + 1.0 / fp.lifetime;
    let denom = gamma_up + gamma_dn;
    (gamma_dn / denom, gamma_up / denom)
}

pub fn gain(fs: FieldState, fp: FibreParams) -> (f64, f64) {
    let (g, e) = pops(fs, fp);
    (
        fp.density * (-g * fp.pump_ab + e * fp.pump_em),
        fp.density * (-g * fp.sgnl_ab + e * fp.sgnl_em),
    )
}

#[derive(Copy, Clone)]
pub struct Pump {
    pub forward: f64,
    pub backward: f64,
}

impl Default for Pump {
    fn default() -> Self {
        Self {
            forward: 10.0,
            backward: 0.0,
        }
    }
}

impl Pump {
    pub fn from_total_and_balance(total: f64, balance: f64) -> Self {
        assert!(balance >= -1.0 && balance <= 1.0 && total >= 0.0);
        let p = (balance + 1.0) * 0.5;
        Self {
            forward: (p * total).sqrt(),
            backward: ((1.0 - p) * total).sqrt(),
        }
    }
}

pub fn find_threshold_and_slope(
    pump_start: f64,
    pump_step: f64,
    ip: IterationConfig,
    mut output_power: impl FnMut(f64) -> (f64, f64, bool),
) -> Result<(f64, f64, f64), SolverError> {
    let mut current_pump = pump_start;
    let mut total_diff = -1.0;
    let mut sf = 0.0;
    let mut sb = 0.0;
    for _ in 0..ip.max {
        let (new_sf, new_sb, success) = output_power(current_pump);
        if !success {
            current_pump += pump_step;
            continue;
        }

        let new_total_diff = (new_sf + new_sb) - (sb + sf);

        if relative_diff(new_total_diff, total_diff) < ip.tol && new_total_diff > 0.0 {
            let slope_f = (new_sf - sf) / pump_step;
            let slope_b = (new_sb - sb) / pump_step;
            let threshold = current_pump - (new_sf + new_sb) / (slope_b + slope_f);
            return Ok((slope_f, slope_b, threshold));
        } else {
            current_pump += pump_step;
            total_diff = new_total_diff;
            sb = new_sb;
            sf = new_sf;
        }
    }
    Err(SolverError::ThresholdNotFound)
}

#[cfg(test)]
mod tests {
    use crate::dfb::transfer;
    use super::*;
    #[test]
    fn check_transfer() {
        let (a, b, c, d) = transfer(1.0, 0.0, 1.0);
        println!("Transfer {:?}", (a, b, c, d));
        assert_eq!(a ,(0.5_f64).exp());
    }
}