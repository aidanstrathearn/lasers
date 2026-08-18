use crate::rootfind::{RootFindConfig, RootFindError, rootfind_1d};
use crate::utils::{linspace, relative_diff};

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
            length: 10.0,
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
        Self(500)
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
    pub kappa_max: f64,
    pub pi_shift_position: f64,
}

impl Default for GratingProfile {
    fn default() -> Self {
        Self {
            kappa_max: 1.0,
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
                    self.kappa_max
                } else {
                    -self.kappa_max
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
    pub fn propagate(self, fp: FibreParams, kappa: f64, dz: f64) -> Self {
        self.general_step(self, fp, kappa, dz)
    }

    pub fn general_step(self, other: Self, fp: FibreParams, kappa: f64, dz: f64) -> Self {
        let (gp, gs) = gain(other, fp);
        let (a, b, c, d) = transfer(gs, kappa, dz);
        let expg = (0.5 * gp * dz).exp();

        FieldState {
            sgnl_f: a * self.sgnl_f + b * self.sgnl_b,
            sgnl_b: c * self.sgnl_f + d * self.sgnl_b,
            pump_f: self.pump_f * expg,
            pump_b: self.pump_b / expg,
        }
    }
}

pub fn field_max_diff(f1: FieldState, f2: FieldState) -> f64 {
    let diffs = [
        relative_diff(f1.pump_f, f2.pump_f),
        relative_diff(f1.pump_b, f2.pump_b),
        relative_diff(f1.sgnl_f, f2.sgnl_f),
        relative_diff(f1.sgnl_b, f2.sgnl_b),
    ];
    diffs.into_iter().reduce(f64::max).unwrap_or(f64::NAN)
}
pub fn profile_max_diff(p1: &Vec<FieldState>, p2: &Vec<FieldState>) -> f64 {
    p1.iter()
        .zip(p2.iter())
        .map(|(&f1, &f2)| field_max_diff(f1, f2))
        .reduce(f64::max)
        .unwrap_or(f64::NAN)
}

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

    pub fn plotpoints(&self, field: &str) -> Vec<[f64; 2]> {
        let z = self.z();
        match field {
            "sgnl_b" => {
                let field = self.sgnl_b();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }

            "sgnl_f" => {
                let field = self.sgnl_f();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }

            "pump_f" => {
                let field = self.pump_f();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }

            "pump_b" => {
                let field = self.pump_b();
                z.zip(field).map(|(x, y)| [x, y]).collect()
            }
            _ => panic!(),
        }
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

pub fn solve_profile(fs: FieldState, fp: FibreParams, dz: f64, kappas: &[f64]) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(kappas.len() + 1);
    result.push(current);
    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
        result.push(current);
    }
    result
}

pub fn out_field(fs: FieldState, fp: FibreParams, dz: f64, kappas: &[f64]) -> FieldState {
    let mut current = fs;
    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
    }
    current
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

pub fn dfb_solve(
    pu: Pump,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    full_profile: bool,
    config: impl Into<RootFindConfig>,
) -> Result<FieldProfile, RootFindError> {
    let kappas = kp.grid(gp.0);
    let dz = gp.dz(fp.length);
    let trial = |sgnl_b| FieldState {
        sgnl_f: 0.0,
        sgnl_b: sgnl_b,
        pump_f: pu.forward,
        pump_b: 0.0, // bidirectional pump not implemented yet
    };
    let f = |sgnl_b| out_field(trial(sgnl_b), fp, dz, &kappas).sgnl_b;
    let sgnl_b = rootfind_1d(f, config)?;

    if full_profile {
        let z = gp.grid(fp.length);
        let fields = solve_profile(trial(sgnl_b), fp, dz, &kappas);
        Ok(FieldProfile::new(z, fields))
    } else {
        let z = vec![0.0_f64, fp.length];
        let out_left = trial(sgnl_b);
        let fields = vec![out_left, out_field(out_left, fp, dz, &kappas)];
        Ok(FieldProfile::new(z, fields))
    }
}

pub fn dfb_threshold_curve(
    pumps: Vec<f64>,
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Vec<Result<f64, RootFindError>> {
    pumps
        .iter()
        .map(|&pmp| {
            let pu = Pump {
                forward: pmp,
                backward: 0.0,
            };
            let result = dfb_solve(pu, fp, gp, kp, false, config)?;
            Ok(result.sgnl_f().last().unwrap())
        })
        .collect()
}

pub fn dfb_threshold_curve_with_zeros(
    pumps: &[f64],
    fp: FibreParams,
    gp: GridPoints,
    kp: GratingProfile,
    config: impl Into<RootFindConfig> + Copy,
) -> Vec<(f64, f64)> {
    let full_profile = false;
    pumps
        .iter()
        .map(|&pmp| {
            let pu = Pump {
                forward: pmp,
                backward: 0.0,
            };
            let result = dfb_solve(pu, fp, gp, kp, full_profile, config)?;
            Ok((
                result.sgnl_f().last().unwrap(),
                result.sgnl_b().next().unwrap(),
            ))
        })
        .map(|result: Result<(f64, f64), RootFindError>| result.unwrap_or((0.0, 0.0)))
        .collect()

    // pumps
    //     .iter()
    //     .map(|&pmp| {
    //         let pu = Pump {
    //             forward: pmp,
    //             backward: 0.0,
    //         };
    //
    //         dfb_solve(pu, fp, gp, kp, full_profile, config).map_or((0.0, 0.0), |result| {
    //             (
    //                 result.sgnl_f().last().unwrap(),
    //                 result.sgnl_b().next().unwrap(),
    //             )
    //         })
    //     })
    //     .collect()
}
