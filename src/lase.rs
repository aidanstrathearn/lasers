use crate::rootfind::{newton1d, BisectionConfig, Newton1dConfig, Newton1dError, BisectionError, geometric_mid, bisection};

pub fn linspace(start: f64, stop: f64, nsteps: usize) -> Vec<f64> {
    let step: f64 = (stop - start) / (nsteps as f64);
    (0..=nsteps).map(|x| start + (x as f64) * step).collect()
}

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

#[derive(Copy, Clone)]
pub struct Reflectivities {
    left: f64,
    right: f64,
}

#[derive(Copy, Clone)]
pub struct GridPoints(usize);

#[derive(Copy, Clone)]
pub struct Discretisation {
    pub length: f64,
    pub nz: usize,
}

impl Discretisation {
    pub fn grid(self) -> Vec<f64> {
        linspace(0.0, self.length, self.nz)
    }

    pub fn step(self) -> f64 {
        self.length / self.nz as f64
    }
}

#[derive(Copy, Clone)]
pub struct GratingProfile {
    pub kappa_max: f64,
    pub pi_shift_position: f64,
}

impl GratingProfile {
    pub fn grid(self, n: usize) -> Vec<f64> {
        (0..=n)
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

#[derive(Copy, Clone, Debug)]
pub struct FieldState {
    pub sgnl_f: f64,
    pub sgnl_b: f64,
    pub pump_f: f64,
    pub pump_b: f64,
}

impl FieldState {
    fn propagate(self, fp: FibreParams, kappa: f64, dz: f64) -> Self {
        let (gp, gs) = gain(self, fp);
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

pub struct FieldProfile {
    data: Vec<FieldState>,
}

impl FieldProfile {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            data: Vec::with_capacity(n),
        }
    }

    pub fn push(&mut self, item: FieldState) {
        self.data.push(item);
    }

    pub fn sgnl_f(&self) -> impl Iterator<Item = f64> + '_ {
        // borrows from self so '_ lifetime needs to match self
        // but looks like rust can infer this so doesnt need to be explicit
        self.data.iter().map(|x| x.sgnl_f)
    }

    pub fn sgnl_b(&self) -> impl Iterator<Item = f64> {
        self.data.iter().map(|x| x.sgnl_b)
    }

    pub fn pump_f(&self) -> impl Iterator<Item = f64> {
        self.data.iter().map(|x| x.pump_f)
    }

    pub fn pump_b(&self) -> impl Iterator<Item = f64> {
        self.data.iter().map(|x| x.pump_b)
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

pub fn solve_profile(
    fs: FieldState,
    fp: FibreParams,
    ds: Discretisation,
    kappas: &[f64],
) -> FieldProfile {
    let dz = ds.step();

    let mut current = fs;
    let mut result = FieldProfile::with_capacity(ds.nz + 1);

    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
        result.push(current);
    }

    result
}

pub fn residual(fs: FieldState, fp: FibreParams, ds: Discretisation, kappas: &[f64]) -> f64 {
    let dz = ds.step();

    let mut current = fs;

    for &kappa in kappas {
        current = current.propagate(fp, kappa, dz);
    }

    current.sgnl_b
}



pub fn find_lasing(
    fs: FieldState,
    fp: FibreParams,
    ds: Discretisation,
    kp: GratingProfile,
    config: BisectionConfig,
) -> Result<FieldProfile, BisectionError> {
    let kappas = kp.grid(ds.nz);
    let trial = |sgnl_b| FieldState { sgnl_b, ..fs };
    let f = |sgnl_b| residual(trial(sgnl_b), fp, ds, &kappas);
    let sgnl_b = bisection(f, geometric_mid, config)?;
    Ok(solve_profile(trial(sgnl_b), fp, ds, &kappas))
}


pub fn find_lasing_newton(
    fs: FieldState,
    fp: FibreParams,
    ds: Discretisation,
    kp: GratingProfile,
    config: Newton1dConfig,
) -> Result<FieldProfile, Newton1dError> {
    let kappas = kp.grid(ds.nz);
    let trial = |sgnl_b| FieldState { sgnl_b, ..fs };
    let f = |sgnl_b| residual(trial(sgnl_b), fp, ds, &kappas);
    let sgnl_b = newton1d(f, config)?;
    Ok(solve_profile(trial(sgnl_b), fp, ds, &kappas))
}

pub enum FibreKind {
    Amplifier,
    DFB(GratingProfile),
    DBR(Reflectivities),
    Hybrid(GratingProfile, Reflectivities),
}

pub struct Solver {
    params: FibreParams,
    kind: FibreKind,
}
