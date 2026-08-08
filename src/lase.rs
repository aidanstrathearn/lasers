pub fn linspace(start: f64, stop: f64, nsteps: usize) -> Vec<f64> {
    let step: f64 = (stop - start) / (nsteps as f64);
    (0..=nsteps).map(|x| start + (x as f64) * step).collect()
}

#[derive(Copy, Clone)]
pub struct Dopant {
    pub density: f64,
    pub lifetime: f64,
    pub pump_ab: f64,
    pub pump_em: f64,
    pub sgnl_ab: f64,
    pub sgnl_em: f64,
    pub length: f64
}

#[derive(Copy, Clone)]
pub struct Reflectivities {
left: f64,
right: f64
}

//#[derive(Copy, Clone)]
//pub struct Fibre {
//	pub length: f64,
//	pub left_reflection: f64,
//	pub right_reflection: f64,
//	}

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
        (0..=n).map(|j| { let z = j as f64 / n as f64;
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
    fn propagate(self, dp: Dopant, kappa: f64, dz: f64) -> Self {
        let (gp, gs) = gain(self, dp);
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

pub fn pops(fs: FieldState, dp: Dopant) -> (f64, f64) {
    let pump_flux = fs.pump_f * fs.pump_f + fs.pump_b * fs.pump_b;
    let sgnl_flux = fs.sgnl_f * fs.sgnl_f + fs.sgnl_b * fs.sgnl_b;
    let gamma_up = pump_flux * dp.pump_ab + sgnl_flux * dp.sgnl_ab;
    let gamma_dn = pump_flux * dp.pump_em + sgnl_flux * dp.sgnl_em + 1.0 / dp.lifetime;
    let denom = gamma_up + gamma_dn;
    (gamma_dn / denom, gamma_up / denom)
}

pub fn gain(fs: FieldState, dp: Dopant) -> (f64, f64) {
    let (g, e) = pops(fs, dp);
    (
        dp.density * (-g * dp.pump_ab + e * dp.pump_em),
        dp.density * (-g * dp.sgnl_ab + e * dp.sgnl_em),
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
    dp: Dopant,
    ds: Discretisation,
    kappas: &[f64],
) -> FieldProfile {
    let dz = ds.step();

    let mut current = fs;
    let mut result = FieldProfile::with_capacity(ds.nz + 1);

    for &kappa in kappas {
        current = current.propagate(dp, kappa, dz);
        result.push(current);
    }

    result
}

pub fn residual(fs: FieldState, dp: Dopant, ds: Discretisation, kappas: &[f64]) -> f64 {
    let dz = ds.step();

    let mut current = fs;

    for &kappa in kappas {
        current = current.propagate(dp, kappa, dz);
    }

    current.sgnl_b
}

#[derive(Copy, Clone)]
pub struct BisectionConfig {
    pub tolerance: f64,
    pub max_iters: usize,
    pub upper: f64,
    pub lower: f64,
}

#[derive(Debug)]
pub enum BisectionError {
    DidNotConverge,
    RootNotBracketed,
}

pub fn arithmetic_mid(a: f64, b: f64) -> f64 {
    0.5 * (a + b)
}

pub fn geometric_mid(a: f64, b: f64) -> f64 {
    (a * b).sqrt()
}

fn bisection(
    f: impl Fn(f64) -> f64,
    mid: impl Fn(f64, f64) -> f64,
    config: BisectionConfig,
) -> Result<f64, BisectionError> {
    let mut lower = config.lower;
    let mut upper = config.upper;
    let mut f_lower = f(lower);
    let f_upper = f(upper);

    if f_lower * f_upper > 0.0 {
        return Err(BisectionError::RootNotBracketed);
    }

    for _ in 0..config.max_iters {
        let midpoint = mid(upper, lower);
        let f_midpoint = f(midpoint);

        if f_midpoint.abs() < config.tolerance {
            return Ok(midpoint);
        }

        if f_lower * f_midpoint <= 0.0 {
            upper = midpoint;
        } else {
            lower = midpoint;
            f_lower = f_midpoint;
        }
    }

    Err(BisectionError::DidNotConverge)
}

pub fn find_lasing(
    fs: FieldState,
    dp: Dopant,
    ds: Discretisation,
    kp: GratingProfile,
    config: BisectionConfig,
) -> Result<FieldProfile, BisectionError> {
    let kappas = kp.grid(ds.nz);
    let trial = |sgnl_b| FieldState { sgnl_b, ..fs };
    let f = |sgnl_b| residual(trial(sgnl_b), dp, ds, &kappas);
    let sgnl_b = bisection(f, geometric_mid, config)?;
    Ok(solve_profile(trial(sgnl_b), dp, ds, &kappas))
}

#[derive(Debug)]
pub enum Newton1dError {
    DidNotConverge,
}

#[derive(Copy, Clone)]
pub struct Newton1dConfig {
    pub tolerance: f64,
    pub max_iters: usize,
    pub initial: f64,
    pub dx: f64,
}

fn newton1d(f: impl Fn(f64) -> f64, config: Newton1dConfig) -> Result<f64, Newton1dError> {
    let dx = config.dx;
    let mut x = config.initial;
    for _ in 0..config.max_iters {
        let fx = f(x);
        if fx.abs() < config.tolerance {
            return Ok(x);
        }
        let dfdx = (f(x + dx) - fx) / dx;
        x -= fx / dfdx;
    }
    Err(Newton1dError::DidNotConverge)
}

pub fn find_lasing_newton(
    fs: FieldState,
    dp: Dopant,
    ds: Discretisation,
    kp: GratingProfile,
    config: Newton1dConfig,
) -> Result<FieldProfile, Newton1dError> {
    let kappas = kp.grid(ds.nz);
    let trial = |sgnl_b| FieldState { sgnl_b, ..fs };
    let f = |sgnl_b| residual(trial(sgnl_b), dp, ds, &kappas);
    let sgnl_b = newton1d(f, config)?;
    Ok(solve_profile(trial(sgnl_b), dp, ds, &kappas))
}

pub enum FibreKind {
	Amplifier,
	DFB(GratingProfile),
	DBR(Reflectivities),
	Hybrid(GratingProfile, Reflectivities)
	}

pub struct Solver {
params: Dopant,
kind: FibreKind
}


