//! One-dimensional root-finding algorithms used by the laser solver.

use crate::utils::IterationConfig;
use std::fmt;

#[derive(Debug)]
pub enum RootFindError {
    RootNotBracketed,
    DidNotConverge,
}

impl fmt::Display for RootFindError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RootNotBracketed => {
                write!(formatter, "root is not bracketed")
            }
            Self::DidNotConverge => {
                write!(formatter, "root finder did not converge")
            }
        }
    }
}

impl std::error::Error for RootFindError {}

// pub struct RfConfig<G> {
//     iteration: IterationConfig,
//     method: G
// }
//
// pub struct NConfig {
//     initial: f64,
//     dx: f64
// }
//
// impl Default for NConfig {
//     fn default() -> Self {
//         Self { initial: 0.0, dx: 1e-6}
//     }
// }
//
// pub struct BConfig {
//     upper: f64,
//     lower: f64,
//     midpoint: Midpoint
// }
//
// impl Default for BConfig {
//     fn default() -> Self {
//         Self { upper: 2.0, lower: 1e-6, midpoint: Midpoint::Geometric}
//     }
// }
//
// impl<G> RfConfig<G> {
//     pub fn max_iters(&mut self, max: usize) -> &mut Self {
//         self.iteration.max = max;
//         self
//     }
//     pub fn tolerance(&mut self, tol: f64) -> &mut Self {
//         self.iteration.tol = tol;
//         self
//     }
// }
//
// impl RfConfig<NConfig> {
//     pub fn newton() -> Self {
//         RfConfig {iteration: IterationConfig::default(), method: NConfig::default()}
//     }
//     pub fn initial(&mut self, initial: f64) -> &mut Self {
//         self.method.initial = initial;
//         self
//     }
//     pub fn dx(&mut self, dx: f64) -> &mut Self {
//         self.method.dx = dx;
//         self
//     }
// }
//
// impl RfConfig<BConfig> {
//     pub fn bisection() -> Self {
//         RfConfig {iteration: IterationConfig::default(), method: BConfig::default()}
//     }
//     pub fn upper(&mut self, upper: f64) -> &mut Self {
//         self.method.upper = upper;
//         self
//     }
//     pub fn lower(&mut self, lower: f64) -> &mut Self {
//         self.method.lower = lower;
//         self
//     }
//     pub fn midpoint(&mut self, midpoint: Midpoint) -> &mut Self {
//         self.method.midpoint = midpoint;
//         self
//     }
// }

#[derive(Copy, Clone)]
pub enum RootFindConfig {
    Newton1d(Newton1dConfig),
    Bisection(BisectionConfig),
}

impl From<Newton1dConfig> for RootFindConfig {
    fn from(config: Newton1dConfig) -> Self {
        Self::Newton1d(config)
    }
}

impl From<BisectionConfig> for RootFindConfig {
    fn from(config: BisectionConfig) -> Self {
        Self::Bisection(config)
    }
}

#[derive(Copy, Clone)]
pub struct Newton1dConfig {
    pub iteration: IterationConfig,
    pub initial: f64,
    pub dx: f64,
}

impl Default for Newton1dConfig {
    fn default() -> Self {
        Self {
            iteration: IterationConfig::default(),
            initial: 1.0,
            dx: 1e-8,
        }
    }
}

pub fn newton1d(f: impl Fn(f64) -> f64, config: Newton1dConfig) -> Result<f64, RootFindError> {
    let dx = config.dx;
    let mut x = config.initial;
    for _ in 0..config.iteration.max {
        let fx = f(x);
        if fx.abs() < config.iteration.tol {
            return Ok(x);
        }
        let dfdx = (f(x + dx) - fx) / dx;
        x -= fx / dfdx;
    }
    Err(RootFindError::DidNotConverge)
}

pub fn try_newton1d<F, E>(f: F, config: Newton1dConfig) -> Result<f64, E>
where
    F: Fn(f64) -> Result<f64, E>,
    E: From<RootFindError>,
{
    let mut x = config.initial;
    for _ in 0..config.iteration.max {
        let fx = f(x)?;
        if fx.abs() < config.iteration.tol {
            return Ok(x);
        }
        let dfdx = (f(x + config.dx)? - fx) / config.dx;
        x -= fx / dfdx;
    }

    Err(RootFindError::DidNotConverge.into())
}

#[derive(Copy, Clone)]
pub struct BisectionConfig {
    pub iteration: IterationConfig,
    pub upper: f64,
    pub lower: f64,
    pub midpoint: Midpoint,
}

impl Default for BisectionConfig {
    fn default() -> Self {
        Self {
            iteration: IterationConfig::default(),
            upper: 1.0,
            lower: 1e-8,
            midpoint: Midpoint::Geometric,
        }
    }
}

#[derive(Copy, Clone)]
pub enum Midpoint {
    Arithmetic,
    Geometric,
}

pub fn bisection(
    f: impl Fn(f64) -> f64,
    mid: impl Fn(f64, f64) -> f64,
    config: BisectionConfig,
) -> Result<f64, RootFindError> {
    let mut lower = config.lower;
    let mut upper = config.upper;
    let mut f_lower = f(lower);
    let f_upper = f(upper);

    if f_lower * f_upper > 0.0 {
        return Err(RootFindError::RootNotBracketed);
    }

    for _ in 0..config.iteration.max {
        let midpoint = mid(upper, lower);
        let f_midpoint = f(midpoint);

        if f_midpoint.abs() < config.iteration.tol {
            return Ok(midpoint);
        }

        if f_lower * f_midpoint <= 0.0 {
            upper = midpoint;
        } else {
            lower = midpoint;
            f_lower = f_midpoint;
        }
    }

    Err(RootFindError::DidNotConverge)
}

pub fn try_bisection<F, E, M>(f: F, mid: M, config: BisectionConfig) -> Result<f64, E>
where
    F: Fn(f64) -> Result<f64, E>,
    E: From<RootFindError>,
    M: Fn(f64, f64) -> f64,
{
    let mut lower = config.lower;
    let mut upper = config.upper;
    let mut f_lower = f(lower)?;
    let f_upper = f(upper)?;

    if f_lower * f_upper > 0.0 {
        return Err(RootFindError::RootNotBracketed.into());
    }

    for _ in 0..config.iteration.max {
        let midpoint = mid(upper, lower);
        let f_midpoint = f(midpoint)?;

        if f_midpoint.abs() < config.iteration.tol {
            return Ok(midpoint);
        }

        if f_lower * f_midpoint <= 0.0 {
            upper = midpoint;
        } else {
            lower = midpoint;
            f_lower = f_midpoint;
        }
    }

    Err(RootFindError::DidNotConverge.into())
}

pub fn arithmetic_mid(a: f64, b: f64) -> f64 {
    0.5 * (a + b)
}

pub fn geometric_mid(a: f64, b: f64) -> f64 {
    (a * b).sqrt()
}

pub fn bracket_bisection(
    f: impl Fn(f64) -> f64,
    config: BisectionConfig,
) -> Result<f64, RootFindError> {
    match config.midpoint {
        Midpoint::Arithmetic => bisection(f, arithmetic_mid, config),
        Midpoint::Geometric => bisection(f, geometric_mid, config),
    }
}

pub fn try_bracket_bisection<F, E>(f: F, config: BisectionConfig) -> Result<f64, E>
where
    F: Fn(f64) -> Result<f64, E>,
    E: From<RootFindError>,
{
    match config.midpoint {
        Midpoint::Arithmetic => try_bisection(f, arithmetic_mid, config),
        Midpoint::Geometric => try_bisection(f, geometric_mid, config),
    }
}

pub fn rootfind_1d(
    f: impl Fn(f64) -> f64,
    config: impl Into<RootFindConfig>,
) -> Result<f64, RootFindError> {
    match config.into() {
        // seems a bit silly to into() and then immediately destructure
        RootFindConfig::Newton1d(n_config) => newton1d(f, n_config),
        RootFindConfig::Bisection(b_config) => bracket_bisection(f, b_config),
    }
}

pub fn try_rootfind_1d<F, E>(f: F, config: impl Into<RootFindConfig>) -> Result<f64, E>
where
    F: Fn(f64) -> Result<f64, E>,
    E: From<RootFindError>,
{
    match config.into() {
        RootFindConfig::Newton1d(config) => try_newton1d(f, config),
        RootFindConfig::Bisection(config) => try_bracket_bisection(f, config),
    }
}
