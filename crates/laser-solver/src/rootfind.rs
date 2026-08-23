//! One-dimensional root-finding algorithms used by the laser solver.

use crate::utils::IterationConfig;
use std::fmt;

#[derive(Debug)]
pub enum RootFindError {
    RootNotBracketed,
    DidNotConverge,
    NotImplemented
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

            Self::NotImplemented => {
                write!(formatter, "method not implemented")
            }
        }
    }
}

impl std::error::Error for RootFindError {}


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

pub fn try_newton1d<F, E>(mut f: F, config: Newton1dConfig) -> Result<f64, E>
where
    F: FnMut(f64) -> Result<f64, E>,
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
    Interp
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

pub fn bisection_lininterp(
    f: impl Fn(f64) -> f64,
    config: BisectionConfig,
) -> Result<f64, RootFindError> {
    let mut lower = config.lower;
    let mut upper = config.upper;
    let mut f_lower = f(lower);
    let mut f_upper = f(upper);

    if f_lower * f_upper > 0.0 {
        return Err(RootFindError::RootNotBracketed);
    }

    for _ in 0..config.iteration.max {
        //let midpoint = (f_upper * lower - f_lower * upper) / (f_upper - f_lower) ;
        let midpoint = ((f_upper * lower.ln() - f_lower * upper.ln())
            / (f_upper - f_lower))
            .exp();
        let f_midpoint = f(midpoint);

        if f_midpoint.abs() < config.iteration.tol {
            return Ok(midpoint);
        }

        if f_lower * f_midpoint <= 0.0 {
            upper = midpoint;
            f_upper = f_midpoint
        } else {
            lower = midpoint;
            f_lower = f_midpoint;
        }
    }

    Err(RootFindError::DidNotConverge)
}

pub fn try_bisection<F, E, M>(mut f: F, mid: M, config: BisectionConfig) -> Result<f64, E>
// foo(mut f: F) {} roughly similar to foo(f: F) { let mut f = f ... }
where
    F: FnMut(f64) -> Result<f64, E>,
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
        Midpoint::Interp => bisection_lininterp(f, config)
    }
}

pub fn try_bracket_bisection<F, E>(f: F, config: BisectionConfig) -> Result<f64, E>
where
    F: FnMut(f64) -> Result<f64, E>,
    E: From<RootFindError>,
{
    match config.midpoint {
        Midpoint::Arithmetic => try_bisection(f, arithmetic_mid, config),
        Midpoint::Geometric => try_bisection(f, geometric_mid, config),
        Midpoint::Interp => Err(RootFindError::NotImplemented.into())
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
    F: FnMut(f64) -> Result<f64, E>,
    E: From<RootFindError>,
{
    match config.into() {
        RootFindConfig::Newton1d(config) => try_newton1d(f, config),
        RootFindConfig::Bisection(config) => try_bracket_bisection(f, config),
    }
}
