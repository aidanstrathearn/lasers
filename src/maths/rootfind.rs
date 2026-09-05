//! One-dimensional root-finding algorithms used by the laser solver.

use super::utils::IterationConfig;
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

fn try_newton1d<F, E>(mut f: F, config: Newton1dConfig) -> Result<f64, E>
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
    Interp,
}

#[derive(Copy, Clone)]
struct Bracket {
    lower: f64,
    upper: f64,
    f_lower: f64,
    f_upper: f64,
}

fn try_bisection<F, E, M>(mut f: F, mid: M, config: BisectionConfig) -> Result<f64, E>
where
    F: FnMut(f64) -> Result<f64, E>,
    E: From<RootFindError>,
    M: Fn(Bracket) -> f64,
{
    let mut bracket = Bracket {
        lower: config.lower,
        upper: config.upper,
        f_lower: f(config.lower)?,
        f_upper: f(config.upper)?,
    };

    if bracket.f_lower * bracket.f_upper > 0.0 {
        return Err(RootFindError::RootNotBracketed.into());
    }

    for _ in 0..config.iteration.max {
        let midpoint = mid(bracket);
        let f_midpoint = f(midpoint)?;

        if f_midpoint.abs() < config.iteration.tol {
            return Ok(midpoint);
        }

        if bracket.f_lower * f_midpoint <= 0.0 {
            bracket.upper = midpoint;
            bracket.f_upper = f_midpoint;
        } else {
            bracket.lower = midpoint;
            bracket.f_lower = f_midpoint;
        }
    }

    Err(RootFindError::DidNotConverge.into())
}

fn arithmetic_mid(bracket: Bracket) -> f64 {
    0.5 * (bracket.lower + bracket.upper)
}

fn geometric_mid(bracket: Bracket) -> f64 {
    (bracket.lower * bracket.upper).sqrt()
}

fn interp_mid(bracket: Bracket) -> f64 {
    ((bracket.f_upper * bracket.lower.ln() - bracket.f_lower * bracket.upper.ln())
        / (bracket.f_upper - bracket.f_lower))
        .exp()
}

fn try_bracket_bisection<F, E>(f: F, config: BisectionConfig) -> Result<f64, E>
where
    F: FnMut(f64) -> Result<f64, E>,
    E: From<RootFindError>,
{
    let mid = match config.midpoint {
        Midpoint::Arithmetic => arithmetic_mid,
        Midpoint::Geometric => geometric_mid,
        Midpoint::Interp => interp_mid,
    };
    try_bisection(f, mid, config)
}

pub fn rootfind_1d<F>(mut f: F, config: impl Into<RootFindConfig>) -> Result<f64, RootFindError>
where
    F: FnMut(f64) -> f64,
{
    try_rootfind_1d(|x| Ok::<_, RootFindError>(f(x)), config)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn bisection_config(midpoint: Midpoint) -> BisectionConfig {
        BisectionConfig {
            iteration: IterationConfig {
                max: 100,
                tol: 1e-10,
            },
            lower: 1.0,
            upper: 4.0,
            midpoint,
        }
    }

    #[test]
    fn infallible_wrapper_supports_all_bisection_midpoints() {
        for midpoint in [Midpoint::Arithmetic, Midpoint::Geometric, Midpoint::Interp] {
            let root = rootfind_1d(|x| x - 2.0, bisection_config(midpoint)).unwrap();
            assert!((root - 2.0).abs() < 1e-8);
        }
    }

    #[test]
    fn infallible_wrapper_supports_newton_and_fn_mut() {
        let mut calls = 0;
        let root = rootfind_1d(
            |x| {
                calls += 1;
                x * x - 2.0
            },
            Newton1dConfig::default(),
        )
        .unwrap();

        assert!((root - 2.0_f64.sqrt()).abs() < 1e-8);
        assert!(calls > 0);
    }

    #[test]
    fn fallible_rootfinder_supports_all_bisection_midpoints() {
        for midpoint in [Midpoint::Arithmetic, Midpoint::Geometric, Midpoint::Interp] {
            let root = try_rootfind_1d(
                |x| Ok::<_, RootFindError>(x - 2.0),
                bisection_config(midpoint),
            )
            .unwrap();
            assert!((root - 2.0).abs() < 1e-8);
        }
    }

    #[derive(Debug, PartialEq)]
    enum TestError {
        Callback,
        RootFind,
    }

    impl From<RootFindError> for TestError {
        fn from(_: RootFindError) -> Self {
            Self::RootFind
        }
    }

    #[test]
    fn fallible_rootfinder_propagates_callback_errors() {
        let result = try_rootfind_1d(
            |_| Err::<f64, _>(TestError::Callback),
            bisection_config(Midpoint::Arithmetic),
        );

        assert_eq!(result, Err(TestError::Callback));
    }
}
