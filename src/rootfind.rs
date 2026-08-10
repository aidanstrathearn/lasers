const DEFAULT_TOLERANCE: f64 = 1e-8;
const DEFAULT_MAX_ITERS: usize = 100;

#[derive(Debug)]
pub enum RootFindError {
    RootNotBracketed,
    DidNotConverge,
}

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
    pub tolerance: f64,
    pub max_iters: usize,
    pub initial: f64,
    pub dx: f64,
}

impl Default for Newton1dConfig{
    fn default() -> Self {
        Self {
            tolerance: DEFAULT_TOLERANCE,
            max_iters: DEFAULT_MAX_ITERS,
            initial: 1.0,
            dx: 1e-8,
        }
    }
}

pub fn newton1d(f: impl Fn(f64) -> f64, config: Newton1dConfig) -> Result<f64, RootFindError> {
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
    Err(RootFindError::DidNotConverge)
}

#[derive(Copy, Clone)]
pub struct BisectionConfig {
    pub tolerance: f64,
    pub max_iters: usize,
    pub upper: f64,
    pub lower: f64,
    pub midpoint: Midpoint,
}

impl Default for BisectionConfig {
    fn default() -> Self {
        Self {
            tolerance: DEFAULT_TOLERANCE,
            max_iters: DEFAULT_MAX_ITERS,
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

    Err(RootFindError::DidNotConverge)
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