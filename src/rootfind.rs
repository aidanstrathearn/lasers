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

#[derive(Copy, Clone)]
pub struct Newton1dConfig {
    pub tolerance: f64,
    pub max_iters: usize,
    pub initial: f64,
    pub dx: f64,
}

impl Into<RootFindConfig> for Newton1dConfig {
    fn into(self) -> RootFindConfig {
        RootFindConfig::Newton1d(self)
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
    pub midpoint: Midpoint
}

#[derive(Copy, Clone)]
pub enum Midpoint {
    Arithmetic,
    Geometric
}

impl Into<RootFindConfig> for BisectionConfig {
    fn into(self) -> RootFindConfig {
        RootFindConfig::Bisection(self)
    }
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
        Midpoint::Geometric => bisection(f, geometric_mid, config)
    }
    
}

pub fn geometric_bisection(
    f: impl Fn(f64) -> f64,
    config: BisectionConfig,
) -> Result<f64, RootFindError> {
    bisection(f, arithmetic_mid, config)
}

pub fn arithmetic_bisection(
    f: impl Fn(f64) -> f64,
    config: BisectionConfig,
) -> Result<f64, RootFindError> {
    bisection(f, arithmetic_mid, config)
}
