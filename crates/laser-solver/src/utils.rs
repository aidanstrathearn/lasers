pub fn linspace(start: f64, stop: f64, nsteps: usize) -> Vec<f64> {
    let step: f64 = (stop - start) / (nsteps as f64);
    (0..=nsteps).map(|x| start + (x as f64) * step).collect()
}

pub fn geomspace(start: f64, stop: f64, nsteps: usize) -> Vec<f64> {
    linspace(start, stop, nsteps)
        .iter()
        .map(|&x| 10.0_f64.powf(x))
        .collect()
}

pub fn relative_diff(x1: f64, x2: f64) -> f64 {
    let scale = x1.abs().max(x2.abs());
    if scale == 0.0 {
        0.0
    } else {
        (x1 - x2).abs() / scale
    }
}
