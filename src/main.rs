use lasers::lase::{
    FibreParams, FieldState, GratingProfile, GridPoints, find_lasing, gain, pops, transfer,
};

use lasers::myplotlib::Plotter;
use lasers::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
use std::hint::black_box;
use std::time::Instant;

fn main() -> eframe::Result {
    let fs = FieldState {
        sgnl_f: 0.0,
        sgnl_b: 0.0001,
        pump_f: 10.0,
        pump_b: 0.0,
    };

    let fp = FibreParams {
        density: 1.0,
        lifetime: 1.0,
        pump_ab: 0.01,
        pump_em: 0.0,
        sgnl_ab: 0.0,
        sgnl_em: 1.0,
        length: 10.0,
    };

    let gp = GridPoints(500);

    let kp = GratingProfile {
        kappa_max: 1.0,
        pi_shift_position: 0.5,
    };

    let bc = BisectionConfig {
        tolerance: 1e-8f64,
        max_iters: 100usize,
        upper: fs.pump_f,
        lower: fs.sgnl_b,
        midpoint: Midpoint::Geometric,
    };

    let nc = Newton1dConfig {
        tolerance: 1e-8f64,
        max_iters: 100usize,
        initial: fs.pump_f,
        dx: 1e-6,
    };

    let kappa: Vec<f64> = kp.grid(gp.0);

    let start = Instant::now();
    let result = find_lasing(fs, fp, gp, kp, nc).unwrap();
    let elapsed = start.elapsed();
    println!("{:?}", elapsed);

    let runs = 1000usize;
    let start = Instant::now();
    for _ in 0..runs {
        let result = find_lasing(fs, fp, gp, kp, nc).unwrap();
        black_box(result);
    }
    let elapsed = start.elapsed();
    println!(
        "average: {:.3} µs",
        elapsed.as_secs_f64() * 1_000_000.0 / runs as f64
    );

    let x = gp.grid(fp.length);
    let pump_f: Vec<f64> = result.pump_f().collect();
    let pump_b: Vec<f64> = result.pump_b().collect();
    let sgnl_f: Vec<f64> = result
        .sgnl_f()
        .map(|x| (x.powi(2) + 1e-6).log10())
        .collect();
    let sgnl_b: Vec<f64> = result
        .sgnl_b()
        .map(|x| (x.powi(2) + 1e-6).log10())
        .collect();

    let mut plt = Plotter::new();
    plt.plot(&x, &kappa).label("Kappa");
    plt.plot(&x, &pump_f).label("Forward Pump");
    plt.plot(&x, &pump_b).label("Backward Pump");
    plt.plot(&x, &sgnl_f).label("Forward signal");
    plt.plot(&x, &sgnl_b).label("Backward signal");

    plt.xlabel("z");
    plt.ylabel("Amplitude");
    plt.title("Fields");
    plt.show()?;

    let (p1, p2) = pops(fs, fp);
    println!("Populations {}  {}", p1, p2);

    let (g1, g2) = gain(fs, fp);
    println!("Gain {}  {}", g1, g2);

    let (a, b, c, d) = transfer(1.0, 0.0, 1.0);
    println!("Transfer {:?}", (a, b, c, d));
    println!("{}", a == (0.5_f64).exp());

    Ok(())
}
