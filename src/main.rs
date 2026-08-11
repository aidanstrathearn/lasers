use lasers::lase::{
    FibreParams, FieldState, GratingProfile, GridPoints, Pump, find_lasing_profile, gain, pops,
    transfer,
};

use lasers::myplotlib::Plotter;
use lasers::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
use std::hint::black_box;
use std::time::Instant;

fn main() -> eframe::Result {
    let pu = Pump {
        forward: 10.0,
        backward: 0.0,
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
        upper: pu.forward,
        lower: 1e-8,
        midpoint: Midpoint::Geometric,
    };

    let nc = Newton1dConfig {
        tolerance: 1e-8f64,
        max_iters: 100usize,
        initial: pu.forward,
        dx: 1e-6,
    };

    let runs = 1000usize;
    let start = Instant::now();
    for _ in 0..runs {
        let result = find_lasing_profile(pu, fp, gp, kp, nc).unwrap();
        black_box(result);
    }
    let elapsed = start.elapsed();
    println!(
        "average: {:.3} µs",
        elapsed.as_secs_f64() * 1_000_000.0 / runs as f64
    );

    let result = find_lasing_profile(pu, fp, gp, kp, nc).unwrap();
    result.show()?;

    let mut plt = Plotter::new();
    let x: Vec<f64> = gp.grid(fp.length);
    let kappas: Vec<f64> = kp.grid(gp.0 + 1);
    plt.plot(&x, &kappas);
    plt.xlabel("z");
    plt.ylabel("Kappa");
    plt.title("Coupling Profile");
    plt.show()?;

    // let (p1, p2) = pops(fs, fp);
    // println!("Populations {}  {}", p1, p2);
    //
    // let (g1, g2) = gain(fs, fp);
    // println!("Gain {}  {}", g1, g2);

    let (a, b, c, d) = transfer(1.0, 0.0, 1.0);
    println!("Transfer {:?}", (a, b, c, d));
    println!("{}", a == (0.5_f64).exp());

    Ok(())
}
