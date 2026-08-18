mod myplotlib;
mod plots;

use laser_solver::dfb::{dfb_pump_scan, dfb_solve, transfer};
use laser_solver::lase::{FibreParams, GratingProfile, GridPoints, Pump};
use laser_solver::picard::{find_pump_b, initial_profile};
use laser_solver::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
use laser_solver::utils::geomspace;
use myplotlib::Plotter;
use plots::show_field_profile;
use std::hint::black_box;
use std::time::Instant;

fn main() -> eframe::Result {
    let pu = Pump {
        forward: 100.0,
        backward: 10.0,
    };

    let fp = FibreParams {
        density: 1.0,
        lifetime: 1.0,
        pump_ab: 0.01 * 100.0,
        pump_em: 0.0,
        sgnl_ab: 0.0,
        sgnl_em: 1.0,
        length: 10.0,
    };

    let gp = GridPoints(500);
    let full_profile = true;

    let kp = GratingProfile {
        kappa_left: 1.0,
        kappa_right: 1.0,
        pi_shift_position: 0.45,
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
        let result = dfb_solve(pu, fp, gp, kp, full_profile, nc).unwrap();
        black_box(result);
    }
    let elapsed = start.elapsed();
    println!(
        "average: {:.3} µs",
        elapsed.as_secs_f64() * 1_000_000.0 / runs as f64
    );

    let result = dfb_solve(pu, fp, gp, kp, full_profile, nc).unwrap();
    show_field_profile(&result)?;

    let profile = initial_profile(pu, fp, gp);
    //let calc = profile.fields[0].pump_b;
    show_field_profile(&profile)?;
    // let pump_b_left = find_pump_b_out(pu.backward, profile, fp, gp.dz(fp.length));
    // println!("caclulate {}   actual {}", pump_b_left, calc);

    let mut plt = Plotter::new();
    let mut pumps = geomspace(-1.0, 1.0, 200);
    let start = Instant::now();
    let threshold = dfb_pump_scan(&pumps, fp, gp, kp, bc);
    let sgnl_f: Vec<f64> = threshold.iter().map(|x| x.0).collect();
    let sgnl_b: Vec<f64> = threshold.iter().map(|x| x.1).collect();
    let elapsed = start.elapsed();
    println!("pump sweep {:.3}", elapsed.as_secs_f64());

    pumps = pumps.iter().map(|x| x.log10()).collect();
    //threshold = threshold.iter().map(|x| x.max(1e-10).log10()).collect();
    plt.plot(&pumps, &sgnl_f).label("Forward");
    plt.plot(&pumps, &sgnl_b).label("Back");
    plt.show()?;

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
