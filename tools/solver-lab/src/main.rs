mod myplotlib;
mod plots;

use laser_solver::dfb::{dfb_pump_scan, dfb_solve, solve_profile, transfer};
use laser_solver::lase::{
    FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump, profile_max_diff,
};
use laser_solver::picard::{initial_profile, solve_profile_picard};
use laser_solver::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
use laser_solver::utils::{IterationConfig, geomspace};
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

    let ic = IterationConfig {
        max: 100usize,
        tol: 1e-10f64,
    };

    let bc = BisectionConfig {
        iteration: ic,
        upper: pu.forward,
        lower: 1e-8,
        midpoint: Midpoint::Geometric,
    };

    let nc = Newton1dConfig {
        iteration: ic,
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

    let comparison_pump = Pump {
        backward: 0.0,
        ..pu
    };
    let comparison_sgnl_b = 1.0;
    let comparison_kappas = kp.grid(gp.0);
    let comparison_boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b: comparison_sgnl_b,
        pump_f: comparison_pump.forward,
        pump_b: 0.0,
    };

    let direct_profile = FieldProfile::new(
        gp.grid(fp.length),
        solve_profile(
            comparison_boundary,
            fp,
            gp.dz(fp.length),
            &comparison_kappas,
        ),
    );
    let picard_profile = solve_profile_picard(
        comparison_sgnl_b,
        comparison_pump,
        fp,
        gp,
        ic,
        &comparison_kappas,
    )
    .expect("Picard profile comparison did not converge");

    let max_diff = profile_max_diff(&direct_profile.fields, &picard_profile.fields);
    println!("Picard/direct profile max diff: {max_diff:e}");

    let z: Vec<f64> = direct_profile.z().collect();
    let direct_sgnl_f: Vec<f64> = direct_profile.sgnl_f().collect();
    let picard_sgnl_f: Vec<f64> = picard_profile.sgnl_f().collect();
    let direct_sgnl_b: Vec<f64> = direct_profile.sgnl_b().collect();
    let picard_sgnl_b: Vec<f64> = picard_profile.sgnl_b().collect();
    let diff_f: Vec<f64> = direct_sgnl_f
        .iter()
        .zip(picard_sgnl_f.iter())
        .map(|(a, b)| a - b)
        .collect();

    let diff_b: Vec<f64> = direct_sgnl_b
        .iter()
        .zip(picard_sgnl_b.iter())
        .map(|(a, b)| a - b)
        .collect();

    let mut plt = Plotter::new();
    plt.plot(&z, &direct_sgnl_f).label("Direct forward signal");
    plt.plot(&z, &picard_sgnl_f).label("Picard forward signal");
    plt.plot(&z, &direct_sgnl_b).label("Direct backward signal");
    plt.plot(&z, &picard_sgnl_b).label("Picard backward signal");
    plt.plot(&z, &diff_b).label("backward diff");
    plt.plot(&z, &diff_f).label("forward diff");
    plt.xlabel("z");
    plt.ylabel("Field amplitude");
    plt.title("Direct vs Picard profile (zero backward pump)");
    plt.show()?;

    Ok(())
}
