mod myplotlib;
mod plots;

use crate::plots::plot_profile_diff;
use laser_solver::dfb::{
    dfb_find_threshold_and_slope, dfb_pump_scan, dfb_solve_shooting, out_field, solve_profile, GratingProfile,
};
use laser_solver::lase::{
    FibreParams, FieldProfile, FieldState, GridPoints, Pump, profile_max_diff,
};
use laser_solver::picard::{PicardConfig, PicardDfbSolver, dfb_solve_picard, initial_profile};
use laser_solver::rootfind::{BisectionConfig, Midpoint, Newton1dConfig, rootfind_1d};
use laser_solver::utils::{IterationConfig, linspace};
use myplotlib::Plotter;
use plots::show_field_profile;
use std::time::Instant;

const PUMP_FORWARD_AMPLITUDE: f64 = 100.0;
const PUMP_BACKWARD_AMPLITUDE: f64 = 10.0;
const PUMP: Pump = Pump {
    total: PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
        + PUMP_BACKWARD_AMPLITUDE * PUMP_BACKWARD_AMPLITUDE,
    balance: (PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
        - PUMP_BACKWARD_AMPLITUDE * PUMP_BACKWARD_AMPLITUDE)
        / (PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
            + PUMP_BACKWARD_AMPLITUDE * PUMP_BACKWARD_AMPLITUDE),
};
const FORWARD_PUMP: Pump = Pump {
    total: PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE,
    balance: 1.0,
};

const FIBRE: FibreParams = FibreParams {
    density: 1.0,
    lifetime: 1.0,
    pump_ab: 0.01 * 100.0,
    pump_em: 0.0,
    sgnl_ab: 0.0,
    sgnl_em: 1.0,
    length: 10.0,
};

const GRID: GridPoints = GridPoints(500);
const FULL_PROFILE: bool = true;

const GRATING: GratingProfile = GratingProfile {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.45,
};

const ITERATION: IterationConfig = IterationConfig {
    max: 500,
    tol: 1e-10,
};

const PICARD: PicardConfig = PicardConfig {
    max_iterations: 500,
    relative_tolerance: 1e-10,
    absolute_tolerance: 1e-12,
};

const BISECTION: BisectionConfig = BisectionConfig {
    iteration: ITERATION,
    upper: PUMP_FORWARD_AMPLITUDE,
    lower: 1e-8,
    midpoint: Midpoint::Geometric,
};

const NEWTON: Newton1dConfig = Newton1dConfig {
    iteration: ITERATION,
    initial: PUMP_FORWARD_AMPLITUDE,
    dx: 1e-6,
};

const SHOW_PLOTS: bool = true;

fn main() -> eframe::Result {
    inspect_resiudal_curve(SHOW_PLOTS)?;
    inspect_field_profiles(SHOW_PLOTS)?;
    run_pump_scan(SHOW_PLOTS)?;
    plot_pump_scan_derivatives(SHOW_PLOTS)?;
    inspect_grating(SHOW_PLOTS)?;
    compare_profile_solvers(SHOW_PLOTS)?;
    compare_dfb_solvers(SHOW_PLOTS)?;

    Ok(())
}
fn inspect_resiudal_curve(show_plots: bool) -> eframe::Result {
    let kappas = GRATING.grid(GRID.0);
    let dz = GRID.dz(FIBRE.length);
    let trial = |sgnl_b| FieldState {
        sgnl_f: 0.0,
        sgnl_b,
        pump_f: 2.0,
        pump_b: 0.0, // shooting method requires zero backward pump amplitude
    };
    let f = |sgnl_b| out_field(trial(sgnl_b), FIBRE, dz, &kappas).sgnl_b / sgnl_b;
    let root = rootfind_1d(f, BISECTION).expect("root not found");
    println!("root is at {}", root);
    println!("residual at 0 {}", f(0.0));

    let sgnl_bs = linspace(1e-8, 5.0 * root, 1000);
    let residuals: Vec<f64> = sgnl_bs.iter().map(|&s| f(s).abs().log10()).collect();
    if show_plots {
        let mut plot = Plotter::new();
        plot.plot(&sgnl_bs, &residuals);
        plot.title("Residuals");
        plot.show()?;
    }
    Ok(())
}
fn inspect_field_profiles(show_plots: bool) -> eframe::Result {
    let result =
        dfb_solve_shooting(FORWARD_PUMP, FIBRE, GRID, GRATING, FULL_PROFILE, NEWTON).unwrap();
    show_field_profile(&result, show_plots)?;

    let profile = initial_profile(PUMP, FIBRE, GRID);
    show_field_profile(&profile, show_plots)?;

    Ok(())
}

fn run_pump_scan(show_plots: bool) -> eframe::Result {
    let pumps = linspace(0.0, 10.0, 200);
    let start = Instant::now();
    let threshold = dfb_pump_scan(&pumps, 1.0, FIBRE, GRID, GRATING, BISECTION, PICARD)
        .expect("pump scan failed");
    let elapsed = start.elapsed();
    println!("pump sweep {:.3}", elapsed.as_secs_f64());

    if show_plots {
        let sgnl_f: Vec<f64> = threshold
            .iter()
            .map(|output| output.as_ref().map_or(0.0, |output| output.0))
            .collect();
        let sgnl_b: Vec<f64> = threshold
            .iter()
            .map(|output| output.as_ref().map_or(0.0, |output| output.1))
            .collect();

        let mut plot = Plotter::new();
        plot.plot(&pumps, &sgnl_f).label("Forward");
        plot.plot(&pumps, &sgnl_b).label("Back");
        plot.title("threshold shooting");
        plot.show()?;
    }

    Ok(())
}

fn plot_pump_scan_derivatives(show_plot: bool) -> eframe::Result {
    if !show_plot {
        return Ok(());
    }
    let pc = PicardConfig {
        max_iterations: 500,
        relative_tolerance: 1e-6,
        absolute_tolerance: 1e-6,
    };
    let pumps = linspace(0.0, 10.0, 200);
    let balance = 0.95;
    let outputs = dfb_pump_scan(&pumps, balance, FIBRE, GRID, GRATING, BISECTION, pc)
        .expect("pump scan failed");
    let outputs: Vec<(f64, f64)> = outputs
        .into_iter()
        .map(|output| output.unwrap_or((0.0, 0.0)))
        .collect();

    let derivative_pumps: Vec<f64> = pumps
        .windows(2)
        .map(|window| 0.5 * (window[0] + window[1]))
        .collect();
    let forward_derivative: Vec<f64> = pumps
        .windows(2)
        .zip(outputs.windows(2))
        .map(|(pump, output)| (output[1].0 - output[0].0) / (pump[1] - pump[0]))
        .collect();
    let backward_derivative: Vec<f64> = pumps
        .windows(2)
        .zip(outputs.windows(2))
        .map(|(pump, output)| (output[1].1 - output[0].1) / (pump[1] - pump[0]))
        .collect();

    let threshold_config = IterationConfig { tol: 1e-3, max: 20 };
    let (forward_slope, backward_slope, threshold) = dfb_find_threshold_and_slope(
        Pump {
            total: 2.0,
            balance,
        },
        1.0,
        threshold_config,
        FIBRE,
        GRID,
        GRATING,
        BISECTION,
        pc,
    )
    .expect("threshold not found");

    println!(
        "forward slope {forward_slope}, backward slope {backward_slope}, threshold {threshold}"
    );

    let mut plot = Plotter::new();
    plot.plot(&derivative_pumps, &forward_derivative)
        .label("Forward derivative");
    plot.plot(&derivative_pumps, &backward_derivative)
        .label("Backward derivative");
    plot.axhline(forward_slope).label("Forward slope");
    plot.axhline(backward_slope).label("Backward slope");
    plot.axvline(threshold).label("Threshold");
    plot.xlabel("Pump power");
    plot.ylabel("d(output power) / d(pump power)");
    plot.title("Pump-scan slope efficiency");
    plot.show()
}

fn inspect_grating(show_plot: bool) -> eframe::Result {
    if !show_plot {
        return Ok(());
    }

    let z = GRID.grid(FIBRE.length);
    let kappas = GRATING.grid(GRID.0 + 1);
    let mut plot = Plotter::new();
    plot.plot(&z, &kappas);
    plot.xlabel("z");
    plot.ylabel("Kappa");
    plot.title("Coupling Profile");
    plot.show()
}

fn compare_profile_solvers(show_plots: bool) -> eframe::Result {
    let comparison_pump = FORWARD_PUMP;
    let comparison_sgnl_b = 1.0;
    let comparison_kappas = GRATING.grid(GRID.0);
    let comparison_boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b: comparison_sgnl_b,
        pump_f: comparison_pump.forward_amplitude(),
        pump_b: 0.0,
    };

    let direct_profile = FieldProfile::new(
        GRID.grid(FIBRE.length),
        solve_profile(
            comparison_boundary,
            FIBRE,
            GRID.dz(FIBRE.length),
            &comparison_kappas,
        ),
    );

    let mut picard_solver = PicardDfbSolver::new(comparison_pump, FIBRE, GRID);
    let picard_fields = picard_solver
        .solve_profile_picard(
            comparison_sgnl_b,
            comparison_pump,
            FIBRE,
            PICARD,
            &comparison_kappas,
            GRID.dz(FIBRE.length),
        )
        .expect("Picard profile comparison did not converge")
        .to_vec();
    let picard_profile = FieldProfile::new(direct_profile.z.clone(), picard_fields);

    let max_diff = profile_max_diff(&direct_profile.fields, &picard_profile.fields);
    println!("Picard/direct profile max diff: {max_diff:e}");

    if show_plots {
        plot_profile_diff(
            &direct_profile,
            &picard_profile,
            "shooting vs picard profile",
        )?;
    }

    Ok(())
}

fn compare_dfb_solvers(show_plots: bool) -> eframe::Result {
    let comparison_pump = FORWARD_PUMP;

    let start = Instant::now();
    let shooting_profile =
        dfb_solve_shooting(comparison_pump, FIBRE, GRID, GRATING, FULL_PROFILE, NEWTON)
            .expect("shooting DFB solve failed");
    let shooting_elapsed = start.elapsed();

    let start = Instant::now();
    let picard_profile = dfb_solve_picard(
        comparison_pump,
        FIBRE,
        GRID,
        GRATING,
        FULL_PROFILE,
        NEWTON,
        PICARD,
    )
    .expect("Picard DFB solve failed");
    let picard_elapsed = start.elapsed();

    let max_diff = profile_max_diff(&shooting_profile.fields, &picard_profile.fields);
    println!(
        "shooting DFB solve: {:.3} ms",
        shooting_elapsed.as_secs_f64() * 1_000.0
    );
    println!(
        "Picard DFB solve: {:.3} ms",
        picard_elapsed.as_secs_f64() * 1_000.0
    );
    println!("shooting/Picard profile max diff: {max_diff:e}");

    if show_plots {
        plot_profile_diff(
            &shooting_profile,
            &picard_profile,
            "shooting vs picard DFB solution",
        )?;
    }

    Ok(())
}

fn difference(left: &[f64], right: &[f64]) -> Vec<f64> {
    left.iter()
        .zip(right)
        .map(|(left, right)| left - right)
        .collect()
}
