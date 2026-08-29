mod myplotlib;
mod plots;

use crate::plots::plot_profile_diff;
use laser_solver::dfb::picard::solve_profile_picard;
use laser_solver::dfb::{DfbLaser, DfbSolveConfig, Grating};
use laser_solver::lase::{
    BidirectionalAmplitude, Fibre, FibreGeometry, FieldMode, FieldProfile, FieldState, Pump,
    ResolvedFibre, TwoLevelCrossSections, TwoLevelDopant, UniformGrid, profile_max_diff,
};
use laser_solver::maths::picard::{PicardConfig, PicardSolver};
use laser_solver::maths::rootfind::{
    BisectionConfig, Midpoint, Newton1dConfig, RootFindConfig, rootfind_1d,
};
use laser_solver::maths::utils::{IterationConfig, linspace};
use laser_solver::propagation::{out_field_coupled, solve_profile_coupled};
use myplotlib::Plotter;
use plots::show_field_profile;
use std::time::Instant;

const PUMP_FORWARD_AMPLITUDE: f64 = 100.0;
const FORWARD_PUMP: Pump = Pump {
    total: PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE,
    balance: 1.0,
};

static FIBRE: Fibre = Fibre {
    geometry: FibreGeometry {
        core_radius: 4e-6,
        numerical_aperture: 0.1,
        length: 10.0,
    },
    dopant: TwoLevelDopant {
        density: 1.0,
        lifetime: 1.0,
    },
};
const PUMP_MODE: FieldMode = FieldMode::new(970e-9);
const SGNL_MODE: FieldMode = FieldMode::new(1060e-9);
const PUMP_INTERACTION: TwoLevelCrossSections = TwoLevelCrossSections::new(0.01 * 100.0, 0.0);
const SGNL_INTERACTION: TwoLevelCrossSections = TwoLevelCrossSections::new(0.0, 1.0);

const STEPS: usize = 500;
const FULL_PROFILE: bool = true;

const GRATING: Grating = Grating {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.45,
};

fn resolved_fibre() -> ResolvedFibre<'static> {
    FIBRE.resolve_with_interactions(PUMP_MODE, PUMP_INTERACTION, SGNL_MODE, SGNL_INTERACTION)
}

fn dfb_laser() -> DfbLaser<'static> {
    DfbLaser {
        fibre: resolved_fibre(),
        grating: GRATING,
    }
}

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

const NEWTON_SOLVE_CONFIG: DfbSolveConfig = DfbSolveConfig {
    steps: STEPS,
    root_find: RootFindConfig::Newton1d(NEWTON),
    picard: PICARD,
};

const BISECTION_SOLVE_CONFIG: DfbSolveConfig = DfbSolveConfig {
    steps: STEPS,
    root_find: RootFindConfig::Bisection(BISECTION),
    picard: PICARD,
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
    let fibre = resolved_fibre();
    let grid = UniformGrid::new(fibre.length(), STEPS);
    let kappas = GRATING.grid(grid.steps());
    let dz = grid.dz();
    let trial = |sgnl_b| FieldState {
        signal: BidirectionalAmplitude {
            forward: 0.0,
            backward: sgnl_b,
        },
        pump: BidirectionalAmplitude {
            forward: 2.0,
            backward: 0.0, // shooting method requires zero backward pump amplitude
        },
    };
    let f = |sgnl_b| {
        out_field_coupled(trial(sgnl_b), |fields| fibre.gain(fields), dz, &kappas)
            .signal
            .backward
            / sgnl_b
    };
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
    let result = dfb_laser()
        .solve_shooting(FORWARD_PUMP, NEWTON_SOLVE_CONFIG, FULL_PROFILE)
        .unwrap();
    show_field_profile(&result, show_plots)?;

    Ok(())
}

fn run_pump_scan(show_plots: bool) -> eframe::Result {
    let pumps = linspace(0.0, 10.0, 200);
    let start = Instant::now();
    let threshold = dfb_laser()
        .pump_scan(&pumps, 1.0, BISECTION_SOLVE_CONFIG)
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
    let solve_config = DfbSolveConfig {
        picard: pc,
        ..BISECTION_SOLVE_CONFIG
    };
    let outputs = dfb_laser()
        .pump_scan(&pumps, balance, solve_config)
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
    let (forward_slope, backward_slope, threshold) = dfb_laser()
        .find_threshold_and_slope(
            Pump {
                total: 2.0,
                balance,
            },
            1.0,
            threshold_config,
            solve_config,
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

    let grid = UniformGrid::new(FIBRE.geometry.length, STEPS);
    let z = grid.positions().collect::<Vec<_>>();
    let kappas = GRATING.grid(grid.points());
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
    let grid = UniformGrid::new(FIBRE.geometry.length, STEPS);
    let comparison_kappas = GRATING.grid(grid.steps());
    let comparison_boundary = FieldState {
        signal: BidirectionalAmplitude {
            forward: 0.0,
            backward: comparison_sgnl_b,
        },
        pump: BidirectionalAmplitude {
            forward: comparison_pump.forward_amplitude(),
            backward: 0.0,
        },
    };
    let fibre = resolved_fibre();

    let direct_profile = FieldProfile::new(
        grid.positions().collect(),
        solve_profile_coupled(
            comparison_boundary,
            |fields| fibre.gain(fields),
            grid.dz(),
            &comparison_kappas,
        ),
    );

    let mut picard_solver = PicardSolver::filled(grid.points(), comparison_boundary);
    let picard_fields = solve_profile_picard(
        &mut picard_solver,
        comparison_sgnl_b,
        comparison_pump,
        &fibre,
        PICARD,
        &comparison_kappas,
        grid.dz(),
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
    let laser = dfb_laser();
    let shooting_profile = laser
        .solve_shooting(comparison_pump, NEWTON_SOLVE_CONFIG, FULL_PROFILE)
        .expect("shooting DFB solve failed");
    let shooting_elapsed = start.elapsed();

    let start = Instant::now();
    let picard_profile = laser
        .solve_picard(comparison_pump, NEWTON_SOLVE_CONFIG, FULL_PROFILE)
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
