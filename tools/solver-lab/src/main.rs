mod myplotlib;
mod plots;

use laser_solver::dfb::{dfb_pump_scan_shooting, dfb_solve_shooting, solve_profile, transfer};
use laser_solver::lase::{
    FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump, profile_max_diff,
};
use laser_solver::picard::{PicardConfig, dfb_solve_picard, initial_profile, solve_profile_picard};
use laser_solver::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
use laser_solver::utils::{IterationConfig, geomspace};
use myplotlib::Plotter;
use plots::show_field_profile;
use std::time::Instant;
use crate::plots::plot_profile_diff;

const PUMP: Pump = Pump {
    forward: 100.0,
    backward: 10.0,
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
    upper: PUMP.forward,
    lower: 1e-8,
    midpoint: Midpoint::Geometric,
};

const NEWTON: Newton1dConfig = Newton1dConfig {
    iteration: ITERATION,
    initial: PUMP.forward,
    dx: 1e-6,
};

const SHOW_PLOTS: bool = true;

fn main() -> eframe::Result {
    inspect_field_profiles(SHOW_PLOTS)?;
    run_pump_scan(SHOW_PLOTS)?;
    inspect_grating(SHOW_PLOTS)?;
    compare_profile_solvers(SHOW_PLOTS)?;
    compare_dfb_solvers(SHOW_PLOTS)?;

    Ok(())
}

fn inspect_field_profiles(show_plots: bool) -> eframe::Result {
    let result = dfb_solve_shooting(PUMP, FIBRE, GRID, GRATING, FULL_PROFILE, NEWTON).unwrap();
    show_field_profile(&result, show_plots)?;

    let profile = initial_profile(PUMP, FIBRE, GRID);
    show_field_profile(&profile, show_plots)?;

    Ok(())
}

fn run_pump_scan(show_plots: bool) -> eframe::Result {
    let pumps = geomspace(-1.0, 1.0, 200);
    let start = Instant::now();
    let threshold = dfb_pump_scan_shooting(&pumps, 1.0, FIBRE, GRID, GRATING, BISECTION, PICARD);
    let elapsed = start.elapsed();

    println!("pump sweep {:.3}", elapsed.as_secs_f64());

    if show_plots {
        let pump_log10: Vec<f64> = pumps.iter().map(|x| x.log10()).collect();
        let sgnl_f: Vec<f64> = threshold.iter().map(|x| x.0).collect();
        let sgnl_b: Vec<f64> = threshold.iter().map(|x| x.1).collect();

        let mut plot = Plotter::new();
        plot.plot(&pump_log10, &sgnl_f).label("Forward");
        plot.plot(&pump_log10, &sgnl_b).label("Back");
        plot.show()?;
    }

    Ok(())
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
    let comparison_pump = Pump {
        backward: 0.0,
        ..PUMP
    };
    let comparison_sgnl_b = 1.0;
    let comparison_kappas = GRATING.grid(GRID.0);
    let comparison_boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b: comparison_sgnl_b,
        pump_f: comparison_pump.forward,
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

    let current = initial_profile(comparison_pump, FIBRE, GRID);
    let picard_fields = solve_profile_picard(
        comparison_sgnl_b,
        current.fields.clone(),
        comparison_pump,
        FIBRE,
        PICARD,
        &comparison_kappas,
        GRID.dz(FIBRE.length),
    )
    .expect("Picard profile comparison did not converge");
    let picard_profile = FieldProfile::new(direct_profile.z.clone(), picard_fields);

    let max_diff = profile_max_diff(&direct_profile.fields, &picard_profile.fields);
    println!("Picard/direct profile max diff: {max_diff:e}");

    if show_plots {
        plot_profile_diff(&direct_profile, &picard_profile, "shooting vs picard profile")?;
    }

    Ok(())
}

fn compare_dfb_solvers(show_plots: bool) -> eframe::Result {
    let comparison_pump = Pump {
        backward: 0.0,
        ..PUMP
    };

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
        plot_profile_diff(&shooting_profile, &picard_profile, "shooting vs picard DFB solution")?;
    }

    Ok(())
}


fn difference(left: &[f64], right: &[f64]) -> Vec<f64> {
    left.iter()
        .zip(right)
        .map(|(left, right)| left - right)
        .collect()
}
