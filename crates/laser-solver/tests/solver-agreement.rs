use laser_solver::dfb::{dfb_solve_shooting, solve_profile};
use laser_solver::lase::{
    FibreParams, FieldProfile, FieldState, GratingProfile, GridPoints, Pump, profile_max_diff,
};
use laser_solver::picard::{
    PicardConfig, PicardDfbSolver, dfb_solve_picard, initial_profile, solve_profile_picard,
};
use laser_solver::rootfind::Newton1dConfig;
use laser_solver::utils::IterationConfig;

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

const NEWTON: Newton1dConfig = Newton1dConfig {
    iteration: ITERATION,
    initial: PUMP.forward,
    dx: 1e-6,
};

const MAX_RELATIVE_DIFFERENCE: f64 = 1e-16;

#[test]
fn direct_and_picard_profile_solvers_agree() {
    let pump = Pump {
        backward: 0.0,
        ..PUMP
    };
    let sgnl_b = 1.0;
    let kappas = GRATING.grid(GRID.0);
    let boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b,
        pump_f: pump.forward,
        pump_b: 0.0,
    };

    let direct_profile = FieldProfile::new(
        GRID.grid(FIBRE.length),
        solve_profile(boundary, FIBRE, GRID.dz(FIBRE.length), &kappas),
    );
    let picard_fields = solve_profile_picard(
        sgnl_b,
        initial_profile(pump, FIBRE, GRID).fields,
        pump,
        FIBRE,
        PICARD,
        &kappas,
        GRID.dz(FIBRE.length),
    )
    .expect("Picard profile solve failed");
    let picard_profile = FieldProfile::new(direct_profile.z.clone(), picard_fields);

    assert_profiles_agree(
        "direct and Picard profile solvers",
        &direct_profile,
        &picard_profile,
    );
}

#[test]
fn direct_and_buffered_picard_profile_solvers_agree() {
    let pump = Pump {
        backward: 0.0,
        ..PUMP
    };
    let sgnl_b = 1.0;
    let kappas = GRATING.grid(GRID.0);
    let boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b,
        pump_f: pump.forward,
        pump_b: 0.0,
    };

    let direct_profile = FieldProfile::new(
        GRID.grid(FIBRE.length),
        solve_profile(boundary, FIBRE, GRID.dz(FIBRE.length), &kappas),
    );
    let mut picard_solver = PicardDfbSolver::init(initial_profile(pump, FIBRE, GRID).fields);
    let picard_fields = picard_solver
        .solve_profile_picard(sgnl_b, pump, FIBRE, PICARD, &kappas, GRID.dz(FIBRE.length))
        .expect("buffered Picard profile solve failed")
        .to_vec();
    let picard_profile = FieldProfile::new(direct_profile.z.clone(), picard_fields);

    assert_profiles_agree(
        "direct and buffered Picard profile solvers",
        &direct_profile,
        &picard_profile,
    );
}

#[test]
fn shooting_and_picard_dfb_solvers_agree() {
    let pump = Pump {
        backward: 0.0,
        ..PUMP
    };

    let shooting_profile = dfb_solve_shooting(pump, FIBRE, GRID, GRATING, true, NEWTON)
        .expect("shooting DFB solve failed");
    let picard_profile = dfb_solve_picard(pump, FIBRE, GRID, GRATING, true, NEWTON, PICARD)
        .expect("Picard DFB solve failed");

    assert_profiles_agree(
        "shooting and Picard DFB solvers",
        &shooting_profile,
        &picard_profile,
    );
}

fn assert_profiles_agree(label: &str, left: &FieldProfile, right: &FieldProfile) {
    assert_eq!(
        left.z, right.z,
        "{label} returned profiles on different grids"
    );
    assert_eq!(
        left.fields.len(),
        right.fields.len(),
        "{label} returned profiles of different lengths"
    );

    let max_diff = profile_max_diff(&left.fields, &right.fields);
    assert!(
        max_diff <= MAX_RELATIVE_DIFFERENCE,
        "{label} differ by {max_diff:e}, exceeding {MAX_RELATIVE_DIFFERENCE:e}"
    );

    assert_eq!(max_diff, 0.0, "{label} differ by {max_diff:e}, not bitwise equal");
}
