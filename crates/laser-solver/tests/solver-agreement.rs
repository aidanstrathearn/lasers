use laser_solver::amplifier::{AmplifierSolveConfig, solve_amp_picard, solve_shooting};
use laser_solver::dfb::picard::solve_profile_picard;
use laser_solver::dfb::{DfbLaser, DfbSolveConfig, Grating, initial_profile};
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, FieldProfile, FieldState, GridPoints, Pump, ResolvedFibre,
    Signal, TwoLevelDopant, profile_max_diff,
};
use laser_solver::maths::picard::{PicardConfig, PicardSolver};
use laser_solver::maths::rootfind::{
    BisectionConfig, Midpoint, Newton1dConfig, RootFindConfig,
};
use laser_solver::maths::utils::IterationConfig;
use laser_solver::propagation::solve_profile_coupled;

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
        pump_ab: 0.01 * 100.0,
        pump_em: 0.0,
        sgnl_ab: 0.0,
        sgnl_em: 1.0,
    },
};
const PUMP_MODE: FieldMode = FieldMode::new(970e-9);
const SGNL_MODE: FieldMode = FieldMode::new(1060e-9);

const GRID: GridPoints = GridPoints(500);
// Gain is sampled at the left edge of each step, so reversal symmetry converges
// with grid refinement rather than holding bit-for-bit on a coarse grid.
const SYMMETRY_GRID: GridPoints = GridPoints(5_000);

const GRATING: Grating = Grating {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.45,
};

const SYMMETRIC_GRATING: Grating = Grating {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.5,
};

fn resolved_fibre() -> ResolvedFibre<'static> {
    FIBRE.resolve(PUMP_MODE, SGNL_MODE)
}

fn dfb_laser() -> DfbLaser<'static> {
    DfbLaser {
        fibre: resolved_fibre(),
        grating: GRATING,
    }
}

fn symmetric_dfb_laser() -> DfbLaser<'static> {
    DfbLaser {
        fibre: resolved_fibre(),
        grating: SYMMETRIC_GRATING,
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

// A nonzero backward pump converges asymptotically rather than reaching a
// bitwise fixed point, so use tolerances appropriate for the convergence metric.
const SYMMETRY_PICARD: PicardConfig = PicardConfig {
    max_iterations: 500,
    relative_tolerance: 1e-6,
    absolute_tolerance: 1e-10,
};

const NEWTON: Newton1dConfig = Newton1dConfig {
    iteration: ITERATION,
    initial: PUMP_FORWARD_AMPLITUDE,
    dx: 1e-6,
};

const BISECTION: BisectionConfig = BisectionConfig {
    iteration: ITERATION,
    upper: PUMP_FORWARD_AMPLITUDE,
    lower: 1e-8,
    midpoint: Midpoint::Geometric,
};

const NEWTON_SOLVE_CONFIG: DfbSolveConfig = DfbSolveConfig {
    grid_points: GRID,
    root_find: RootFindConfig::Newton1d(NEWTON),
    picard: PICARD,
};

const BISECTION_SOLVE_CONFIG: DfbSolveConfig = DfbSolveConfig {
    grid_points: GRID,
    root_find: RootFindConfig::Bisection(BISECTION),
    picard: PICARD,
};

const MAX_RELATIVE_DIFFERENCE: f64 = 1e-16;
const MIRRORED_PUMP: f64 = PUMP_FORWARD_AMPLITUDE;
const SYMMETRY_ABSOLUTE_TOLERANCE: f64 = 1e-8;
const SYMMETRY_RELATIVE_TOLERANCE: f64 = 5e-3;
const AMPLIFIER_PROFILE_TOLERANCE: f64 = 1e-8;

#[test]
fn shooting_and_picard_amplifier_profiles_agree() {
    let signal = Signal {
        total: 1.0,
        balance: 1.0,
    };
    let pump = Pump {
        total: 1.0,
        balance: 0.0,
    };
    let config = AmplifierSolveConfig {
        grid_points: GRID,
        root_find: RootFindConfig::Bisection(BISECTION),
        picard: SYMMETRY_PICARD,
    };

    let fibre = resolved_fibre();
    let shooting_profile = solve_shooting(&fibre, signal.forward_amplitude(), pump, config, true)
        .expect("shooting amplifier solve failed");
    let picard_profile = solve_amp_picard(&fibre, signal, pump, config, true)
        .expect("Picard amplifier solve failed");

    assert_eq!(shooting_profile.z, picard_profile.z);
    let max_diff = profile_max_diff(&shooting_profile.fields, &picard_profile.fields);
    assert!(
        max_diff <= AMPLIFIER_PROFILE_TOLERANCE,
        "shooting and Picard amplifier profiles differ by {max_diff:e}, exceeding {AMPLIFIER_PROFILE_TOLERANCE:e}"
    );
}

#[test]
fn direct_and_buffered_picard_profile_solvers_agree() {
    let pump = FORWARD_PUMP;
    let sgnl_b = 1.0;
    let kappas = GRATING.grid(GRID.0);
    let boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b,
        pump_f: pump.forward_amplitude(),
        pump_b: 0.0,
    };
    let fibre = resolved_fibre();

    let direct_profile = FieldProfile::new(
        GRID.grid(fibre.length()),
        solve_profile_coupled(
            boundary,
            |fields| fibre.gain(fields),
            GRID.dz(fibre.length()),
            &kappas,
        ),
    );
    let initial = initial_profile(pump, &fibre, GRID);
    let mut picard_solver = PicardSolver::from_initial(initial.fields);
    let picard_fields = solve_profile_picard(
        &mut picard_solver,
        sgnl_b,
        pump,
        &fibre,
        PICARD,
        &kappas,
        GRID.dz(fibre.length()),
    )
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
fn shooting_and_picard_dfb_solvers_agree_newton() {
    let pump = FORWARD_PUMP;

    let laser = dfb_laser();
    let shooting_profile = laser
        .solve_shooting(pump, NEWTON_SOLVE_CONFIG, true)
        .expect("shooting DFB solve failed");
    let picard_profile = laser
        .solve_picard(pump, NEWTON_SOLVE_CONFIG, true)
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

    assert_eq!(
        max_diff, 0.0,
        "{label} differ by {max_diff:e}, not bitwise equal"
    );
}

#[test]
fn shooting_and_picard_dfb_solvers_agree_bisection() {
    let pump = FORWARD_PUMP;

    let laser = dfb_laser();
    let shooting_profile = laser
        .solve_shooting(pump, BISECTION_SOLVE_CONFIG, true)
        .expect("shooting DFB solve failed");
    let picard_profile = laser
        .solve_picard(pump, BISECTION_SOLVE_CONFIG, true)
        .expect("Picard DFB solve failed");

    assert_profiles_agree(
        "shooting and Picard DFB solvers",
        &shooting_profile,
        &picard_profile,
    );
}

#[test]
fn backward_pumped_picard_is_reverse_of_forward_pumped_shooting() {
    let shooting_pump = Pump {
        total: MIRRORED_PUMP * MIRRORED_PUMP,
        balance: 1.0,
    };
    let picard_pump = Pump {
        total: MIRRORED_PUMP * MIRRORED_PUMP,
        balance: -1.0,
    };

    let laser = symmetric_dfb_laser();
    let shooting_profile = laser
        .solve_shooting(
            shooting_pump,
            DfbSolveConfig {
                grid_points: SYMMETRY_GRID,
                root_find: RootFindConfig::Bisection(BISECTION),
                picard: SYMMETRY_PICARD,
            },
            true,
        )
        .expect("forward-pumped shooting DFB solve failed");
    let picard_profile = laser
        .solve_picard(
            picard_pump,
            DfbSolveConfig {
                grid_points: SYMMETRY_GRID,
                root_find: RootFindConfig::Bisection(BISECTION),
                picard: SYMMETRY_PICARD,
            },
            true,
        )
        .expect("backward-pumped Picard DFB solve failed");

    assert_nontrivial_signal("shooting", &shooting_profile);
    assert_nontrivial_signal("Picard", &picard_profile);
    assert_mirrored_profiles_agree(&picard_profile, &shooting_profile);
}

fn assert_nontrivial_signal(label: &str, profile: &FieldProfile) {
    let max_signal = profile
        .fields
        .iter()
        .flat_map(|field| [field.sgnl_f.abs(), field.sgnl_b.abs()])
        .fold(0.0_f64, f64::max);
    assert!(
        max_signal > BISECTION.lower,
        "{label} solver converged to the trivial zero-signal solution"
    );
}

fn assert_mirrored_profiles_agree(picard: &FieldProfile, shooting: &FieldProfile) {
    assert_eq!(picard.fields.len(), shooting.fields.len());

    for (index, ((&picard_z, picard_field), (&shooting_z, shooting_field))) in picard
        .z
        .iter()
        .zip(&picard.fields)
        .zip(shooting.z.iter().rev().zip(shooting.fields.iter().rev()))
        .enumerate()
    {
        assert_close(
            index,
            "z",
            picard_z,
            FIBRE.geometry.length - shooting_z,
            1e-12,
            0.0,
        );
        assert_close(
            index,
            "sgnl_f",
            picard_field.sgnl_f,
            shooting_field.sgnl_b,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "sgnl_b",
            picard_field.sgnl_b,
            shooting_field.sgnl_f,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "pump_f",
            picard_field.pump_f,
            shooting_field.pump_b,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "pump_b",
            picard_field.pump_b,
            shooting_field.pump_f,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
    }
}

fn assert_close(
    index: usize,
    component: &str,
    actual: f64,
    expected: f64,
    absolute_tolerance: f64,
    relative_tolerance: f64,
) {
    let difference = (actual - expected).abs();
    let scale = actual.abs().max(expected.abs());
    let tolerance = absolute_tolerance + relative_tolerance * scale;
    assert!(
        difference <= tolerance,
        "field {index} {component} was {actual:e}, expected {expected:e}; difference {difference:e} exceeds {tolerance:e}"
    );
}
