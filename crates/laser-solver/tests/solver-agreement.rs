use laser_solver::dfb::picard::solve_profile_picard;
use laser_solver::dfb::{DfbLaser, DfbSolveConfig};
use laser_solver::grating::{NoGrating, PiShift};
use laser_solver::lase::{
    BidirectionalAmplitude, Fibre, FibreGeometry, FieldMode, FieldProfile, FieldState, Pump,
    ResolvedFibre, Signal, TwoLevelCrossSections, TwoLevelDopant, profile_max_diff,
};
use laser_solver::maths::picard::{PicardConfig, PicardSolver};
use laser_solver::maths::rootfind::{BisectionConfig, Midpoint, Newton1dConfig, RootFindConfig};
use laser_solver::maths::utils::IterationConfig;
use laser_solver::two_mode::TwoModeSolver;
use laser_solver::two_mode::amplifier::{
    Amplifier, AmplifierSolveConfig, solve_amp_picard, solve_shooting,
};
use laser_solver::two_mode::propagation::solve_profile_coupled;

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
    grating: NoGrating,
};
const PUMP_MODE: FieldMode = FieldMode::new(970e-9);
const SGNL_MODE: FieldMode = FieldMode::new(1060e-9);
const PUMP_INTERACTION: TwoLevelCrossSections = TwoLevelCrossSections::new(0.01 * 100.0, 0.0);
const SGNL_INTERACTION: TwoLevelCrossSections = TwoLevelCrossSections::new(0.0, 1.0);

const STEPS: usize = 500;
// Gain is sampled at the left edge of each step, so reversal symmetry converges
// with grid refinement rather than holding bit-for-bit on a coarse grid.
const SYMMETRY_STEPS: usize = 5_000;

const GRATING: PiShift = PiShift {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.45,
};

const SYMMETRIC_GRATING: PiShift = PiShift {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.5,
};

static DFB_FIBRE: Fibre<TwoLevelDopant, PiShift> = Fibre {
    geometry: FibreGeometry {
        core_radius: 4e-6,
        numerical_aperture: 0.1,
        length: 10.0,
    },
    dopant: TwoLevelDopant {
        density: 1.0,
        lifetime: 1.0,
    },
    grating: GRATING,
};

static SYMMETRIC_DFB_FIBRE: Fibre<TwoLevelDopant, PiShift> = Fibre {
    geometry: FibreGeometry {
        core_radius: 4e-6,
        numerical_aperture: 0.1,
        length: 10.0,
    },
    dopant: TwoLevelDopant {
        density: 1.0,
        lifetime: 1.0,
    },
    grating: SYMMETRIC_GRATING,
};

fn resolved_fibre() -> ResolvedFibre<'static> {
    FIBRE.resolve_with_interactions(
        PUMP_MODE,
        PUMP_INTERACTION,
        SGNL_MODE,
        SGNL_INTERACTION,
        STEPS,
    )
}

fn dfb_laser() -> DfbLaser<'static> {
    DfbLaser {
        fibre: DFB_FIBRE.resolve_with_interactions(
            PUMP_MODE,
            PUMP_INTERACTION,
            SGNL_MODE,
            SGNL_INTERACTION,
            STEPS,
        ),
    }
}

fn symmetric_dfb_laser() -> DfbLaser<'static> {
    DfbLaser {
        fibre: SYMMETRIC_DFB_FIBRE.resolve_with_interactions(
            PUMP_MODE,
            PUMP_INTERACTION,
            SGNL_MODE,
            SGNL_INTERACTION,
            SYMMETRY_STEPS,
        ),
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
    root_find: RootFindConfig::Newton1d(NEWTON),
    picard: PICARD,
};

const BISECTION_SOLVE_CONFIG: DfbSolveConfig = DfbSolveConfig {
    root_find: RootFindConfig::Bisection(BISECTION),
    picard: PICARD,
};

const MAX_RELATIVE_DIFFERENCE: f64 = 1e-16;
const MIRRORED_PUMP: f64 = PUMP_FORWARD_AMPLITUDE;
const SYMMETRY_ABSOLUTE_TOLERANCE: f64 = 1e-8;
const SYMMETRY_RELATIVE_TOLERANCE: f64 = 5e-3;
const AMPLIFIER_PROFILE_TOLERANCE: f64 = 1e-11;

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
fn injected_solver_agrees_with_existing_amplifier() {
    let signal = Signal {
        total: 1.0,
        balance: 1.0,
    };
    let pump = Pump {
        total: 1.0,
        balance: 0.0,
    };
    let config = AmplifierSolveConfig {
        root_find: RootFindConfig::Bisection(BISECTION),
        picard: SYMMETRY_PICARD,
    };
    let fibre = resolved_fibre();

    let existing_profile = Amplifier {
        fibre: fibre.clone(),
    }
    .solve(signal, pump, config, true)
    .expect("existing amplifier solve failed");
    let injected_profile = TwoModeSolver::new(&fibre)
        .solve_injected(
            pump,
            signal,
            RootFindConfig::Bisection(BISECTION),
            SYMMETRY_PICARD,
        )
        .expect("new injected-signal solve failed");

    assert_eq!(existing_profile.z, injected_profile.z);
    let max_diff = profile_max_diff(&existing_profile.fields, &injected_profile.fields);
    assert!(
        max_diff <= AMPLIFIER_PROFILE_TOLERANCE,
        "existing amplifier and injected solver profiles differ by {max_diff:e}, exceeding {AMPLIFIER_PROFILE_TOLERANCE:e}"
    );
}

#[test]
fn injected_solver_satisfies_active_fibre_boundaries() {
    let signal = Signal {
        total: 1.0,
        balance: 0.0,
    };
    let pump = Pump {
        total: 1.0,
        balance: 0.0,
    };
    let fibre = resolved_fibre();

    let profile = TwoModeSolver::new(&fibre)
        .solve_injected(
            pump,
            signal,
            RootFindConfig::Bisection(BISECTION),
            SYMMETRY_PICARD,
        )
        .expect("active-fibre injected solve failed");
    let left = profile.fields.first().unwrap();
    let right = profile.fields.last().unwrap();

    assert_close(
        0,
        "left signal input",
        left.signal.forward,
        signal.forward_amplitude(),
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert_close(
        profile.fields.len() - 1,
        "right signal input",
        right.signal.backward,
        signal.backward_amplitude(),
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert_close(
        0,
        "left pump input",
        left.pump.forward,
        pump.forward_amplitude(),
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert_close(
        profile.fields.len() - 1,
        "right pump input",
        right.pump.backward,
        pump.backward_amplitude(),
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert!(
        profile.fields.iter().all(|field| {
            field.signal.forward.is_finite()
                && field.signal.backward.is_finite()
                && field.pump.forward.is_finite()
                && field.pump.backward.is_finite()
        }),
        "active-fibre injected solve returned non-finite fields"
    );
}

#[test]
fn direct_and_buffered_picard_profile_solvers_agree() {
    let pump = FORWARD_PUMP;
    let sgnl_b = 1.0;
    let laser = dfb_laser();
    let fibre = &laser.fibre;
    let grid = fibre.grid();
    let kappas = fibre.kappas();
    let boundary = FieldState {
        signal: BidirectionalAmplitude {
            forward: 0.0,
            backward: sgnl_b,
        },
        pump: BidirectionalAmplitude {
            forward: pump.forward_amplitude(),
            backward: 0.0,
        },
    };
    let direct_profile = FieldProfile::new(
        grid.positions().collect(),
        solve_profile_coupled(boundary, |fields| fibre.gain(fields), grid.dz(), kappas),
    );
    let mut picard_solver = PicardSolver::filled(grid.points(), boundary);
    let picard_fields = solve_profile_picard(
        &mut picard_solver,
        sgnl_b,
        pump,
        fibre,
        PICARD,
        kappas,
        grid.dz(),
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
        .flat_map(|field| [field.signal.forward.abs(), field.signal.backward.abs()])
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
            picard_field.signal.forward,
            shooting_field.signal.backward,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "sgnl_b",
            picard_field.signal.backward,
            shooting_field.signal.forward,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "pump_f",
            picard_field.pump.forward,
            shooting_field.pump.backward,
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "pump_b",
            picard_field.pump.backward,
            shooting_field.pump.forward,
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
