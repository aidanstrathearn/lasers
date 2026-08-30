use laser_solver::grating::{NoGrating, PiShift};
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, FieldProfile, Pump, ResolvedFibre, Signal,
    TwoLevelCrossSections, TwoLevelDopant, pump_scan as scan_pump_totals,
};
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::{BisectionConfig, Midpoint, RootFindConfig};
use laser_solver::maths::utils::IterationConfig;
use laser_solver::two_mode::TwoModeSolver;

const PUMP_FORWARD_AMPLITUDE: f64 = 100.0;
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
    )
}

fn resolved_dfb_fibre() -> ResolvedFibre<'static, TwoLevelDopant, PiShift> {
    DFB_FIBRE.resolve_with_interactions(
        PUMP_MODE,
        PUMP_INTERACTION,
        SGNL_MODE,
        SGNL_INTERACTION,
    )
}

fn resolved_symmetric_dfb_fibre() -> ResolvedFibre<'static, TwoLevelDopant, PiShift> {
    SYMMETRIC_DFB_FIBRE.resolve_with_interactions(
        PUMP_MODE,
        PUMP_INTERACTION,
        SGNL_MODE,
        SGNL_INTERACTION,
    )
}

const ITERATION: IterationConfig = IterationConfig {
    max: 500,
    tol: 1e-10,
};

// A nonzero backward pump converges asymptotically rather than reaching a
// bitwise fixed point, so use tolerances appropriate for the convergence metric.
const SYMMETRY_PICARD: PicardConfig = PicardConfig {
    max_iterations: 500,
    relative_tolerance: 1e-6,
    absolute_tolerance: 1e-10,
};

const BISECTION: BisectionConfig = BisectionConfig {
    iteration: ITERATION,
    upper: PUMP_FORWARD_AMPLITUDE,
    lower: 1e-8,
    midpoint: Midpoint::Geometric,
};

const MIRRORED_PUMP: f64 = PUMP_FORWARD_AMPLITUDE;
const SYMMETRY_ABSOLUTE_TOLERANCE: f64 = 1e-8;
const SYMMETRY_RELATIVE_TOLERANCE: f64 = 5e-3;

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

    let profile = TwoModeSolver::new(&fibre, STEPS)
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
fn pump_scan_matches_independent_lasing_solves() {
    let pump_totals = [2.0, 4.0, 6.0, 8.0, 10.0];
    let balance = 0.95;
    let fibre = resolved_dfb_fibre();
    let solver = TwoModeSolver::new(&fibre, STEPS);
    let scan = solver
        .pump_scan(
            &pump_totals,
            balance,
            RootFindConfig::Bisection(BISECTION),
            SYMMETRY_PICARD,
        )
        .expect("pump scan failed");
    let independent = scan_pump_totals(&pump_totals, |total| {
        solver
            .solve_lasing(
                Pump { total, balance },
                RootFindConfig::Bisection(BISECTION),
                SYMMETRY_PICARD,
                false,
            )
            .map(|profile| profile.output_powers())
    })
    .expect("independent lasing solves failed");

    assert!(scan.iter().any(Option::is_some), "scan never reached threshold");
    assert_eq!(scan, independent);
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

    let fibre = resolved_symmetric_dfb_fibre();
    let solver = TwoModeSolver::new(&fibre, SYMMETRY_STEPS);
    let shooting_profile = solver
        .solve_lasing(
            shooting_pump,
            RootFindConfig::Bisection(BISECTION),
            SYMMETRY_PICARD,
            true,
        )
        .expect("forward-pumped shooting DFB solve failed");
    let picard_profile = solver
        .solve_lasing(
            picard_pump,
            RootFindConfig::Bisection(BISECTION),
            SYMMETRY_PICARD,
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
