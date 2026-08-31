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

const MIRRORED_PUMP_FLUX: f64 = PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE;
const SYMMETRY_ABSOLUTE_TOLERANCE: f64 = 1e-8;
const SYMMETRY_RELATIVE_TOLERANCE: f64 = 5e-3;

#[test]
fn injected_solver_satisfies_active_fibre_boundaries() {
    let fibre = resolved_fibre();
    let signal = Signal {
        total: fibre.signal_power(1.0),
        balance: 0.0,
    };
    let pump = Pump {
        total: fibre.pump_power(1.0),
        balance: 0.0,
    };

    let profile = TwoModeSolver::new(&fibre, STEPS)
        .solve_injected(
            pump,
            signal,
            RootFindConfig::Bisection(BISECTION),
            SYMMETRY_PICARD,
        )
        .expect("active-fibre injected solve failed");
    let left_signal_forward = profile.signal_forward_power().next().unwrap();
    let right_signal_backward = profile.signal_backward_power().last().unwrap();
    let left_pump_forward = profile.pump_forward_power().next().unwrap();
    let right_pump_backward = profile.pump_backward_power().last().unwrap();

    assert_close(
        0,
        "left signal input",
        left_signal_forward,
        0.5 * signal.total,
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert_close(
        profile.len() - 1,
        "right signal input",
        right_signal_backward,
        0.5 * signal.total,
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert_close(
        0,
        "left pump input",
        left_pump_forward,
        0.5 * pump.total,
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert_close(
        profile.len() - 1,
        "right pump input",
        right_pump_backward,
        0.5 * pump.total,
        SYMMETRY_ABSOLUTE_TOLERANCE,
        SYMMETRY_RELATIVE_TOLERANCE,
    );
    assert!(
        profile
            .signal_forward_power()
            .chain(profile.signal_backward_power())
            .chain(profile.pump_forward_power())
            .chain(profile.pump_backward_power())
            .all(f64::is_finite),
        "active-fibre injected solve returned non-finite fields"
    );
}

#[test]
fn pump_scan_matches_independent_lasing_solves() {
    let fibre = resolved_dfb_fibre();
    let pump_totals = [2.0, 4.0, 6.0, 8.0, 10.0].map(|flux| fibre.pump_power(flux));
    let balance = 0.95;
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
    let fibre = resolved_symmetric_dfb_fibre();
    let pump_power = fibre.pump_power(MIRRORED_PUMP_FLUX);
    let shooting_pump = Pump {
        total: pump_power,
        balance: 1.0,
    };
    let picard_pump = Pump {
        total: pump_power,
        balance: -1.0,
    };

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
        .signal_forward_power()
        .chain(profile.signal_backward_power())
        .fold(0.0_f64, f64::max);
    assert!(
        max_signal > 0.0,
        "{label} solver converged to the trivial zero-signal solution"
    );
}

fn assert_mirrored_profiles_agree(picard: &FieldProfile, shooting: &FieldProfile) {
    assert_eq!(picard.len(), shooting.len());

    let picard_z = picard.z().collect::<Vec<_>>();
    let shooting_z = shooting.z().collect::<Vec<_>>();
    let picard_signal_forward = picard.signal_forward_power().collect::<Vec<_>>();
    let picard_signal_backward = picard.signal_backward_power().collect::<Vec<_>>();
    let picard_pump_forward = picard.pump_forward_power().collect::<Vec<_>>();
    let picard_pump_backward = picard.pump_backward_power().collect::<Vec<_>>();
    let shooting_signal_forward = shooting.signal_forward_power().collect::<Vec<_>>();
    let shooting_signal_backward = shooting.signal_backward_power().collect::<Vec<_>>();
    let shooting_pump_forward = shooting.pump_forward_power().collect::<Vec<_>>();
    let shooting_pump_backward = shooting.pump_backward_power().collect::<Vec<_>>();

    for index in 0..picard.len() {
        let reverse = shooting.len() - 1 - index;
        assert_close(
            index,
            "z",
            picard_z[index],
            FIBRE.geometry.length - shooting_z[reverse],
            1e-12,
            0.0,
        );
        assert_close(
            index,
            "sgnl_f",
            picard_signal_forward[index],
            shooting_signal_backward[reverse],
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "sgnl_b",
            picard_signal_backward[index],
            shooting_signal_forward[reverse],
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "pump_f",
            picard_pump_forward[index],
            shooting_pump_backward[reverse],
            SYMMETRY_ABSOLUTE_TOLERANCE,
            SYMMETRY_RELATIVE_TOLERANCE,
        );
        assert_close(
            index,
            "pump_b",
            picard_pump_backward[index],
            shooting_pump_forward[reverse],
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
