use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use laser_solver::dfb::{DfbLaser, DfbSolveConfig, Grating};
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, Pump, ResolvedFibre, TwoLevelCrossSections,
    TwoLevelDopant,
};
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::{BisectionConfig, Midpoint, Newton1dConfig, RootFindConfig};
use laser_solver::maths::utils::IterationConfig;
use std::hint::black_box;

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

fn benchmark_dfb_solver(c: &mut Criterion) {
    c.bench_function("dfb/shooting", |b| {
        b.iter(|| {
            let result = dfb_laser()
                .solve_shooting(FORWARD_PUMP, NEWTON_SOLVE_CONFIG, FULL_PROFILE)
                .expect("shooting DFB solve failed");
            black_box(result);
        });
    });
}

fn benchmark_bisection_midpoints(c: &mut Criterion) {
    let mut group = c.benchmark_group("dfb/shooting/bisection");

    for (name, midpoint) in [
        ("arithmetic", Midpoint::Arithmetic),
        ("geometric", Midpoint::Geometric),
        ("interp", Midpoint::Interp),
    ] {
        group.bench_with_input(
            BenchmarkId::from_parameter(name),
            &midpoint,
            |b, &midpoint| {
                let config = BisectionConfig {
                    // Linear interpolation converges much more slowly for this DFB residual.
                    // Keep the tolerance identical while allowing every variant to finish.
                    iteration: IterationConfig {
                        max: 100_000,
                        ..ITERATION
                    },
                    upper: PUMP_FORWARD_AMPLITUDE,
                    lower: 1e-8,
                    midpoint,
                };

                b.iter(|| {
                    let result = dfb_laser()
                        .solve_shooting(
                            FORWARD_PUMP,
                            DfbSolveConfig {
                                steps: STEPS,
                                root_find: RootFindConfig::Bisection(config),
                                picard: PICARD,
                            },
                            FULL_PROFILE,
                        )
                        .expect("shooting DFB solve failed");
                    black_box(result);
                });
            },
        );
    }

    group.finish();
}

fn benchmark_picard_solvers(c: &mut Criterion) {
    let backward_amplitude = 20.0;
    let pump = Pump {
        total: PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
            + backward_amplitude * backward_amplitude,
        balance: (PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
            - backward_amplitude * backward_amplitude)
            / (PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
                + backward_amplitude * backward_amplitude),
    };
    let mut group = c.benchmark_group("dfb/picard");

    group.bench_function("buffered", |b| {
        b.iter(|| {
            let result = dfb_laser()
                .solve_picard(pump, NEWTON_SOLVE_CONFIG, FULL_PROFILE)
                .expect("buffered Picard DFB solve failed");
            black_box(result);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_dfb_solver,
    benchmark_bisection_midpoints,
    benchmark_picard_solvers
);
criterion_main!(benches);
