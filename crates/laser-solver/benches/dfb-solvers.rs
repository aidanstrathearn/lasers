use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use laser_solver::grating::PiShift;
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, Pump, ResolvedFibre, TwoLevelCrossSections,
    TwoLevelDopant,
};
use laser_solver::maths::picard::PicardConfig;
use laser_solver::maths::rootfind::{BisectionConfig, Midpoint, Newton1dConfig, RootFindConfig};
use laser_solver::maths::utils::IterationConfig;
use laser_solver::two_mode::TwoModeSolver;
use std::hint::black_box;

const PUMP_FORWARD_AMPLITUDE: f64 = 100.0;
const FORWARD_PUMP_FLUX: f64 = PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE;

static FIBRE: Fibre<TwoLevelDopant, PiShift> = Fibre {
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
const PUMP_MODE: FieldMode = FieldMode::new(970e-9);
const SGNL_MODE: FieldMode = FieldMode::new(1060e-9);
const PUMP_INTERACTION: TwoLevelCrossSections = TwoLevelCrossSections::new(0.01 * 100.0, 0.0);
const SGNL_INTERACTION: TwoLevelCrossSections = TwoLevelCrossSections::new(0.0, 1.0);

const STEPS: usize = 500;
const FULL_PROFILE: bool = true;

const GRATING: PiShift = PiShift {
    kappa_left: 1.0,
    kappa_right: 1.0,
    pi_shift_position: 0.45,
};

fn resolved_fibre() -> ResolvedFibre<'static, TwoLevelDopant, PiShift> {
    FIBRE.resolve_with_interactions(
        PUMP_MODE,
        PUMP_INTERACTION,
        SGNL_MODE,
        SGNL_INTERACTION,
    )
}

fn pump_for_flux(
    fibre: &ResolvedFibre<'_, TwoLevelDopant, PiShift>,
    total_flux: f64,
    balance: f64,
) -> Pump {
    Pump {
        total: fibre.pump_power(total_flux),
        balance,
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

fn benchmark_dfb_solver(c: &mut Criterion) {
    c.bench_function("dfb/shooting", |b| {
        b.iter(|| {
            let fibre = resolved_fibre();
            let result = TwoModeSolver::new(&fibre, STEPS)
                .solve_lasing(
                    pump_for_flux(&fibre, FORWARD_PUMP_FLUX, 1.0),
                    RootFindConfig::Newton1d(NEWTON),
                    PICARD,
                    FULL_PROFILE,
                )
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
                    let fibre = resolved_fibre();
                    let result = TwoModeSolver::new(&fibre, STEPS)
                        .solve_lasing(
                            pump_for_flux(&fibre, FORWARD_PUMP_FLUX, 1.0),
                            RootFindConfig::Bisection(config),
                            PICARD,
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
    let total_flux = PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
        + backward_amplitude * backward_amplitude;
    let balance = (PUMP_FORWARD_AMPLITUDE * PUMP_FORWARD_AMPLITUDE
        - backward_amplitude * backward_amplitude)
        / total_flux;
    let mut group = c.benchmark_group("dfb/picard");

    group.bench_function("lasing", |b| {
        b.iter(|| {
            let fibre = resolved_fibre();
            let pump = pump_for_flux(&fibre, total_flux, balance);
            let result = TwoModeSolver::new(&fibre, STEPS)
                .solve_lasing(
                    pump,
                    RootFindConfig::Newton1d(NEWTON),
                    PICARD,
                    FULL_PROFILE,
                )
                .expect("Picard DFB solve failed");
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
