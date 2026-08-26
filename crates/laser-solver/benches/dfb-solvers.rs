use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use laser_solver::dfb::{dfb_solve_shooting, GratingProfile};
use laser_solver::lase::{FibreParams, GridPoints, Pump};
use laser_solver::picard::{PicardConfig, dfb_solve_picard};
use laser_solver::rootfind::{BisectionConfig, Midpoint, Newton1dConfig};
use laser_solver::utils::IterationConfig;
use std::hint::black_box;

const PUMP_FORWARD_AMPLITUDE: f64 = 100.0;
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

const NEWTON: Newton1dConfig = Newton1dConfig {
    iteration: ITERATION,
    initial: PUMP_FORWARD_AMPLITUDE,
    dx: 1e-6,
};

fn benchmark_dfb_solver(c: &mut Criterion) {
    c.bench_function("dfb/shooting", |b| {
        b.iter(|| {
            let result =
                dfb_solve_shooting(FORWARD_PUMP, FIBRE, GRID, GRATING, FULL_PROFILE, NEWTON)
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
                    let result = dfb_solve_shooting(
                        FORWARD_PUMP,
                        FIBRE,
                        GRID,
                        GRATING,
                        FULL_PROFILE,
                        config,
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
            let result = dfb_solve_picard(pump, FIBRE, GRID, GRATING, FULL_PROFILE, NEWTON, PICARD)
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
