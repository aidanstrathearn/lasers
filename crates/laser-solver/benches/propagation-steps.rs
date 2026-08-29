use criterion::{Criterion, criterion_group, criterion_main};
use laser_solver::lase::{BidirectionalAmplitude, FieldState, Gain};
use std::hint::black_box;

const FIELDS: FieldState = FieldState {
    signal: BidirectionalAmplitude {
        forward: 2.0,
        backward: 3.0,
    },
    pump: BidirectionalAmplitude {
        forward: 5.0,
        backward: 7.0,
    },
};

const GAIN: Gain = Gain {
    pump: -0.2,
    signal: 0.1,
};

const DZ: f64 = 0.01;

fn benchmark_propagation_steps(c: &mut Criterion) {
    let mut group = c.benchmark_group("propagation/step");

    group.bench_function("uncoupled", |b| {
        b.iter(|| {
            black_box(black_box(FIELDS).uncoupled_step(black_box(GAIN), black_box(DZ)))
        });
    });

    group.bench_function("coupled_kappa_zero", |b| {
        b.iter(|| {
            black_box(black_box(FIELDS).coupled_step(
                black_box(GAIN),
                black_box(0.0),
                black_box(DZ),
            ))
        });
    });

    group.finish();
}

criterion_group!(benches, benchmark_propagation_steps);
criterion_main!(benches);
