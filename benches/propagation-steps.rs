use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
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
const KAPPA: f64 = 1.0;
const BATCH_SIZE: usize = 256;

#[derive(Copy, Clone)]
enum Coupling {
    None,
    Bragg(f64),
}

#[derive(Copy, Clone)]
enum SectionCoupling<'a> {
    None(usize),
    Bragg(&'a [f64]),
}

fn step_if(fields: FieldState, gain: Gain, kappa: f64, dz: f64) -> FieldState {
    if kappa == 0.0 {
        fields.uncoupled_step(gain, dz)
    } else {
        fields.coupled_step(gain, kappa, dz)
    }
}

fn step_match(fields: FieldState, gain: Gain, coupling: Coupling, dz: f64) -> FieldState {
    match coupling {
        Coupling::None => fields.uncoupled_step(gain, dz),
        Coupling::Bragg(kappa) => fields.coupled_step(gain, kappa, dz),
    }
}

fn blocked_pattern() -> Vec<f64> {
    (0..BATCH_SIZE)
        .map(|i| if i < BATCH_SIZE / 2 { 0.0 } else { KAPPA })
        .collect()
}

fn shuffled_pattern() -> Vec<f64> {
    let mut pattern = blocked_pattern();
    let mut state = 0x9e37_79b9_u32;
    for i in (1..pattern.len()).rev() {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        pattern.swap(i, state as usize % (i + 1));
    }
    pattern
}

fn as_coupling(pattern: &[f64]) -> Vec<Coupling> {
    pattern
        .iter()
        .map(|&kappa| {
            if kappa == 0.0 {
                Coupling::None
            } else {
                Coupling::Bragg(kappa)
            }
        })
        .collect()
}

fn run_if(pattern: &[f64]) -> FieldState {
    pattern.iter().fold(FIELDS, |fields, &kappa| {
        step_if(fields, black_box(GAIN), kappa, black_box(DZ))
    })
}

fn run_match(pattern: &[Coupling]) -> FieldState {
    pattern.iter().fold(FIELDS, |fields, &coupling| {
        step_match(fields, black_box(GAIN), coupling, black_box(DZ))
    })
}

fn run_direct_split(kappas: &[f64]) -> FieldState {
    let fields = (0..BATCH_SIZE / 2).fold(FIELDS, |fields, _| {
        fields.uncoupled_step(black_box(GAIN), black_box(DZ))
    });
    kappas.iter().fold(fields, |fields, &kappa| {
        fields.coupled_step(black_box(GAIN), kappa, black_box(DZ))
    })
}

fn run_section(fields: FieldState, coupling: SectionCoupling<'_>) -> FieldState {
    match coupling {
        SectionCoupling::None(nsteps) => (0..nsteps).fold(fields, |fields, _| {
            fields.uncoupled_step(black_box(GAIN), black_box(DZ))
        }),
        SectionCoupling::Bragg(kappas) => kappas.iter().fold(fields, |fields, &kappa| {
            fields.coupled_step(black_box(GAIN), kappa, black_box(DZ))
        }),
    }
}

fn run_outer_match(kappas: &[f64]) -> FieldState {
    let fields = run_section(
        FIELDS,
        black_box(SectionCoupling::None(BATCH_SIZE / 2)),
    );
    run_section(fields, black_box(SectionCoupling::Bragg(kappas)))
}

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

    group.bench_function("if_zero", |b| {
        b.iter(|| {
            black_box(step_if(
                black_box(FIELDS),
                black_box(GAIN),
                black_box(0.0),
                black_box(DZ),
            ))
        });
    });

    group.bench_function("match_none", |b| {
        b.iter(|| {
            black_box(step_match(
                black_box(FIELDS),
                black_box(GAIN),
                black_box(Coupling::None),
                black_box(DZ),
            ))
        });
    });

    group.finish();
}

fn benchmark_branch_patterns(c: &mut Criterion) {
    let blocked = blocked_pattern();
    let shuffled = shuffled_pattern();
    let blocked_match = as_coupling(&blocked);
    let shuffled_match = as_coupling(&shuffled);
    let bragg_section = vec![KAPPA; BATCH_SIZE / 2];
    let mut group = c.benchmark_group("propagation/branch_pattern");
    group.throughput(Throughput::Elements(BATCH_SIZE as u64));

    group.bench_function("direct_split", |b| {
        b.iter(|| black_box(run_direct_split(black_box(&bragg_section))));
    });

    group.bench_function("outer_match_sections", |b| {
        b.iter(|| black_box(run_outer_match(black_box(&bragg_section))));
    });

    for (name, if_pattern, match_pattern) in [
        ("blocked", &blocked, &blocked_match),
        ("shuffled", &shuffled, &shuffled_match),
    ] {
        group.bench_with_input(
            BenchmarkId::new("if", name),
            if_pattern,
            |b, pattern| b.iter(|| black_box(run_if(black_box(pattern)))),
        );
        group.bench_with_input(
            BenchmarkId::new("match", name),
            match_pattern,
            |b, pattern| b.iter(|| black_box(run_match(black_box(pattern)))),
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_propagation_steps,
    benchmark_branch_patterns
);
criterion_main!(benches);
