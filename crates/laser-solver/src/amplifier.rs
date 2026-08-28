use crate::error::SolverError;
use crate::lase::{
    FieldProfile, FieldState, GridPoints, OutputPower, Pump, ResolvedFibre, Signal,
    profile_convergence_error,
};
use crate::maths::picard::{PicardConfig, PicardError, PicardSolver};
use crate::maths::rootfind::{RootFindConfig, rootfind_1d};
use crate::propagation::{out_field_uncoupled, solve_profile_uncoupled};

#[derive(Copy, Clone)]
pub struct AmplifierSolveConfig {
    pub grid_points: GridPoints,
    pub root_find: RootFindConfig,
    pub picard: PicardConfig,
}

#[derive(Clone)]
pub struct Amplifier<'a> {
    pub fibre: ResolvedFibre<'a>,
}

impl Amplifier<'_> {
    pub fn solve(
        &self,
        signal: Signal,
        pump: Pump,
        config: AmplifierSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let (forward_signal, backward_signal) = signal.amplitudes();

        if backward_signal > 0.0 {
            solve_amp_picard(&self.fibre, signal, pump, config, full_profile)
        } else {
            solve_shooting(&self.fibre, forward_signal, pump, config, full_profile)
        }
    }

    pub(crate) fn output_power(
        &self,
        signal: Signal,
        pump: Pump,
        config: AmplifierSolveConfig,
    ) -> Result<OutputPower, SolverError> {
        let profile = self.solve(signal, pump, config, false)?;
        Ok(profile.output_powers())
    }
}

pub fn solve_shooting(
    fibre: &ResolvedFibre<'_>,
    forward_signal: f64,
    pump: Pump,
    config: AmplifierSolveConfig,
    full_profile: bool,
) -> Result<FieldProfile, SolverError> {
    let gp = config.grid_points;
    let nsteps = gp.0;
    let dz = gp.dz(fibre.length());
    let (pump_forward, pump_backward) = pump.amplitudes();
    let trial = |pump_b| FieldState {
        sgnl_f: forward_signal,
        sgnl_b: 0.0,
        pump_f: pump_forward,
        pump_b,
    };
    let pump_b = if pump_backward == 0.0 {
        0.0
    } else {
        let f = |pump_b| {
            out_field_uncoupled(trial(pump_b), |fields| fibre.gain(fields), dz, nsteps).pump_b
                / pump_backward
                - 1.0
        };
        rootfind_1d(f, config.root_find)?
    };

    if full_profile {
        let z = gp.grid(fibre.length());
        let fields =
            solve_profile_uncoupled(trial(pump_b), |fields| fibre.gain(fields), dz, nsteps);
        Ok(FieldProfile::new(z, fields))
    } else {
        let z = vec![0.0_f64, fibre.length()];
        let out_left = trial(pump_b);
        let fields = vec![
            out_left,
            out_field_uncoupled(out_left, |fields| fibre.gain(fields), dz, nsteps),
        ];
        Ok(FieldProfile::new(z, fields))
    }
}

pub fn find_b_fields(
    signal_b_right: f64,
    pump_b_right: f64,
    profile: &[FieldState],
    fp: &ResolvedFibre<'_>,
    dz: f64,
) -> (f64, f64) {
    let (pump_od, signal_od): (f64, f64) =
        profile[..profile.len() - 1]
            .iter()
            .fold((0.0, 0.0), |acc, &field| {
                let gain = fp.gain(field);
                let new = (0.5 * gain.pump * dz, 0.5 * gain.signal * dz);
                (acc.0 + new.0, acc.1 + new.1)
            });

    (
        signal_b_right * signal_od.exp(),
        pump_b_right * pump_od.exp(),
    )
}

pub fn solve_amp_profile_picard<'a>(
    solver: &'a mut PicardSolver<FieldState>,
    signal: Signal,
    pump: Pump,
    fp: &ResolvedFibre<'_>,
    config: PicardConfig,
    dz: f64,
) -> Result<&'a [FieldState], PicardError> {
    let (pump_forward, pump_backward) = pump.amplitudes();
    let (sgnl_forward, sgnl_backward) = signal.amplitudes();
    let boundary = FieldState {
        sgnl_f: sgnl_forward,
        sgnl_b: 0.0,
        pump_f: pump_forward,
        pump_b: 0.0,
    };

    let set_boundary = |current: &[FieldState]| {
        let (sgnl_b, pump_b) = find_b_fields(sgnl_backward, pump_backward, current, fp, dz);
        FieldState {
            sgnl_b,
            pump_b,
            ..boundary
        }
    };

    let step = |new_previous: &FieldState, old_current: &FieldState, _i| {
        new_previous.uncoupled_step(fp.gain(*old_current), dz)
    };

    let error = |current: &[FieldState], previous: &[FieldState]| {
        profile_convergence_error(
            current,
            previous,
            config.absolute_tolerance,
            config.relative_tolerance,
        )
    };

    solver.solve(config.max_iterations, set_boundary, step, error)
}

pub fn initial_profile(
    signal: Signal,
    pump: Pump,
    fp: &ResolvedFibre<'_>,
    gp: GridPoints,
) -> FieldProfile {
    let gain = fp.initial_gain();
    let gpump = gain.pump;
    let gsignal = gain.signal;
    let zs = gp.grid(fp.length());
    let pump_end_factor = (0.5 * gpump * fp.length()).exp();
    let signal_end_factor = (0.5 * gsignal * fp.length()).exp();
    let (pump_forward, pump_backward) = pump.amplitudes();
    let (signal_forward, signal_backward) = signal.amplitudes();

    let fields = zs
        .iter()
        .map(|z| {
            let pump_factor = (0.5 * gpump * z).exp();
            let signal_factor = (0.5 * gsignal * z).exp();

            FieldState {
                sgnl_f: signal_factor * signal_forward,
                sgnl_b: signal_end_factor / signal_factor * signal_backward,
                pump_f: pump_factor * pump_forward,
                pump_b: pump_end_factor / pump_factor * pump_backward,
            }
        })
        .collect();
    FieldProfile::new(zs, fields)
}

pub fn solve_amp_picard(
    fibre: &ResolvedFibre<'_>,
    signal: Signal,
    pump: Pump,
    config: AmplifierSolveConfig,
    full_profile: bool,
) -> Result<FieldProfile, SolverError> {
    let gp = config.grid_points;
    let dz = gp.dz(fibre.length());
    let initial = initial_profile(signal, pump, fibre, gp);
    let mut solver = PicardSolver::from_initial(initial.fields);
    solve_amp_profile_picard(&mut solver, signal, pump, fibre, config.picard, dz)?;

    let fields = solver.profile();
    if full_profile {
        Ok(FieldProfile::new(gp.grid(fibre.length()), fields.to_vec()))
    } else {
        Ok(FieldProfile::new(
            vec![0.0, fibre.length()],
            vec![fields[0], fields.last().copied().unwrap()],
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::maths::rootfind::BisectionConfig;

    #[test]
    fn backward_signal_uses_picard_boundary_conditions() {
        let fibre = crate::lase::Fibre {
            dopant: crate::lase::TwoLevelDopant {
                density: 0.0,
                lifetime: 1.0,
                pump_ab: 1.0,
                pump_em: 0.0,
                sgnl_ab: 0.0,
                sgnl_em: 1.0,
            },
            ..crate::lase::Fibre::default()
        };
        let fibre = fibre.resolve(
            crate::lase::FieldMode::new(970e-9),
            crate::lase::FieldMode::new(1060e-9),
        );
        let signal = Signal {
            total: 2.0,
            balance: 0.0,
        };
        let pump = Pump {
            total: 2.0,
            balance: 0.0,
        };
        let config = AmplifierSolveConfig {
            grid_points: GridPoints(10),
            root_find: BisectionConfig::default().into(),
            picard: PicardConfig::default(),
        };

        let profile = Amplifier { fibre }
            .solve(signal, pump, config, true)
            .expect("zero-gain amplifier Picard solve should converge");

        let left = profile.fields.first().unwrap();
        let right = profile.fields.last().unwrap();
        assert_eq!(left.sgnl_f, 1.0);
        assert_eq!(left.pump_f, 1.0);
        assert_eq!(right.sgnl_b, 1.0);
        assert_eq!(right.pump_b, 1.0);
    }
}
