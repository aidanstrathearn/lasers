use crate::dfb::picard::find_pump_b;
use crate::error::SolverError;
use crate::lase::{gain, Fibre, FieldProfile, FieldState, GridPoints, OutputPower, Pump};
use crate::picard::{PicardConfig, PicardError, PicardSolver};
use crate::propagation::{out_field_uncoupled, solve_profile_uncoupled};
use crate::rootfind::{RootFindConfig, rootfind_1d};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Signal {
    pub total: f64,
    pub balance: f64,
}

impl Default for Signal {
    fn default() -> Self {
        Self {
            total: 1.0,
            balance: 1.0,
        }
    }
}

impl Signal {
    pub fn amplitudes(self) -> (f64, f64) {
        assert!(
            self.total >= 0.0 && (-1.0..=1.0).contains(&self.balance),
            "pump total must be non-negative and balance must be between -1 and 1"
        );
        let forward_fraction = (self.balance + 1.0) * 0.5;
        (
            (forward_fraction * self.total).sqrt(),
            ((1.0 - forward_fraction) * self.total).sqrt(),
        )
    }

    pub fn forward_amplitude(self) -> f64 {
        self.amplitudes().0
    }

    pub fn backward_amplitude(self) -> f64 {
        self.amplitudes().1
    }
}

#[derive(Copy, Clone)]
pub struct AmplifierSolveConfig {
    pub grid_points: GridPoints,
    pub root_find: RootFindConfig,
}

#[derive(Copy, Clone)]
pub struct Amplifier {
    pub fibre: Fibre,
}

impl Amplifier {
    pub fn solve(
        &self,
        signal: Signal,
        pump: Pump,
        config: AmplifierSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let forward_signal = signal.forward_amplitude();
        assert!(
            forward_signal >= 0.0,
            "input signal power must be greater than or equal to zero"
        );
        let (pump_forward, pump_backward) = pump.amplitudes();
        let gp = config.grid_points;
        let nsteps = gp.0;
        let dz = gp.dz(self.fibre.length);

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
                out_field_uncoupled(trial(pump_b), self.fibre, dz, nsteps).pump_b / pump_backward
                    - 1.0
            };
            rootfind_1d(f, config.root_find)?
        };

        if full_profile {
            let z = gp.grid(self.fibre.length);
            let fields = solve_profile_uncoupled(trial(pump_b), self.fibre, dz, nsteps);
            Ok(FieldProfile::new(z, fields))
        } else {
            let z = vec![0.0_f64, self.fibre.length];
            let out_left = trial(pump_b);
            let fields = vec![
                out_left,
                out_field_uncoupled(out_left, self.fibre, dz, nsteps),
            ];
            Ok(FieldProfile::new(z, fields))
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

pub fn find_b_fields(signal_b_right: f64,
                     pump_b_right: f64, profile: &[FieldState], fp: Fibre, dz: f64) -> (f64, f64) {
    let (pump_od, signal_od): (f64, f64) = profile[..profile.len() - 1]
        .iter()
        .fold((0.0, 0.0), |acc, &field| {
            let (gp, gs) = gain(field, fp);
            let new = (0.5 * gp * dz, 0.5 * gs * dz);
            (acc.0 + new.0, acc.1 + new.1)
        });

    (
        signal_b_right * signal_od.exp(),
        pump_b_right * pump_od.exp(),
    )
}

pub fn solve_amp_profile_picard<'a>(
    solver: &'a mut PicardSolver,
    signal: Signal,
    pump: Pump,
    fp: Fibre,
    config: PicardConfig,
    kappas: &[f64],
    dz: f64,
) -> Result<&'a [FieldState], PicardError> {
    assert_eq!(kappas.len() + 1, solver.profile().len());
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



    let step = |new_previous: FieldState, old_current: FieldState, _i| {
        new_previous.uncoupled_step_general(old_current, fp, dz)
    };

    solver.solve(
        config,
        set_boundary,
        step,
    )
}
