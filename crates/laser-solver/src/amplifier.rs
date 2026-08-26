use crate::error::SolverError;
use crate::lase::{Fibre, FieldProfile, FieldState, GridPoints, OutputPower, Pump};
use crate::rootfind::{RootFindConfig, rootfind_1d};

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
        input_signal_power: f64,
        pump: Pump,
        config: AmplifierSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        assert!(
            input_signal_power >= 0.0,
            "input signal power must be greater than or equal to zero"
        );
        let (pump_forward, pump_backward) = pump.amplitudes();
        let gp = config.grid_points;
        let nsteps = gp.0;
        let dz = gp.dz(self.fibre.length);

        let trial = |pump_b| FieldState {
            sgnl_f: input_signal_power.sqrt(),
            sgnl_b: 0.0,
            pump_f: pump_forward,
            pump_b,
        };
        let pump_b = if pump_backward == 0.0 {
            0.0
        } else {
            let f =
                |pump_b| out_field_uncoupled(trial(pump_b), self.fibre, dz, nsteps).pump_b / pump_backward - 1.0;
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
        input_signal_power: f64,
        pump: Pump,
        config: AmplifierSolveConfig,
    ) -> Result<OutputPower, SolverError> {
        let profile = self.solve(input_signal_power, pump, config, false)?;
        Ok(profile.output_powers())
    }
}

pub fn solve_profile_uncoupled(
    fs: FieldState,
    fp: Fibre,
    dz: f64,
    nsteps: usize,
) -> Vec<FieldState> {
    let mut current = fs;
    let mut result = Vec::with_capacity(nsteps + 1);
    result.push(current);
    for _ in 0..nsteps {
        current = current.uncoupled_step_shooting(fp, dz);
        result.push(current);
    }
    result
}

pub fn out_field_uncoupled(fs: FieldState, fp: Fibre, dz: f64, nsteps: usize) -> FieldState {
    let mut current = fs;
    for _ in 0..nsteps {
        current = current.uncoupled_step_shooting(fp, dz);
    }
    current
}
