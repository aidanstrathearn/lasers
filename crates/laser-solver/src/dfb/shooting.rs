use super::{DfbLaser, DfbSolveConfig, out_field, solve_profile};
use crate::error::SolverError;
use crate::lase::{FieldProfile, FieldState, OutputPower, Pump};
use crate::rootfind::rootfind_1d;

impl DfbLaser {
    pub fn solve_shooting(
        &self,
        pump: Pump,
        config: DfbSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let (pump_forward, pump_backward) = pump.amplitudes();
        assert_eq!(
            pump_backward, 0.0,
            "shooting solver requires a forward-only pump"
        );
        let gp = config.grid_points;
        let kappas = self.grating.grid(gp.0);
        let dz = gp.dz(self.fibre.length);
        let trial = |sgnl_b| FieldState {
            sgnl_f: 0.0,
            sgnl_b,
            pump_f: pump_forward,
            pump_b: 0.0, // shooting method requires zero backward pump amplitude
        };
        let f = |sgnl_b| out_field(trial(sgnl_b), self.fibre, dz, &kappas).sgnl_b / sgnl_b;
        let sgnl_b = rootfind_1d(f, config.root_find)?;

        if full_profile {
            let z = gp.grid(self.fibre.length);
            let fields = solve_profile(trial(sgnl_b), self.fibre, dz, &kappas);
            Ok(FieldProfile::new(z, fields))
        } else {
            let z = vec![0.0_f64, self.fibre.length];
            let out_left = trial(sgnl_b);
            let fields = vec![out_left, out_field(out_left, self.fibre, dz, &kappas)];
            Ok(FieldProfile::new(z, fields))
        }
    }

    pub(crate) fn output_power_shooting(
        &self,
        pump: Pump,
        config: DfbSolveConfig,
    ) -> Result<OutputPower, SolverError> {
        let profile = self.solve_shooting(pump, config, false)?;
        Ok(profile.output_powers())
    }
}
