use super::{DfbLaser, DfbSolveConfig};
use crate::dopant::DopantModel;
use crate::error::SolverError;
use crate::lase::{
    BidirectionalAmplitude, FieldProfile, FieldState, OutputPower, Pump, UniformGrid,
};
use crate::maths::rootfind::rootfind_1d;
use crate::propagation::{out_field_coupled, solve_profile_coupled};

impl<D: DopantModel> DfbLaser<'_, D> {
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
        let grid = UniformGrid::new(self.fibre.length(), config.steps);
        let kappas = self.grating.grid(grid.steps());
        let dz = grid.dz();
        let trial = |sgnl_b| FieldState {
            signal: BidirectionalAmplitude {
                forward: 0.0,
                backward: sgnl_b,
            },
            pump: BidirectionalAmplitude {
                forward: pump_forward,
                backward: 0.0, // shooting method requires zero backward pump amplitude
            },
        };
        let f = |sgnl_b| {
            out_field_coupled(trial(sgnl_b), |fields| self.fibre.gain(fields), dz, &kappas)
                .signal
                .backward
                / sgnl_b
        };
        let sgnl_b = rootfind_1d(f, config.root_find)?;

        if full_profile {
            let z = grid.positions().collect();
            let fields =
                solve_profile_coupled(trial(sgnl_b), |fields| self.fibre.gain(fields), dz, &kappas);
            Ok(FieldProfile::new(z, fields))
        } else {
            let z = vec![0.0_f64, self.fibre.length()];
            let out_left = trial(sgnl_b);
            let fields = vec![
                out_left,
                out_field_coupled(out_left, |fields| self.fibre.gain(fields), dz, &kappas),
            ];
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
