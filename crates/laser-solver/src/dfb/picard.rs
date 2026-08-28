use super::{DfbLaser, DfbSolveConfig};
use crate::error::SolverError;
use crate::lase::{FieldProfile, FieldState, OutputPower, Pump, ResolvedFibre};
use crate::picard::{PicardConfig, PicardError, PicardSolver};
use crate::rootfind::try_rootfind_1d;

pub fn find_pump_b(pump: Pump, profile: &[FieldState], fp: &ResolvedFibre<'_>, dz: f64) -> f64 {
    let expg: f64 = profile[..profile.len() - 1]
        .iter()
        .map(|&field| 0.5 * fp.gain(field).pump * dz)
        .sum::<f64>() // dont know why it couldnt infer f64 here
        .exp();
    pump.backward_amplitude() * expg
}

// 0. if a function sig borrows nothing but returns a borrow &T
// then it basically must be borrowing from static so 'a='static
// fn foo<T>(x: T, y: T) -> &'static str { "hello" }
//
// 1. If there is one borrowed argument the compiler will
// assume the returned borrow comes from that, so
// fn foo<T, V>(x: &T, y: T) -> &V   implicitly becomes
// fn foo<'a, T, V>(x: &'a T, y: T) -> &'a V
// but still possible to write
// fn foo<T, V>(x: &T, y: T) -> &'static V
//
// 2. If there is more than one borrowed argument (as below) then the returned
// borrow must be 'static or have come from one of the argument borrows
// but compiler cant know which from looking at the args/return type alone
// so must annotate it
// fn foo<'a, T, V>(x: &'a T, y: &T) -> &'a V
// or
// fn foo<T, V>(x: &T, y: &T) -> &'static V
//
// here we are returning a borrowed slice of the input borrowed PicardSolver
// but we also pass a borrowed slice of kappas, so must annotate
// to show that &[FieldState] comes from &mut PicardSolver
pub fn solve_profile_picard<'a>(
    solver: &'a mut PicardSolver,
    sgnl_b: f64,
    pump: Pump,
    fp: &ResolvedFibre<'_>,
    config: PicardConfig,
    kappas: &[f64],
    dz: f64,
) -> Result<&'a [FieldState], PicardError> {
    assert_eq!(kappas.len() + 1, solver.profile().len());
    let (pump_forward, _) = pump.amplitudes();
    let boundary = FieldState {
        sgnl_f: 0.0,
        sgnl_b,
        pump_f: pump_forward,
        pump_b: 0.0,
    };

    let set_boundary = |current: &[FieldState]| FieldState {
        pump_b: find_pump_b(pump, current, fp, dz),
        ..boundary
    };

    let step = |new_previous: FieldState, old_current: FieldState, i| {
        new_previous.coupled_step_general(old_current, fp, kappas[i], dz)
    };

    solver.solve(config, set_boundary, step)
}

impl DfbLaser<'_> {
    fn solve_with_picard_solver(
        &self,
        pump: Pump,
        config: DfbSolveConfig,
        full_profile: bool,
        solver: &mut PicardSolver,
    ) -> Result<FieldProfile, SolverError> {
        let gp = config.grid_points;
        let kappas = self.grating.grid(gp.0);
        let dz = gp.dz(self.fibre.length());
        let f = |sgnl_b| -> Result<f64, SolverError> {
            let fields = solve_profile_picard(
                solver,
                sgnl_b,
                pump,
                &self.fibre,
                config.picard,
                &kappas,
                dz,
            )?;
            Ok(fields.last().unwrap().sgnl_b / sgnl_b)
        };
        // try_rootfind_1d muts the solver which leaves the lasing solution in the 'current' buffer
        let _sgnl_b = try_rootfind_1d(f, config.root_find)?;
        if full_profile {
            Ok(FieldProfile::new(
                gp.grid(self.fibre.length()),
                solver.profile().to_vec(),
            ))
        } else {
            let fields = solver.profile();
            Ok(FieldProfile::new(
                vec![0.0_f64, self.fibre.length()],
                vec![fields[0], fields.last().copied().unwrap()],
            ))
        }
    }

    pub fn solve_picard(
        &self,
        pump: Pump,
        config: DfbSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let mut solver = self.initial_picard_solver(pump, config.grid_points);
        self.solve_with_picard_solver(pump, config, full_profile, &mut solver)
    }

    pub(crate) fn output_power_picard(
        &self,
        pump: Pump,
        config: DfbSolveConfig,
        solver: &mut PicardSolver,
    ) -> Result<OutputPower, SolverError> {
        let profile = self.solve_with_picard_solver(pump, config, false, solver)?;
        Ok(profile.output_powers())
    }
}
