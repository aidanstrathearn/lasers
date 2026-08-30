use crate::error::SolverError;
use crate::grating::GratingModel;
use crate::lase::{DopantModel, FieldState};
use crate::maths::picard::PicardConfig;
use crate::two_mode::{FieldProfile, Pump, ResolvedFibre, Signal};

pub struct TwoModeSolver<'a, D: DopantModel, G: GratingModel> {
    fibre: &'a ResolvedFibre<'a, D, G>,
}

impl<D: DopantModel, G: GratingModel> TwoModeSolver<'_, D, G> {
    pub fn solve_injected(&self, pump: Pump, signal: Signal, picard_config: PicardConfig) -> Result<FieldProfile, SolverError> {
        let use_picard = pump.balance == 1.0 && signal.balance == 1.0;
        let solution = if use_picard {
            vec![FieldState::default(); self.fibre.grid.points()]
        } else {
            vec![FieldState::default(); self.fibre.grid.points()]
        };
        
        Result::Ok::<FieldProfile, SolverError>(FieldProfile::new(
            self.fibre.grid.positions().collect(),
            solution,
        ))
    }

    pub fn solve_lasing() -> Result<FieldProfile, SolverError> {
        Result::Ok::<FieldProfile, SolverError>(FieldProfile::new(
            vec![0.0],
            vec![FieldState::default()],
        ))
    }
}
