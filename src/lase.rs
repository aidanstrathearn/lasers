use crate::error::SolverError;
use crate::maths::rootfind::RootFindError;

pub use crate::dopant::{
    DopantError, DopantModel, TwoLevelCrossSections, TwoLevelDopant, TwoLevelPopulations,
};
pub use crate::fibre::{BidirectionalAmplitude, Fibre, FibreGeometry, FieldMode};
pub use crate::two_mode::{
    FieldProfile, FieldState, Gain, OutputPower, Pump, ResolvedFibre, Signal, field_avg_diff,
    field_max_diff, profile_avg_diff, profile_convergence_error, profile_max_diff,
};

pub type PumpScan = Vec<Option<OutputPower>>;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct UniformGrid {
    length: f64,
    steps: usize,
    dz: f64,
}

impl UniformGrid {
    pub fn new(length: f64, steps: usize) -> Self {
        assert!(
            length.is_finite() && length > 0.0,
            "grid length must be positive and finite"
        );
        assert!(steps > 0, "grid must contain at least one step");
        Self {
            length,
            steps,
            dz: length / steps as f64,
        }
    }

    pub fn length(self) -> f64 {
        self.length
    }

    pub fn steps(self) -> usize {
        self.steps
    }

    pub fn points(self) -> usize {
        self.steps + 1
    }

    pub fn dz(self) -> f64 {
        self.dz
    }

    pub fn position(self, point: usize) -> f64 {
        assert!(point <= self.steps, "grid point is out of bounds");
        if point == self.steps {
            self.length
        } else {
            point as f64 * self.dz
        }
    }

    pub fn positions(self) -> impl ExactSizeIterator<Item = f64> {
        (0..self.points()).map(move |point| self.position(point))
    }
}

pub fn classify_output(
    result: Result<OutputPower, SolverError>,
) -> Result<Option<OutputPower>, SolverError> {
    match result {
        Ok(output) => Ok(Some(output)),
        // not bracketed error likely means below threshold, not true error
        Err(SolverError::RootFind(RootFindError::RootNotBracketed)) => Ok(None),
        Err(error) => Err(error),
    }
}

pub fn pump_scan(
    pumps: &[f64],
    mut output_power: impl FnMut(f64) -> Result<OutputPower, SolverError>,
) -> Result<PumpScan, SolverError> {
    pumps
        .iter()
        .map(|&pump| classify_output(output_power(pump)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::maths::picard::PicardError;
    use crate::fibre::transfer;

    #[test]
    fn uniform_grid_distinguishes_steps_from_points() {
        let grid = UniformGrid::new(10.0, 4);

        assert_eq!(grid.length(), 10.0);
        assert_eq!(grid.steps(), 4);
        assert_eq!(grid.points(), 5);
        assert_eq!(grid.dz(), 2.5);
        assert_eq!(
            grid.positions().collect::<Vec<_>>(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0]
        );
    }

    #[test]
    fn check_transfer() {
        let (a, b, c, d) = transfer(1.0, 0.0, 1.0);
        println!("Transfer {:?}", (a, b, c, d));
        assert_eq!(a, (0.5_f64).exp());
    }

    #[test]
    fn pump_scan_marks_only_unbracketed_roots_as_below_threshold() {
        let samples = pump_scan(&[0.0, 1.0], |pump| {
            if pump == 0.0 {
                Err(RootFindError::RootNotBracketed.into())
            } else {
                Ok((2.0, 3.0))
            }
        })
        .unwrap();

        assert_eq!(samples, vec![None, Some((2.0, 3.0))]);
    }

    #[test]
    fn pump_scan_propagates_picard_errors() {
        let result = pump_scan(&[0.0], |_| {
            Err(SolverError::Picard(PicardError::DidNotConverge))
        });

        assert!(matches!(
            result,
            Err(SolverError::Picard(PicardError::DidNotConverge))
        ));
    }
}
