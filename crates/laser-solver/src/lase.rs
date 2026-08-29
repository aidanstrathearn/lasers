use crate::error::SolverError;
use crate::maths::rootfind::RootFindError;
use crate::maths::utils::{IterationConfig, relative_diff};

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

pub fn find_threshold_and_slope(
    pump_start: f64,
    pump_step: f64,
    ip: IterationConfig,
    mut output_power: impl FnMut(f64) -> Result<OutputPower, SolverError>,
) -> Result<(f64, f64, f64), SolverError> {
    assert!(pump_step > 0.0, "pump step must be positive");

    let mut current_pump = pump_start;
    let mut previous_output = None;
    let mut previous_total_slope = None;

    for _ in 0..ip.max {
        if let Some((new_sf, new_sb)) = classify_output(output_power(current_pump))? {
            if let Some((previous_pump, sf, sb)) = previous_output {
                let dp = current_pump - previous_pump;
                let slope_f = (new_sf - sf) / dp;
                let slope_b = (new_sb - sb) / dp;
                let total_slope = slope_f + slope_b;

                if previous_total_slope.is_some_and(|previous| {
                    relative_diff(total_slope, previous) < ip.tol && total_slope > 0.0
                }) {
                    let threshold = current_pump - (new_sf + new_sb) / total_slope;
                    return Ok((slope_f, slope_b, threshold));
                }

                previous_total_slope = Some(total_slope);
            }

            previous_output = Some((current_pump, new_sf, new_sb));
        }

        current_pump += pump_step;
    }

    Err(SolverError::ThresholdNotFound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::maths::picard::PicardError;
    use crate::propagation::transfer;

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

    #[test]
    fn threshold_slope_uses_spacing_between_successful_samples() {
        let result =
            find_threshold_and_slope(0.0, 1.0, IterationConfig { max: 6, tol: 1e-12 }, |pump| {
                match pump as usize {
                    0 | 2 | 4 => Err(RootFindError::RootNotBracketed.into()),
                    _ => Ok((2.0 * (pump - 1.0), 3.0 * (pump - 1.0))),
                }
            })
            .unwrap();

        assert_eq!(result.0, 2.0);
        assert_eq!(result.1, 3.0);
        assert!((result.2 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn threshold_slope_propagates_picard_errors() {
        let result =
            find_threshold_and_slope(0.0, 1.0, IterationConfig { max: 2, tol: 1e-3 }, |_| {
                Err(SolverError::Picard(PicardError::DidNotConverge))
            });

        assert!(matches!(
            result,
            Err(SolverError::Picard(PicardError::DidNotConverge))
        ));
    }
}
