use crate::error::SolverError;
use crate::rootfind::RootFindError;
use crate::utils::{IterationConfig, linspace, relative_diff};

pub type OutputPower = (f64, f64);
pub type PumpScan = Vec<Option<OutputPower>>;

#[derive(Copy, Clone)]
pub struct FibreParams {
    pub density: f64,
    pub lifetime: f64,
    pub pump_ab: f64,
    pub pump_em: f64,
    pub sgnl_ab: f64,
    pub sgnl_em: f64,
    pub length: f64,
}

impl Default for FibreParams {
    fn default() -> Self {
        Self {
            density: 1.0,
            lifetime: 1.0,
            pump_ab: 0.01,
            pump_em: 0.0,
            sgnl_ab: 0.0,
            sgnl_em: 1.0,
            length: 5.0,
        }
    }
}

#[allow(dead_code)]
#[derive(Copy, Clone)]
pub struct Reflectivities {
    left: f64,
    right: f64,
}

#[derive(Copy, Clone)]
pub struct GridPoints(pub usize);

impl Default for GridPoints {
    fn default() -> Self {
        Self(100)
    }
}

impl GridPoints {
    pub fn grid(self, length: f64) -> Vec<f64> {
        linspace(0.0, length, self.0)
    }

    pub fn dz(self, length: f64) -> f64 {
        length / self.0 as f64
    }
}

#[derive(Copy, Clone)]
pub struct GratingProfile {
    pub kappa_left: f64,
    pub kappa_right: f64,
    pub pi_shift_position: f64,
}

impl Default for GratingProfile {
    fn default() -> Self {
        Self {
            kappa_left: 1.0,
            kappa_right: 1.0,
            pi_shift_position: 0.45,
        }
    }
}

impl GratingProfile {
    pub fn grid(self, n: usize) -> Vec<f64> {
        (0..n)
            .map(|j| {
                let z = j as f64 / n as f64;
                if z < self.pi_shift_position {
                    self.kappa_left
                } else {
                    -self.kappa_right
                }
            })
            .collect()
    }
}

#[derive(Copy, Clone, Debug, Default)]
pub struct FieldState {
    pub sgnl_f: f64,
    pub sgnl_b: f64,
    pub pump_f: f64,
    pub pump_b: f64,
}

impl FieldState {
    pub(crate) fn field_powers(self) -> [f64; 2] {
        [
            self.sgnl_f * self.sgnl_f + self.sgnl_b * self.sgnl_b,
            self.pump_f * self.pump_f + self.pump_b * self.pump_b,
        ]
    }
}

pub fn field_max_diff(f1: FieldState, f2: FieldState) -> f64 {
    let diffs = [
        relative_diff(f1.pump_f, f2.pump_f),
        relative_diff(f1.pump_b, f2.pump_b),
        relative_diff(f1.sgnl_f, f2.sgnl_f),
        relative_diff(f1.sgnl_b, f2.sgnl_b),
    ];
    diffs.into_iter().fold(0.0, f64::max)
}
pub fn profile_max_diff(p1: &[FieldState], p2: &[FieldState]) -> f64 {
    p1.iter()
        .zip(p2.iter())
        .map(|(&f1, &f2)| field_max_diff(f1, f2))
        .reduce(f64::max)
        .unwrap_or(f64::NAN)
}

pub fn field_avg_diff(f1: FieldState, f2: FieldState) -> f64 {
    let diffs = [
        relative_diff(f1.pump_f, f2.pump_f),
        relative_diff(f1.pump_b, f2.pump_b),
        relative_diff(f1.sgnl_f, f2.sgnl_f),
        relative_diff(f1.sgnl_b, f2.sgnl_b),
    ];
    diffs.into_iter().sum::<f64>() / 4.0
}
pub fn profile_avg_diff(p1: &[FieldState], p2: &[FieldState]) -> f64 {
    p1.iter()
        .zip(p2.iter())
        .map(|(&f1, &f2)| field_max_diff(f1, f2))
        .sum::<f64>()
        / p1.len() as f64
}

#[derive(Clone)]
pub struct FieldProfile {
    pub z: Vec<f64>,
    pub fields: Vec<FieldState>,
}

impl FieldProfile {
    pub fn new(z: Vec<f64>, fields: Vec<FieldState>) -> Self {
        assert_eq!(z.len(), fields.len());
        Self { z, fields }
    }

    pub fn sgnl_f(&self) -> impl Iterator<Item = f64> + '_ {
        // borrows from self so '_ lifetime needs to match self
        // but looks like rust can infer this so doesnt need to be explicit
        self.fields.iter().map(|x| x.sgnl_f)
    }

    pub fn sgnl_b(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.sgnl_b)
    }

    pub fn pump_f(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.pump_f)
    }

    pub fn pump_b(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.pump_b)
    }

    pub fn z(&self) -> impl Iterator<Item = f64> {
        //manual deref needed here and not for self.fields because struct field access e.g. x.sgnl_b derefs implicitly.
        //could also do self.z.iter().copied()
        self.z.iter().copied()
    }

    pub fn output_powers(&self) -> OutputPower {
        let left = self.fields.first().expect("field profile is empty");
        let right = self.fields.last().expect("field profile is empty");
        (right.sgnl_f.powi(2), left.sgnl_b.powi(2))
    }
}

pub fn pops(fs: FieldState, fp: FibreParams) -> (f64, f64) {
    let pump_flux = fs.pump_f * fs.pump_f + fs.pump_b * fs.pump_b;
    let sgnl_flux = fs.sgnl_f * fs.sgnl_f + fs.sgnl_b * fs.sgnl_b;
    let gamma_up = pump_flux * fp.pump_ab + sgnl_flux * fp.sgnl_ab;
    let gamma_dn = pump_flux * fp.pump_em + sgnl_flux * fp.sgnl_em + 1.0 / fp.lifetime;
    let denom = gamma_up + gamma_dn;
    (gamma_dn / denom, gamma_up / denom)
}

pub fn gain(fs: FieldState, fp: FibreParams) -> (f64, f64) {
    let (g, e) = pops(fs, fp);
    (
        fp.density * (-g * fp.pump_ab + e * fp.pump_em),
        fp.density * (-g * fp.sgnl_ab + e * fp.sgnl_em),
    )
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pump {
    pub total: f64,
    pub balance: f64,
}

impl Default for Pump {
    fn default() -> Self {
        Self {
            total: 100.0,
            balance: 1.0,
        }
    }
}

impl Pump {
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
    use crate::dfb::transfer;
    use crate::picard::PicardError;

    #[test]
    fn check_transfer() {
        let (a, b, c, d) = transfer(1.0, 0.0, 1.0);
        println!("Transfer {:?}", (a, b, c, d));
        assert_eq!(a, (0.5_f64).exp());
    }

    #[test]
    fn pump_converts_power_and_balance_to_amplitudes() {
        for (pump, expected) in [
            (
                Pump {
                    total: 100.0,
                    balance: 1.0,
                },
                (10.0, 0.0),
            ),
            (
                Pump {
                    total: 100.0,
                    balance: -1.0,
                },
                (0.0, 10.0),
            ),
            (
                Pump {
                    total: 100.0,
                    balance: 0.0,
                },
                (50.0_f64.sqrt(), 50.0_f64.sqrt()),
            ),
            (
                Pump {
                    total: 0.0,
                    balance: 0.25,
                },
                (0.0, 0.0),
            ),
        ] {
            assert_eq!(pump.amplitudes(), expected);
        }
    }

    #[test]
    #[should_panic(expected = "pump total must be non-negative")]
    fn pump_rejects_negative_total_power() {
        Pump {
            total: -1.0,
            balance: 0.0,
        }
        .amplitudes();
    }

    #[test]
    #[should_panic(expected = "balance must be between -1 and 1")]
    fn pump_rejects_out_of_range_balance() {
        Pump {
            total: 1.0,
            balance: 2.0,
        }
        .amplitudes();
    }

    #[test]
    fn field_profile_reports_boundary_output_powers() {
        let profile = FieldProfile::new(
            vec![0.0, 1.0],
            vec![
                FieldState {
                    sgnl_b: -3.0,
                    ..FieldState::default()
                },
                FieldState {
                    sgnl_f: 2.0,
                    ..FieldState::default()
                },
            ],
        );

        assert_eq!(profile.output_powers(), (4.0, 9.0));
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
