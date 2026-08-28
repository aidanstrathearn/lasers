use crate::error::SolverError;
use crate::rootfind::RootFindError;
use crate::utils::{IterationConfig, linspace, relative_diff};

pub type OutputPower = (f64, f64);
pub type PumpScan = Vec<Option<OutputPower>>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gain {
    pub pump: f64,
    pub signal: f64,
}

#[derive(Clone, Copy)]
pub struct TwoLevelDopant {
    pub density: f64,
    pub lifetime: f64,
    pub pump_ab: f64,
    pub pump_em: f64,
    pub sgnl_ab: f64,
    pub sgnl_em: f64,
}

impl TwoLevelDopant {
    pub fn steady_state(&self, gamma_dn: f64, gamma_up: f64) -> (f64, f64) {
        let gamma_decay = 1.0 / self.lifetime;
        let gamma_dn_total = gamma_decay + gamma_dn;
        let denom = gamma_up + gamma_dn_total;
        (gamma_dn_total / denom, gamma_up / denom)
    }
    pub fn pops(&self, fs: FieldState, pump_overlap: f64, sgnl_overlap: f64) -> (f64, f64) {
        let pump_flux = (fs.pump_f * fs.pump_f + fs.pump_b * fs.pump_b) * pump_overlap;

        let sgnl_flux = (fs.sgnl_f * fs.sgnl_f + fs.sgnl_b * fs.sgnl_b) * sgnl_overlap;

        let gamma_up = pump_flux * self.pump_ab + sgnl_flux * self.sgnl_ab;
        let gamma_dn = pump_flux * self.pump_em + sgnl_flux * self.sgnl_em;
        self.steady_state(gamma_dn, gamma_up)
    }

    pub fn gain(&self, fs: FieldState, pump_overlap: f64, sgnl_overlap: f64) -> Gain {
        let (g, e) = self.pops(fs, pump_overlap, sgnl_overlap);
        self.gain_from_populations(g, e, pump_overlap, sgnl_overlap)
    }

    fn gain_from_populations(
        &self,
        ground: f64,
        excited: f64,
        pump_overlap: f64,
        sgnl_overlap: f64,
    ) -> Gain {
        Gain {
            pump: self.density * (-ground * self.pump_ab + excited * self.pump_em) * pump_overlap,
            signal: self.density * (-ground * self.sgnl_ab + excited * self.sgnl_em) * sgnl_overlap,
        }
    }
}

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
const SPEED_OF_LIGHT_MS: f64 = 299_792_458.0;

fn numerical_aperture(n_core: f64, n_cladding: f64) -> f64 {
    (n_core * n_core - n_cladding * n_cladding).sqrt()
}

fn v_number(numerical_aperture: f64, core_radius: f64, wavelength: f64) -> f64 {
    TWO_PI * core_radius * numerical_aperture / wavelength
}

fn dimensionless_marcuse_radius(v_number: f64) -> f64 {
    // this is an approximation for the ratio mode_radius / core_radius
    // D. Marcuse, “Loss analysis of single-mode fiber splices,” Bell System Technical Journal, vol. 56, no. 5, pp. 703–718, 1977.
    0.65 + 1.619 / v_number.powf(1.5) + 2.879 / v_number.powi(6)
}

#[derive(Copy, Clone)]
pub struct FieldMode {
    wavelength: f64,
}

impl FieldMode {
    pub const fn new(wavelength: f64) -> Self {
        Self { wavelength }
    }

    pub fn wavelength(self) -> f64 {
        self.wavelength
    }
}

impl Default for FieldMode {
    fn default() -> Self {
        Self::new(1.0)
    }
}

#[derive(Clone)]
pub struct FibreGeometry {
    pub core_radius: f64,
    pub numerical_aperture: f64,
    pub length: f64,
}

impl FibreGeometry {
    fn single_mode_cutoff_wavelength(&self) -> f64 {
        TWO_PI * self.core_radius * self.numerical_aperture / 2.405
    }

    fn is_single_spatial_mode(&self, mode: FieldMode) -> bool {
        mode.wavelength > self.single_mode_cutoff_wavelength()
    }

    fn v_number(&self, mode: FieldMode) -> f64 {
        TWO_PI * self.core_radius * self.numerical_aperture / mode.wavelength
    }

    fn mode_overlap(&self, mode: FieldMode) -> f64 {
        //let _ = mode;
        let v = self.v_number(mode);
        let mode_over_core = dimensionless_marcuse_radius(v);
        let gamma = 1.0 - f64::exp(-2.0 / (mode_over_core * mode_over_core));
        gamma * 1e-10 + 1.0
        //1.0
    }
}

#[derive(Clone)]
pub struct Fibre {
    pub geometry: FibreGeometry,
    pub dopant: TwoLevelDopant,
}

impl Fibre {
    pub fn resolve(&self, pump_mode: FieldMode, sgnl_mode: FieldMode) -> ResolvedFibre<'_> {
        ResolvedFibre {
            fibre: self,
            pump_mode,
            sgnl_mode,
            pump_overlap: self.geometry.mode_overlap(pump_mode),
            sgnl_overlap: self.geometry.mode_overlap(sgnl_mode),
        }
    }
}

#[derive(Clone)]
pub struct ResolvedFibre<'a> {
    fibre: &'a Fibre,
    pump_mode: FieldMode,
    sgnl_mode: FieldMode,
    pump_overlap: f64,
    sgnl_overlap: f64,
}

impl ResolvedFibre<'_> {
    pub fn length(&self) -> f64 {
        self.fibre.geometry.length
    }

    pub fn pump_mode(&self) -> FieldMode {
        self.pump_mode
    }

    pub fn sgnl_mode(&self) -> FieldMode {
        self.sgnl_mode
    }

    pub fn pump_overlap(&self) -> f64 {
        self.pump_overlap
    }

    pub fn sgnl_overlap(&self) -> f64 {
        self.sgnl_overlap
    }

    pub fn gain(&self, fs: FieldState) -> Gain {
        self.fibre
            .dopant
            .gain(fs, self.pump_overlap, self.sgnl_overlap)
    }

    pub fn populations(&self, fs: FieldState) -> (f64, f64) {
        self.fibre
            .dopant
            .pops(fs, self.pump_overlap, self.sgnl_overlap)
    }

    pub fn initial_gain(&self) -> Gain {
        self.fibre
            .dopant
            .gain_from_populations(0.5, 0.5, self.pump_overlap, self.sgnl_overlap)
    }
}

impl Default for Fibre {
    fn default() -> Self {
        Self {
            geometry: FibreGeometry {
                core_radius: 1.0,
                numerical_aperture: 0.1,
                length: 5.0,
            },
            dopant: TwoLevelDopant {
                density: 1.0,
                lifetime: 1.0,
                pump_ab: 0.01,
                pump_em: 0.0,
                sgnl_ab: 0.0,
                sgnl_em: 1.0,
            },
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
            "signal total must be non-negative and balance must be between -1 and 1"
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
    use crate::picard::PicardError;
    use crate::propagation::transfer;

    #[test]
    fn check_transfer() {
        let (a, b, c, d) = transfer(1.0, 0.0, 1.0);
        println!("Transfer {:?}", (a, b, c, d));
        assert_eq!(a, (0.5_f64).exp());
    }

    #[test]
    fn two_level_dopant_returns_named_modal_gain() {
        let dopant = TwoLevelDopant {
            density: 2.0,
            lifetime: 1.0,
            pump_ab: 3.0,
            pump_em: 0.0,
            sgnl_ab: 5.0,
            sgnl_em: 0.0,
        };

        let gain = dopant.gain(FieldState::default(), 0.25, 0.5);

        assert_eq!(
            gain,
            Gain {
                pump: -1.5,
                signal: -5.0,
            }
        );
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
