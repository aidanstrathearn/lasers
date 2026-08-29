use crate::dopant::{TwoLevelDopant, TwoLevelPopulations};
use crate::fibre::{BidirectionalAmplitude, Fibre, FieldMode, bidirectional_amplitudes};
use crate::maths::utils::relative_diff;

pub type OutputPower = (f64, f64);

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gain {
    pub pump: f64,
    pub signal: f64,
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

    pub fn mode_fluxes(&self, fs: FieldState) -> (f64, f64) {
        (
            fs.pump.total_power() * self.pump_overlap,
            fs.signal.total_power() * self.sgnl_overlap,
        )
    }

    pub fn gain(&self, fs: FieldState) -> Gain {
        let (pump_flux, sgnl_flux) = self.mode_fluxes(fs);
        let mut gain = self.fibre.dopant.gain(pump_flux, sgnl_flux);

        gain.pump = gain.pump * self.pump_overlap;
        gain.signal = gain.signal * self.sgnl_overlap;
        gain
    }

    pub fn populations(&self, fs: FieldState) -> (f64, f64) {
        let (pump_flux, sgnl_flux) = self.mode_fluxes(fs);
        let data = &[
            (pump_flux, self.fibre.dopant.pump_cross_section()),
            (sgnl_flux, self.fibre.dopant.signal_cross_section()),
        ];
        let pops = self.fibre.dopant.pops(data);
        (pops.ground, pops.excited)

    }

    pub fn initial_gain(&self) -> Gain {
        Gain {pump: 0.0, signal: 0.0}
    }
}

#[derive(Copy, Clone, Debug, Default)]
pub struct FieldState {
    pub signal: BidirectionalAmplitude,
    pub pump: BidirectionalAmplitude,
}

impl FieldState {
    pub(crate) fn field_powers(self) -> [f64; 2] {
        [self.signal.total_power(), self.pump.total_power()]
    }
}

pub fn profile_convergence_error(
    current: &[FieldState],
    previous: &[FieldState],
    absolute_tolerance: f64,
    relative_tolerance: f64,
) -> f64 {
    assert_eq!(current.len(), previous.len());
    let mut max_dif_s = 0.0_f64;
    let mut max_dif_p = 0.0_f64;
    let mut max_mag_s = 0.0_f64;
    let mut max_mag_p = 0.0_f64;

    for (&current, &previous) in current.iter().zip(previous) {
        let current_powers = current.field_powers();
        let previous_powers = previous.field_powers();
        if !current_powers[0].is_finite()
            || !previous_powers[0].is_finite()
            || !current_powers[1].is_finite()
            || !previous_powers[1].is_finite()
        {
            return f64::INFINITY;
        }
        max_dif_s = max_dif_s.max((current_powers[0] - previous_powers[0]).abs().sqrt());
        max_dif_p = max_dif_p.max((current_powers[1] - previous_powers[1]).abs().sqrt());
        max_mag_s = max_mag_s.max(current_powers[0].max(previous_powers[0]).sqrt());
        max_mag_p = max_mag_p.max(current_powers[1].max(previous_powers[1]).sqrt());
    }
    (max_dif_p / (absolute_tolerance + relative_tolerance * max_mag_p))
        .max(max_dif_s / (absolute_tolerance + relative_tolerance * max_mag_s))
}

pub fn field_max_diff(f1: FieldState, f2: FieldState) -> f64 {
    let diffs = [
        relative_diff(f1.pump.forward, f2.pump.forward),
        relative_diff(f1.pump.backward, f2.pump.backward),
        relative_diff(f1.signal.forward, f2.signal.forward),
        relative_diff(f1.signal.backward, f2.signal.backward),
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
        relative_diff(f1.pump.forward, f2.pump.forward),
        relative_diff(f1.pump.backward, f2.pump.backward),
        relative_diff(f1.signal.forward, f2.signal.forward),
        relative_diff(f1.signal.backward, f2.signal.backward),
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
        self.fields.iter().map(|x| x.signal.forward)
    }

    pub fn sgnl_b(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.signal.backward)
    }

    pub fn pump_f(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.pump.forward)
    }

    pub fn pump_b(&self) -> impl Iterator<Item = f64> {
        self.fields.iter().map(|x| x.pump.backward)
    }

    pub fn z(&self) -> impl Iterator<Item = f64> {
        //manual deref needed here and not for self.fields because struct field access e.g. x.sgnl_b derefs implicitly.
        //could also do self.z.iter().copied()
        self.z.iter().copied()
    }

    pub fn output_powers(&self) -> OutputPower {
        let left = self.fields.first().expect("field profile is empty");
        let right = self.fields.last().expect("field profile is empty");
        (
            right.signal.forward_power(),
            left.signal.backward_power(),
        )
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
        bidirectional_amplitudes(self.total, self.balance)
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
        bidirectional_amplitudes(self.total, self.balance)
    }

    pub fn forward_amplitude(self) -> f64 {
        self.amplitudes().0
    }

    pub fn backward_amplitude(self) -> f64 {
        self.amplitudes().1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CONVERGENCE_RELATIVE_TOLERANCE: f64 = 1e-6;
    const CONVERGENCE_ABSOLUTE_TOLERANCE: f64 = 1e-10;

    fn convergence_error(current: &[FieldState], previous: &[FieldState]) -> f64 {
        profile_convergence_error(
            current,
            previous,
            CONVERGENCE_ABSOLUTE_TOLERANCE,
            CONVERGENCE_RELATIVE_TOLERANCE,
        )
    }

    #[test]
    fn identical_profiles_have_zero_convergence_error() {
        let profile = vec![FieldState {
            signal: BidirectionalAmplitude {
                forward: 1.0,
                backward: -1.0,
            },
            pump: BidirectionalAmplitude {
                forward: 100.0,
                backward: 10.0,
            },
        }];

        assert_eq!(convergence_error(&profile, &profile), 0.0);
    }

    #[test]
    fn tiny_zero_crossing_converges() {
        let current = vec![FieldState {
            signal: BidirectionalAmplitude {
                backward: 1e-12,
                ..BidirectionalAmplitude::default()
            },
            ..FieldState::default()
        }];
        let previous = vec![FieldState {
            signal: BidirectionalAmplitude {
                backward: -1e-12,
                ..BidirectionalAmplitude::default()
            },
            ..FieldState::default()
        }];

        let error = convergence_error(&current, &previous);

        assert!(error <= 1.0, "tiny zero crossing error was {error:e}");
    }

    #[test]
    fn localized_convergence_error_is_not_hidden() {
        let current = vec![
            FieldState {
                pump: BidirectionalAmplitude {
                    forward: 1.0,
                    ..BidirectionalAmplitude::default()
                },
                ..FieldState::default()
            };
            100
        ];
        let mut previous = current.clone();
        previous[50].pump.forward = 1.01;

        let error = convergence_error(&current, &previous);

        assert!(error > 1.0, "localized profile error was {error:e}");
    }

    #[test]
    fn convergence_fields_are_scaled_independently() {
        let current = vec![FieldState {
            signal: BidirectionalAmplitude {
                backward: 1e-6,
                ..BidirectionalAmplitude::default()
            },
            pump: BidirectionalAmplitude {
                forward: 1e6,
                ..BidirectionalAmplitude::default()
            },
        }];
        let previous = vec![FieldState {
            signal: BidirectionalAmplitude {
                backward: 2e-6,
                ..BidirectionalAmplitude::default()
            },
            pump: BidirectionalAmplitude {
                forward: 1e6,
                ..BidirectionalAmplitude::default()
            },
        }];

        let error = convergence_error(&current, &previous);

        assert!(error > 1.0, "signal error was hidden by pump scale");
    }

    #[test]
    fn non_finite_values_fail_convergence() {
        let current = vec![FieldState::default()];
        let previous = vec![FieldState {
            pump: BidirectionalAmplitude {
                backward: f64::NAN,
                ..BidirectionalAmplitude::default()
            },
            ..FieldState::default()
        }];

        assert!(convergence_error(&current, &previous).is_infinite());
    }

    #[test]
    fn two_level_dopant_returns_material_gain() {
        let dopant = TwoLevelDopant {
            density: 2.0,
            lifetime: 1.0,
            pump_ab: 3.0,
            pump_em: 0.0,
            sgnl_ab: 1.0,
            sgnl_em: 0.0,
        };

        let gain = dopant.gain(2.0, 3.0);

        assert!((gain.pump - -0.6).abs() < 1e-12);
        assert!((gain.signal - -0.2).abs() < 1e-12);
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
                    signal: BidirectionalAmplitude {
                        backward: -3.0,
                        ..BidirectionalAmplitude::default()
                    },
                    ..FieldState::default()
                },
                FieldState {
                    signal: BidirectionalAmplitude {
                        forward: 2.0,
                        ..BidirectionalAmplitude::default()
                    },
                    ..FieldState::default()
                },
            ],
        );

        assert_eq!(profile.output_powers(), (4.0, 9.0));
    }
}
