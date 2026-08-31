use crate::dopant::{DopantError, DopantModel, TwoLevelDopant};
use crate::fibre::{ActiveMode, Fibre, FieldMode};
use crate::grating::{GratingModel, NoGrating};

use super::fieldstate::{FieldProfile, FieldState};

pub type OutputPower = (f64, f64);

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gain {
    pub pump: f64,
    pub signal: f64,
}

impl<D: DopantModel, G: GratingModel> Fibre<D, G> {
    pub fn resolve(
        &self,
        pump_mode: FieldMode,
        sgnl_mode: FieldMode,
    ) -> Result<ResolvedFibre<'_, D, G>, DopantError> {
        let pump_interaction = self.dopant.interaction(pump_mode.wavelength())?;
        let sgnl_interaction = self.dopant.interaction(sgnl_mode.wavelength())?;
        Ok(
            self.resolve_with_interactions(
                pump_mode,
                pump_interaction,
                sgnl_mode,
                sgnl_interaction,
            ),
        )
    }

    pub fn resolve_with_interactions(
        &self,
        pump_mode: FieldMode,
        pump_interaction: D::Interaction,
        sgnl_mode: FieldMode,
        sgnl_interaction: D::Interaction,
    ) -> ResolvedFibre<'_, D, G> {
        ResolvedFibre {
            fibre: self,
            pump: ActiveMode::new(self, pump_mode, pump_interaction),
            sgnl: ActiveMode::new(self, sgnl_mode, sgnl_interaction),
        }
    }
}



pub struct ResolvedFibre<'a, D: DopantModel = TwoLevelDopant, G: GratingModel = NoGrating> {
    fibre: &'a Fibre<D, G>,
    pump: ActiveMode<D>,
    sgnl: ActiveMode<D>,
}

impl<D: DopantModel, G: GratingModel> ResolvedFibre<'_, D, G> {
    pub fn length(&self) -> f64 {
        self.fibre.geometry.length
    }

    pub fn pump_mode(&self) -> FieldMode {
        self.pump.mode
    }

    pub fn sgnl_mode(&self) -> FieldMode {
        self.sgnl.mode
    }

    pub fn pump_overlap(&self) -> f64 {
        self.pump.overlap
    }

    pub fn sgnl_overlap(&self) -> f64 {
        self.sgnl.overlap
    }

    pub fn grating(&self) -> &G {
        &self.fibre.grating
    }

    pub fn pump_flux(&self, power_watts: f64) -> f64 {
        self.pump.flux_from_power(power_watts)
    }

    pub fn signal_flux(&self, power_watts: f64) -> f64 {
        self.sgnl.flux_from_power(power_watts)
    }

    pub fn pump_power(&self, flux: f64) -> f64 {
        self.pump.power_from_flux(flux)
    }

    pub fn signal_power(&self, flux: f64) -> f64 {
        self.sgnl.power_from_flux(flux)
    }

    pub(crate) fn pump_flux_amplitudes(&self, pump: Pump) -> (f64, f64) {
        pump.validate();
        self.pump.flux_amplitudes(pump.total, pump.balance)
    }

    pub(crate) fn signal_flux_amplitudes(&self, signal: Signal) -> (f64, f64) {
        signal.validate();
        self.sgnl.flux_amplitudes(signal.total, signal.balance)
    }

    pub(crate) fn flux_per_watt(&self) -> (f64, f64) {
        (self.pump.flux_per_watt, self.sgnl.flux_per_watt)
    }

    pub fn mode_fluxes(&self, fs: FieldState) -> (f64, f64) {
        (
            fs.pump.total_flux() * self.pump.overlap,
            fs.signal.total_flux() * self.sgnl.overlap,
        )
    }

    fn rates(&self, fs: FieldState) -> D::Rates {
        let (pump_flux, sgnl_flux) = self.mode_fluxes(fs);
        let mut rates = D::Rates::default();
        self.fibre
            .dopant
            .add_rates(&mut rates, &self.pump.interaction, pump_flux);
        self.fibre
            .dopant
            .add_rates(&mut rates, &self.sgnl.interaction, sgnl_flux);
        rates
    }

    pub fn gain(&self, fs: FieldState) -> Gain {
        let rates = self.rates(fs);
        let populations = self.fibre.dopant.populations(&rates);
        let mut gain = Gain {
            pump: self.fibre.dopant.gain(&self.pump.interaction, &populations),
            signal: self.fibre.dopant.gain(&self.sgnl.interaction, &populations),
        };

        gain.pump *= self.pump.overlap;
        gain.signal *= self.sgnl.overlap;
        gain
    }

    pub fn populations(&self, fs: FieldState) -> D::Populations {
        let rates = self.rates(fs);
        self.fibre.dopant.populations(&rates)
    }

    pub fn profile_populations(&self, profile: &FieldProfile) -> Vec<D::Populations> {
        profile
            .fields
            .iter()
            .map(|&field| self.populations(field))
            .collect()
    }
}

impl<D, G> Clone for ResolvedFibre<'_, D, G>
where
    D: DopantModel,
    G: GratingModel,
    D::Interaction: Clone,
{
    fn clone(&self) -> Self {
        Self {
            fibre: self.fibre,
            pump: self.pump.clone(),
            sgnl: self.sgnl.clone(),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Signal {
    /// Total injected signal power in watts.
    pub total: f64,
    pub balance: f64,
}

impl Default for Signal {
    fn default() -> Self {
        Self {
            total: 1e-4,
            balance: 1.0,
        }
    }
}

impl Signal {
    fn validate(self) {
        assert!(
            self.total >= 0.0 && (-1.0..=1.0).contains(&self.balance),
            "signal total must be non-negative and balance must be between -1 and 1"
        );
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Pump {
    /// Total injected pump power in watts.
    pub total: f64,
    pub balance: f64,
}

impl Default for Pump {
    fn default() -> Self {
        Self {
            total: 1e-2,
            balance: 1.0,
        }
    }
}

impl Pump {
    fn validate(self) {
        assert!(
            self.total >= 0.0 && (-1.0..=1.0).contains(&self.balance),
            "pump total must be non-negative and balance must be between -1 and 1"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fibre::BidirectionalAmplitude;
    use crate::two_mode::fieldstate::{FieldProfile, profile_convergence_error};

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
    fn active_mode_converts_power_and_balance_to_flux_amplitudes() {
        let fibre = Fibre::<TwoLevelDopant, NoGrating>::default();
        let fibre = fibre.resolve_with_interactions(
            FieldMode::new(970e-9),
            crate::dopant::TwoLevelCrossSections::new(1.0, 0.0),
            FieldMode::new(1060e-9),
            crate::dopant::TwoLevelCrossSections::new(0.0, 1.0),
        );
        let hundred_flux_power = fibre.pump_power(100.0);

        for (balance, expected) in [
            (1.0, (10.0, 0.0)),
            (-1.0, (0.0, 10.0)),
            (0.0, (50.0_f64.sqrt(), 50.0_f64.sqrt())),
        ] {
            assert_eq!(
                fibre.pump_flux_amplitudes(Pump {
                    total: hundred_flux_power,
                    balance,
                }),
                expected
            );
        }

        assert_eq!(
            fibre.pump_flux_amplitudes(Pump {
                total: 0.0,
                balance: 0.25,
            }),
            (0.0, 0.0)
        );
    }

    #[test]
    fn active_mode_uses_scaled_photon_flux_per_watt() {
        let fibre = Fibre::<TwoLevelDopant, NoGrating>::default();
        let fibre = fibre.resolve_with_interactions(
            FieldMode::new(970e-9),
            crate::dopant::TwoLevelCrossSections::new(1.0, 0.0),
            FieldMode::new(1060e-9),
            crate::dopant::TwoLevelCrossSections::new(0.0, 1.0),
        );

        let expected = 9_714.604_996_881;
        assert!((fibre.pump_flux(1.0) / expected - 1.0).abs() < 1e-12);
        assert!((fibre.pump_power(expected) - 1.0).abs() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "pump total must be non-negative")]
    fn pump_rejects_negative_total_power() {
        Pump {
            total: -1.0,
            balance: 0.0,
        }
        .validate();
    }

    #[test]
    #[should_panic(expected = "balance must be between -1 and 1")]
    fn pump_rejects_out_of_range_balance() {
        Pump {
            total: 1.0,
            balance: 2.0,
        }
        .validate();
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
            5.0,
            2.0,
        );

        assert_eq!(profile.output_powers(), (2.0, 4.5));
    }
}
