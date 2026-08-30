use crate::dopant::{DopantError, DopantModel, TwoLevelDopant};
use crate::fibre::{Fibre, FieldMode, bidirectional_amplitudes};
use crate::grating::{GratingModel, NoGrating};

use super::fieldstate::FieldState;

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
            pump_mode,
            sgnl_mode,
            pump_overlap: self.geometry.mode_overlap(pump_mode),
            sgnl_overlap: self.geometry.mode_overlap(sgnl_mode),
            pump_interaction,
            sgnl_interaction,
        }
    }
}

pub struct ResolvedFibre<
    'a,
    D: DopantModel = TwoLevelDopant,
    G: GratingModel = NoGrating,
> {
    fibre: &'a Fibre<D, G>,
    pump_mode: FieldMode,
    sgnl_mode: FieldMode,
    pump_overlap: f64,
    sgnl_overlap: f64,
    pump_interaction: D::Interaction,
    sgnl_interaction: D::Interaction,
}

impl<D: DopantModel, G: GratingModel> ResolvedFibre<'_, D, G> {
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

    pub fn grating(&self) -> &G {
        &self.fibre.grating
    }

    pub fn mode_fluxes(&self, fs: FieldState) -> (f64, f64) {
        (
            fs.pump.total_power() * self.pump_overlap,
            fs.signal.total_power() * self.sgnl_overlap,
        )
    }

    fn rates(&self, fs: FieldState) -> D::Rates {
        let (pump_flux, sgnl_flux) = self.mode_fluxes(fs);
        let mut rates = D::Rates::default();
        self.fibre
            .dopant
            .add_rates(&mut rates, &self.pump_interaction, pump_flux);
        self.fibre
            .dopant
            .add_rates(&mut rates, &self.sgnl_interaction, sgnl_flux);
        rates
    }

    pub fn gain(&self, fs: FieldState) -> Gain {
        let rates = self.rates(fs);
        let populations = self.fibre.dopant.populations(&rates);
        let mut gain = Gain {
            pump: self.fibre.dopant.gain(&self.pump_interaction, &populations),
            signal: self.fibre.dopant.gain(&self.sgnl_interaction, &populations),
        };

        gain.pump = gain.pump * self.pump_overlap;
        gain.signal = gain.signal * self.sgnl_overlap;
        gain
    }

    pub fn populations(&self, fs: FieldState) -> D::Populations {
        let rates = self.rates(fs);
        self.fibre.dopant.populations(&rates)
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
            pump_mode: self.pump_mode,
            sgnl_mode: self.sgnl_mode,
            pump_overlap: self.pump_overlap,
            sgnl_overlap: self.sgnl_overlap,
            pump_interaction: self.pump_interaction.clone(),
            sgnl_interaction: self.sgnl_interaction.clone(),
        }
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
    use crate::fibre::BidirectionalAmplitude;
    use crate::two_mode::fieldstate::{profile_convergence_error, FieldProfile};

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
