use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DopantError {
    InteractionUnavailable { wavelength: f64 },
}

impl Display for DopantError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InteractionUnavailable { wavelength } => {
                write!(
                    f,
                    "dopant interaction unavailable at wavelength {wavelength:e} m"
                )
            }
        }
    }
}

impl Error for DopantError {}

pub trait DopantModel {
    type Interaction;
    type Rates: Default;
    type Populations;

    fn interaction(&self, wavelength: f64) -> Result<Self::Interaction, DopantError>;

    fn rates(&self, interaction: &Self::Interaction, flux: f64) -> Self::Rates;

    fn add_rates(&self, rates: &mut Self::Rates, interaction: &Self::Interaction, flux: f64);

    fn populations(&self, rates: &Self::Rates) -> Self::Populations;

    fn gain(&self, interaction: &Self::Interaction, populations: &Self::Populations) -> f64;
}

#[derive(Clone, Copy)]
pub struct TwoLevelDopant {
    pub density: f64,
    pub lifetime: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoLevelPopulations {
    pub ground: f64,
    pub excited: f64,
}

#[derive(Clone, Copy, Default)]
pub struct TwoLevelRates {
    up: f64,
    down: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoLevelCrossSections {
    pub absorption: f64,
    pub emission: f64,
}

impl TwoLevelCrossSections {
    pub const fn new(absorption: f64, emission: f64) -> Self {
        Self {
            absorption,
            emission,
        }
    }
}

impl TwoLevelDopant {
    pub fn steady_state(&self, rates: TwoLevelRates) -> TwoLevelPopulations {
        let gamma_decay = 1.0 / self.lifetime;
        let gamma_dn_total = gamma_decay + rates.down;
        let denom = rates.up + gamma_dn_total;
        TwoLevelPopulations {
            ground: gamma_dn_total / denom,
            excited: rates.up / denom,
        }
    }
}

impl DopantModel for TwoLevelDopant {
    type Interaction = TwoLevelCrossSections;
    type Rates = TwoLevelRates;
    type Populations = TwoLevelPopulations;

    fn interaction(&self, wavelength: f64) -> Result<Self::Interaction, DopantError> {
        Err(DopantError::InteractionUnavailable { wavelength })
    }

    fn rates(&self, interaction: &Self::Interaction, flux: f64) -> Self::Rates {
        TwoLevelRates {
            up: flux * interaction.absorption,
            down: flux * interaction.emission,
        }
    }

    fn add_rates(&self, rates: &mut Self::Rates, interaction: &Self::Interaction, flux: f64) {
        let contribution = self.rates(interaction, flux);
        rates.up += contribution.up;
        rates.down += contribution.down;
    }

    fn populations(&self, rates: &Self::Rates) -> Self::Populations {
        self.steady_state(*rates)
    }

    fn gain(&self, interaction: &Self::Interaction, populations: &Self::Populations) -> f64 {
        self.density
            * (-populations.ground * interaction.absorption
                + populations.excited * interaction.emission)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_level_model_accumulates_rates_and_returns_material_gain() {
        let dopant = TwoLevelDopant {
            density: 2.0,
            lifetime: 1.0,
        };
        let pump = TwoLevelCrossSections::new(3.0, 0.0);
        let signal = TwoLevelCrossSections::new(1.0, 0.0);
        let mut rates = TwoLevelRates::default();
        dopant.add_rates(&mut rates, &pump, 2.0);
        dopant.add_rates(&mut rates, &signal, 3.0);

        let populations = dopant.populations(&rates);

        assert!((dopant.gain(&pump, &populations) - -0.6).abs() < 1e-12);
        assert!((dopant.gain(&signal, &populations) - -0.2).abs() < 1e-12);
    }
}
