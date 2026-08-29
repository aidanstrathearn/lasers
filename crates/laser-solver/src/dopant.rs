use crate::lase::Gain;


struct DopantError;
trait DopantModel {
    type Interaction;
    type Rates: Default;
    type Populations;

    fn interaction(
        &self,
        wavelength: f64,
    ) -> Result<Self::Interaction, DopantError>;

    fn rates(
        &self,
        interaction: &Self::Interaction,
        flux: f64,
    ) -> Self::Rates;

    fn add_rates(
        &self,
        rates: &mut Self::Rates,
        interaction: &Self::Interaction,
        flux: f64,
    );

    fn populations(&self, rates: &Self::Rates) -> Self::Populations;

    fn gain(
        &self,
        interaction: &Self::Interaction,
        populations: &Self::Populations,
    ) -> f64;
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

#[derive(Clone, Copy)]
pub struct TwoLevelPopulations {
    pub(crate) ground: f64,
    pub(crate) excited: f64,
}

#[derive(Clone, Copy, Default)]
pub struct TwoLevelRates {
    up: f64,
    down: f64,
}

#[derive(Clone, Copy)]
pub struct TwoLevelCrossSections {
    absorption: f64,
    emission: f64,
}

impl TwoLevelDopant {
    pub fn signal_cross_section(&self) -> TwoLevelCrossSections {
        TwoLevelCrossSections {
            absorption: self.sgnl_ab,
            emission: self.sgnl_em,
        }
    }

    pub fn pump_cross_section(&self) -> TwoLevelCrossSections {
        TwoLevelCrossSections {
            absorption: self.pump_ab,
            emission: self.pump_em,
        }
    }

    pub fn make_rates(data: &[(f64, TwoLevelCrossSections)]) -> TwoLevelRates {
        let mut rates = TwoLevelRates::default();
        for (flux, sigma) in data.iter() {
            rates = TwoLevelRates {
                up: rates.up + flux * sigma.absorption,
                down: rates.down + flux * sigma.emission,
            };
        }
        rates
    }


    pub fn steady_state(&self, rates: TwoLevelRates) -> TwoLevelPopulations {
        let gamma_decay = 1.0 / self.lifetime;
        let gamma_dn_total = gamma_decay + rates.down;
        let denom = rates.up + gamma_dn_total;
        TwoLevelPopulations {
            ground: gamma_dn_total / denom,
            excited: rates.up / denom,
        }
    }

    pub(crate) fn gain_from_crosssection(
        &self,
        pops: TwoLevelPopulations,
        sigma: TwoLevelCrossSections,
    ) -> f64 {
        self.density * (-pops.ground * sigma.absorption + pops.excited * sigma.emission)
    }



    pub fn pops(&self, data: &[(f64, TwoLevelCrossSections)]) -> TwoLevelPopulations {

        let rates = TwoLevelDopant::make_rates(data);
        self.steady_state(rates)
    }

    // pub fn gain(&self, data: &[(f64, TwoLevelCrossSections)]) -> Gain {
    //     // let data = &[
    //     //     (pump_flux, self.pump_cross_section()),
    //     //     (sgnl_flux, self.signal_cross_section()),
    //     // ];
    //     let pops = self.pops(data);
    //     Gain {
    //         pump: self.gain_from_crosssection(pops, self.pump_cross_section()),
    //         signal: self.gain_from_crosssection(pops, self.signal_cross_section()),
    //     }
    // }




}
