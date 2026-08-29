use crate::lase::Gain;

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
    pub fn steady_state(&self, rates: TwoLevelRates) -> TwoLevelPopulations {
        let gamma_decay = 1.0 / self.lifetime;
        let gamma_dn_total = gamma_decay + rates.down;
        let denom = rates.up + gamma_dn_total;
        TwoLevelPopulations {
            ground: gamma_dn_total / denom,
            excited: rates.up / denom,
        }
    }

    pub fn add_rate(
        flux: f64,
        sigma: TwoLevelCrossSections,
        rates: TwoLevelRates,
    ) -> TwoLevelRates {
        TwoLevelRates {
            up: rates.up + flux * sigma.absorption,
            down: rates.down + flux * sigma.emission,
        }
    }

    pub fn make_rates(data: &[(f64, TwoLevelCrossSections)]) -> TwoLevelRates {
        let mut rates = TwoLevelRates::default();
        for (flux, sigma) in data.iter() {
            rates = TwoLevelDopant::add_rate(*flux, *sigma, rates);
        }
        rates
    }

    pub fn pops(&self, pump_flux: f64, sgnl_flux: f64) -> TwoLevelPopulations {
        // let mut rates = TwoLevelRates::default();
        // rates = TwoLevelDopant::add_rate(pump_flux, self.pump_cross_section(), rates);
        // rates = TwoLevelDopant::add_rate(sgnl_flux, self.signal_cross_section(), rates);
        let rates = TwoLevelDopant::make_rates(&[
            (pump_flux, self.pump_cross_section()),
            (sgnl_flux, self.signal_cross_section()),
        ]);
        self.steady_state(rates)
    }

    pub fn gain(&self, pump_flux: f64, sgnl_flux: f64) -> Gain {
        let pops = self.pops(pump_flux, sgnl_flux);
        self.gain_from_populations(pops)
    }

    pub(crate) fn gain_from_crosssection(
        &self,
        pops: TwoLevelPopulations,
        sigma: TwoLevelCrossSections,
    ) -> f64 {
        self.density * (-pops.ground * sigma.absorption + pops.excited * sigma.emission)
    }

    pub(crate) fn gain_from_populations(&self, pops: TwoLevelPopulations) -> Gain {
        Gain {
            pump: self.gain_from_crosssection(pops, self.pump_cross_section()),
            signal: self.gain_from_crosssection(pops, self.signal_cross_section()),
        }
    }

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
}
