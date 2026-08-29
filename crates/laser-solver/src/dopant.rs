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

pub struct TwoLevelPopulations {
    pub(crate) ground: f64,
    pub(crate) excited: f64
}

struct TwoLevelRates {
    up: f64,
    down: f64
}

impl TwoLevelDopant {
    pub fn steady_state(&self, rates: TwoLevelRates) -> TwoLevelPopulations {
        let gamma_decay = 1.0 / self.lifetime;
        let gamma_dn_total = gamma_decay + rates.down;
        let denom = rates.up + gamma_dn_total;
        TwoLevelPopulations {
            ground: gamma_dn_total / denom,
            excited: rates.up / denom
        }
    }

    pub fn pops(&self, pump_flux: f64, sgnl_flux: f64) -> TwoLevelPopulations {
        let rates = TwoLevelRates {
            up: pump_flux * self.pump_ab + sgnl_flux * self.sgnl_ab,
            down: pump_flux * self.pump_em + sgnl_flux * self.sgnl_em
        };
        self.steady_state(rates)
    }
}

impl TwoLevelDopant {
    pub fn gain(&self, pump_flux: f64, sgnl_flux: f64) -> Gain {
        let pops = self.pops(pump_flux, sgnl_flux);
        let (g, e) = (pops.ground, pops.excited);
        self.gain_from_populations(g, e)
    }

    pub(crate) fn gain_from_populations(&self, ground: f64, excited: f64) -> Gain {
        Gain {
            pump: self.density * (-ground * self.pump_ab + excited * self.pump_em),
            signal: self.density * (-ground * self.sgnl_ab + excited * self.sgnl_em),
        }
    }
}