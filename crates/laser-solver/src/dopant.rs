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

    pub fn pops(&self, pump_flux: f64, sgnl_flux: f64) -> (f64, f64) {
        let gamma_up = pump_flux * self.pump_ab + sgnl_flux * self.sgnl_ab;
        let gamma_dn = pump_flux * self.pump_em + sgnl_flux * self.sgnl_em;
        self.steady_state(gamma_dn, gamma_up)
    }
}
