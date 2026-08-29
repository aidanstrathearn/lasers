use crate::fibre::BidirectionalAmplitude;
use crate::maths::utils::relative_diff;

use super::OutputPower;

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
        (right.signal.forward_power(), left.signal.backward_power())
    }
}
