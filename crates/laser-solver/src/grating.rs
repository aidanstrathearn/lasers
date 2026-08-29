pub trait GratingModel {
    fn kappa(&self, position: f64) -> f64;
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct NoGrating;

impl GratingModel for NoGrating {
    fn kappa(&self, _position: f64) -> f64 {
        0.0
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PiShift {
    pub kappa_left: f64,
    pub kappa_right: f64,
    pub pi_shift_position: f64,
}

impl Default for PiShift {
    fn default() -> Self {
        Self {
            kappa_left: 1.0,
            kappa_right: 1.0,
            pi_shift_position: 0.45,
        }
    }
}

impl GratingModel for PiShift {
    fn kappa(&self, position: f64) -> f64 {
        if position < self.pi_shift_position {
            self.kappa_left
        } else {
            -self.kappa_right
        }
    }
}

pub fn sample_grating<G: GratingModel>(grating: &G, steps: usize) -> Vec<f64> {
    assert!(steps > 0, "grating profile must contain at least one step");
    (0..steps)
        .map(|step| grating.kappa(step as f64 / steps as f64))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_grating_has_zero_coupling() {
        assert_eq!(NoGrating.kappa(0.0), 0.0);
        assert_eq!(NoGrating.kappa(0.5), 0.0);
        assert_eq!(NoGrating.kappa(1.0), 0.0);
    }

    #[test]
    fn pi_shift_changes_sign_at_shift_position() {
        let grating = PiShift {
            kappa_left: 2.0,
            kappa_right: 3.0,
            pi_shift_position: 0.4,
        };

        assert_eq!(grating.kappa(0.39), 2.0);
        assert_eq!(grating.kappa(0.4), -3.0);
    }

    #[test]
    fn sampling_uses_step_left_edges() {
        let grating = PiShift {
            kappa_left: 2.0,
            kappa_right: 3.0,
            pi_shift_position: 0.5,
        };

        assert_eq!(
            sample_grating(&grating, 4),
            vec![2.0, 2.0, -3.0, -3.0]
        );
    }
}
