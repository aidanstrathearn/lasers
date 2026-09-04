use crate::controls::milliwatt_slider;
use crate::dfb::LaserParameters;
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::maths::utils::linspace;
use laser_solver::two_mode::TwoModeSolver;
use plot_app::{AppResult, Plotter, Points, Slider, SliderGroup};

#[derive(Copy, Clone)]
pub struct ThresholdRange {
    lower: f64,
    upper: f64,
    num: usize,
}

impl ThresholdRange {
    pub(crate) fn new_watts(lower: f64, upper: f64, num: usize) -> Self {
        assert!(lower.is_finite() && lower >= 0.0);
        assert!(upper.is_finite() && upper >= lower);
        assert!(num > 1);
        Self { lower, upper, num }
    }

    pub(crate) fn bisection_config(self, config: BisectionConfig) -> BisectionConfig {
        BisectionConfig {
            lower: self.lower,
            upper: self.upper,
            ..config
        }
    }
}

impl LaserParameters {
    pub fn threshold_plot(&mut self) -> AppResult {
        let pumps = linspace(
            self.threshold_range.lower,
            self.threshold_range.upper,
            self.threshold_range.num,
        );
        let threshold = {
            let fibre = self.resolved_fibre();
            let bc = BisectionConfig {
                upper: fibre.pump_flux(self.threshold_range.upper).sqrt(),
                ..self.config
            };
            TwoModeSolver::new(&fibre, self.steps).pump_scan(
                &pumps,
                self.pump.balance,
                bc.into(),
                self.picard_config,
            )
        };
        let threshold = threshold?;
        let sgnl_f = threshold
            .iter()
            .map(|output| output.as_ref().map_or(0.0, |output| output.0));
        let sgnl_b = threshold
            .iter()
            .map(|output| output.as_ref().map_or(0.0, |output| output.1));

        let sgnl_f_points: Points = pumps
            .iter()
            .zip(sgnl_f)
            .map(|(&x, y)| [1_000.0 * x, 1_000.0 * y.abs()])
            .collect();
        let sgnl_b_points: Points = pumps
            .iter()
            .zip(sgnl_b)
            .map(|(&x, y)| [1_000.0 * x, 1_000.0 * y.abs()])
            .collect();

        let mut plt = Plotter::new();
        plt.xlabel("Pump power (mW)");
        plt.ylabel("Output power (mW)");
        plt.add_points(sgnl_f_points).label("Forward");
        plt.add_points(sgnl_b_points).label("Backward");
        plt.axvline(1_000.0 * self.pump.total).label("Current pump");
        plt.xlim(
            1_000.0 * self.threshold_range.lower,
            1_000.0 * self.threshold_range.upper,
        );
        Ok(plt)
    }
}

pub(crate) fn threshold_sliders(tr: &mut ThresholdRange) -> SliderGroup<'_> {
    let upper_mw = 1_000.0 * tr.upper;
    let upper_limit_mw = upper_mw.max(100.0);
    let ThresholdRange { lower, upper, num } = tr;

    SliderGroup::new(
        "Threshold",
        [
            milliwatt_slider("High (mW)", upper, 1e-9..=upper_limit_mw).step_by(0.001),
            milliwatt_slider("Low (mW)", lower, 1e-10..=upper_mw).logarithmic(true),
            Slider::new("num", num, 5..=100).step_by(0.01),
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_threshold_range_preserves_legacy_flux_bounds() {
        let mode = LaserParameters::default();
        let fibre = mode.resolved_fibre();

        let lower_flux = fibre.pump_flux(mode.threshold_range.lower);
        let upper_flux = fibre.pump_flux(mode.threshold_range.upper);

        assert!((lower_flux / 1e-6 - 1.0).abs() < 1e-12);
        assert!((upper_flux / 10.0 - 1.0).abs() < 1e-12);
    }
}
