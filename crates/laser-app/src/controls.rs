use laser_solver::grating::PiShift;
use laser_solver::lase::{FibreGeometry, Pump, TwoLevelCrossSections, TwoLevelDopant};
use laser_solver::maths::rootfind::BisectionConfig;
use plot_app::{Slider, SliderGroup};
use std::ops::RangeInclusive;

pub(crate) fn milliwatt_slider<'a>(
    label: &'static str,
    power_watts: &'a mut f64,
    range_mw: RangeInclusive<f64>,
) -> Slider<'a> {
    Slider::from_get_set(label, range_mw, move |power_mw| {
        if let Some(power_mw) = power_mw {
            *power_watts = power_mw / 1_000.0;
        }
        1_000.0 * *power_watts
    })
}

pub(crate) fn solver_sliders<'a>(
    config: &'a mut BisectionConfig,
    steps: &'a mut usize,
) -> SliderGroup<'a> {
    let iteration = &mut config.iteration;
    SliderGroup::new(
        "Solver",
        [
            Slider::new("Iteration no.", &mut iteration.max, 10..=2000).step_by(10.0),
            Slider::new("Tolerance", &mut iteration.tol, 1e-9..=1e-2)
                // 1e-10 causes the slider field box to resize.
                .logarithmic(true)
                .custom_formatter(|value, _| format!("{value:.1e}")),
            Slider::new("Steps", steps, 10..=1000).step_by(2.0),
        ],
    )
}

pub(crate) fn grating_sliders(grating: &mut PiShift) -> SliderGroup<'_> {
    let PiShift {
        kappa_left,
        kappa_right,
        pi_shift_position,
    } = grating;

    SliderGroup::new(
        "Bragg",
        [
            Slider::new("Kappa left (1/m)", kappa_left, 0.1..=10.0).step_by(0.01),
            Slider::new("Kappa right (1/m)", kappa_right, 0.1..=10.0).step_by(0.01),
            Slider::new("Pi Shift Pos", pi_shift_position, 0.05..=0.95).step_by(0.01),
        ],
    )
}

pub(crate) fn pump_sliders(pump: &mut Pump) -> SliderGroup<'_> {
    let Pump { total, balance } = pump;

    SliderGroup::new(
        "Pump",
        [
            milliwatt_slider("Total power (mW)", total, 0.0..=100.0).step_by(0.001),
            Slider::new("Balance", balance, -1.0..=1.0).step_by(0.01),
        ],
    )
}

pub(crate) fn fibre_param_sliders<'a>(
    geometry: &'a mut FibreGeometry,
    dopant: &'a mut TwoLevelDopant,
    pump_interaction: &'a mut TwoLevelCrossSections,
    signal_interaction: &'a mut TwoLevelCrossSections,
) -> SliderGroup<'a> {
    let FibreGeometry { length, .. } = geometry;
    let TwoLevelDopant { density, lifetime } = dopant;
    let TwoLevelCrossSections {
        absorption: pump_absorption,
        emission: pump_emission,
    } = pump_interaction;
    let TwoLevelCrossSections {
        absorption: signal_absorption,
        emission: signal_emission,
    } = signal_interaction;

    SliderGroup::new(
        "Fibre",
        [
            Slider::new("Pump em. (1e-25 m²)", pump_emission, 0.0..=10.0).step_by(0.01),
            Slider::new("Pump abs. (1e-25 m²)", pump_absorption, 0.05..=10.0).step_by(0.01),
            Slider::new("Signl em. (1e-25 m²)", signal_emission, 0.05..=10.0).step_by(0.01),
            Slider::new("Signl abs. (1e-25 m²)", signal_absorption, 0.0..=10.0).step_by(0.01),
            Slider::new("Dopant density (1e25/m³)", density, 0.1..=10.0).step_by(0.01),
            Slider::new("Lifetime (s)", lifetime, 0.01..=2.0).step_by(0.01),
            Slider::new("Length (m)", length, 0.1..=15.0).step_by(0.01),
        ],
    )
}
