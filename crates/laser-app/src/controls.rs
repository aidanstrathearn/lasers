use eframe::egui;
use eframe::egui::Ui;
use laser_solver::grating::{GratingModel, PiShift};
use laser_solver::lase::{Fibre, Pump, TwoLevelCrossSections, TwoLevelDopant};
use laser_solver::maths::rootfind::BisectionConfig;

pub(crate) fn power_slider_mw(power_watts: &mut f64, ui: &mut Ui) -> bool {
    let mut power_mw = 1_000.0 * *power_watts;
    let changed = ui
        .add(egui::Slider::new(&mut power_mw, 0.0..=10.0).step_by(0.001))
        .changed();
    *power_watts = power_mw / 1_000.0;
    changed
}

pub(crate) fn bisection_slider_grid(config: &mut BisectionConfig, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("bisection").show(ui, |ui| {
        ui.label("Iteration no.");
        changed |= ui
            .add(egui::Slider::new(&mut config.iteration.max, 10..=2000).step_by(10.0))
            .changed();
        ui.end_row();

        ui.label("Tolerance");
        changed |= ui
            .add(
                egui::Slider::new(&mut config.iteration.tol, 1e-9..=1e-2)
                    // 1e-10 causes slider field box to resize
                    .logarithmic(true)
                    .custom_formatter(|value, _| format!("{value:.1e}")),
            )
            .changed();
        ui.end_row();
        
    });

    changed
}

pub(crate) fn steps_slider(steps: &mut usize, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("grid-steps").show(ui, |ui| {
        ui.label("Steps");
        changed |= ui
            .add(egui::Slider::new(steps, 10..=1000).step_by(2.0))
            .changed();
        ui.end_row();
    });

    changed
}

pub(crate) fn grating_slider_grid(grating: &mut PiShift, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("grating").show(ui, |ui| {
        ui.label("Kappa left (m⁻¹)");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_left, 0.1..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Kappa right (m⁻¹)");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_right, 0.1..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Pi Shift Pos");
        changed |= ui
            .add(egui::Slider::new(&mut grating.pi_shift_position, 0.05..=0.95).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}

pub(crate) fn pump_slider_grid(pump: &mut Pump, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("pump").show(ui, |ui| {
        ui.label("Total power (mW)");
        changed |= power_slider_mw(&mut pump.total, ui);
        ui.end_row();

        ui.label("Balance");
        changed |= ui
            .add(egui::Slider::new(&mut pump.balance, -1.0..=1.0).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}

pub(crate) fn fibre_params_slider_grid<G: GratingModel>(
    params: &mut Fibre<TwoLevelDopant, G>,
    pump_interaction: &mut TwoLevelCrossSections,
    signal_interaction: &mut TwoLevelCrossSections,
    ui: &mut Ui,
) -> bool {
    let mut changed = false;
    let dopant = &mut params.dopant;

    egui::Grid::new("params").show(ui, |ui| {
        egui::Grid::new("params1").show(ui, |ui| {
            ui.label("Pump em. (10⁻²⁵ m²)");
            changed |= ui
                .add(egui::Slider::new(&mut pump_interaction.emission, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Pump abs. (10⁻²⁵ m²)");
            changed |= ui
                .add(egui::Slider::new(&mut pump_interaction.absorption, 0.05..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Signl em. (10⁻²⁵ m²)");
            changed |= ui
                .add(egui::Slider::new(&mut signal_interaction.emission, 0.05..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Signl abs. (10⁻²⁵ m²)");
            changed |= ui
                .add(
                    egui::Slider::new(&mut signal_interaction.absorption, 0.0..=10.0).step_by(0.01),
                )
                .changed();
            ui.end_row();
        });

        egui::Grid::new("params2").show(ui, |ui| {
            ui.label("Dopant density (10²⁵ m⁻³)");
            changed |= ui
                .add(egui::Slider::new(&mut dopant.density, 0.1..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Lifetime (s)");
            changed |= ui
                .add(egui::Slider::new(&mut dopant.lifetime, 0.1..=2.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Length (m)");
            changed |= ui
                .add(egui::Slider::new(&mut params.geometry.length, 0.1..=15.0).step_by(0.01))
                .changed();
            ui.end_row();
        });
    });

    changed
}
