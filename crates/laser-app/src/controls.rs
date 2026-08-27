use eframe::egui;
use eframe::egui::Ui;
use laser_solver::dfb::Grating;
use laser_solver::lase::{Fibre, GridPoints, Pump};
use laser_solver::rootfind::BisectionConfig;

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

        ui.label("Lower bracket");
        changed |= ui
            .add(
                egui::Slider::new(&mut config.lower, 1e-9..=1e-2)
                    // 1e-10 causes slider field box to resize
                    .logarithmic(true)
                    .custom_formatter(|value, _| format!("{value:.1e}")),
            )
            .changed();
        ui.end_row();
    });

    changed
}

pub(crate) fn gridpoints_slider(gp: &mut GridPoints, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("grid-points").show(ui, |ui| {
        ui.label("Nz");
        changed |= ui
            .add(egui::Slider::new(&mut gp.0, 10..=1000).step_by(2.0))
            .changed();
        ui.end_row();
    });

    changed
}

pub(crate) fn grating_slider_grid(grating: &mut Grating, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("grating").show(ui, |ui| {
        ui.label("Kappa left");
        changed |= ui
            .add(egui::Slider::new(&mut grating.kappa_left, 0.1..=10.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Kappa right");
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
        ui.label("Total power");
        changed |= ui
            .add(egui::Slider::new(&mut pump.total, 0.0..=100.0).step_by(0.01))
            .changed();
        ui.end_row();

        ui.label("Balance");
        changed |= ui
            .add(egui::Slider::new(&mut pump.balance, -1.0..=1.0).step_by(0.01))
            .changed();
        ui.end_row();
    });

    changed
}

pub(crate) fn fibre_params_slider_grid(params: &mut Fibre, ui: &mut Ui) -> bool {
    let mut changed = false;

    egui::Grid::new("params").show(ui, |ui| {
        egui::Grid::new("params1").show(ui, |ui| {
            ui.label("Pump em.");
            changed |= ui
                .add(egui::Slider::new(&mut params.pump_em, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Pump abs.");
            changed |= ui
                .add(egui::Slider::new(&mut params.pump_ab, 0.05..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Signl em.");
            changed |= ui
                .add(egui::Slider::new(&mut params.sgnl_em, 0.05..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Signl abs.");
            changed |= ui
                .add(egui::Slider::new(&mut params.sgnl_ab, 0.0..=10.0).step_by(0.01))
                .changed();
            ui.end_row();
        });

        egui::Grid::new("params2").show(ui, |ui| {
            ui.label("Dopant density");
            changed |= ui
                .add(egui::Slider::new(&mut params.density, 0.1..=10.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Lifetime");
            changed |= ui
                .add(egui::Slider::new(&mut params.lifetime, 0.1..=2.0).step_by(0.01))
                .changed();
            ui.end_row();

            ui.label("Length");
            changed |= ui
                .add(egui::Slider::new(&mut params.length, 0.1..=15.0).step_by(0.01))
                .changed();
            ui.end_row();
        });
    });

    changed
}
