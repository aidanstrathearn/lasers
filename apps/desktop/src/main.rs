use eframe::egui;
use laser_app::LaserApp;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([800.0, 600.0])
            .with_min_inner_size([400.0, 300.0])
            .with_resizable(true),
        ..Default::default()
    };

    eframe::run_native(
        "Laser Solver",
        options,
        Box::new(|creation_context| Ok(Box::new(LaserApp::new(creation_context)))),
    )
}
