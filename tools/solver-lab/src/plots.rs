use crate::myplotlib::Plotter;
use laser_solver::lase::FieldProfile;

pub fn show_field_profile(profile: &FieldProfile, show_plot: bool) -> eframe::Result {
    if !show_plot {
        return Ok(());
    }

    let x: Vec<f64> = profile.z().collect();
    let clipped_log = |x: f64| x.powi(2).max(1e-6).log10();
    let pump_f: Vec<f64> = profile.pump_f().map(clipped_log).collect();
    let pump_b: Vec<f64> = profile.pump_b().map(clipped_log).collect();
    let sgnl_f: Vec<f64> = profile.sgnl_f().map(clipped_log).collect();
    let sgnl_b: Vec<f64> = profile.sgnl_b().map(clipped_log).collect();

    let mut plot = Plotter::new();
    plot.plot(&x, &pump_f).label("Forward Pump");
    plot.plot(&x, &pump_b).label("Backward Pump");
    plot.plot(&x, &sgnl_f).label("Forward Signal");
    plot.plot(&x, &sgnl_b).label("Backward Signal");
    plot.xlabel("z");
    plot.ylabel("log10(Power)");
    plot.title("Fields");
    plot.show()
}
