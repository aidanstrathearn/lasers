use crate::myplotlib::Plotter;
use laser_solver::lase::FieldProfile;

pub fn show_field_profile(profile: &FieldProfile, show_plot: bool) -> eframe::Result {
    if !show_plot {
        return Ok(());
    }

    let x: Vec<f64> = profile.z().collect();
    let clipped_log = |power_watts: f64| power_watts.max(1e-12).log10();
    let pump_f: Vec<f64> = profile.pump_forward_power().map(clipped_log).collect();
    let pump_b: Vec<f64> = profile.pump_backward_power().map(clipped_log).collect();
    let sgnl_f: Vec<f64> = profile.signal_forward_power().map(clipped_log).collect();
    let sgnl_b: Vec<f64> = profile.signal_backward_power().map(clipped_log).collect();

    let mut plot = Plotter::new();
    plot.plot(&x, &pump_f).label("Forward Pump");
    plot.plot(&x, &pump_b).label("Backward Pump");
    plot.plot(&x, &sgnl_f).label("Forward Signal");
    plot.plot(&x, &sgnl_b).label("Backward Signal");
    plot.xlabel("Position (m)");
    plot.ylabel("log10(Power (W))");
    plot.title("Fields");
    plot.show()
}
