use crate::difference;
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

pub fn plot_profile_diff(
    profile: &FieldProfile,
    other: &FieldProfile,
    title: &str,
) -> eframe::Result {
    let z: Vec<f64> = profile.z().collect();

    let sgnl_f: Vec<f64> = profile.sgnl_f().collect();
    let other_sgnl_f: Vec<f64> = other.sgnl_f().collect();
    let diff_sgnl_f = difference(&sgnl_f, &other_sgnl_f);

    let sgnl_b: Vec<f64> = profile.sgnl_b().collect();
    let other_sgnl_b: Vec<f64> = other.sgnl_b().collect();
    let diff_sgnl_b = difference(&sgnl_b, &other_sgnl_b);

    let pump_b: Vec<f64> = profile.pump_b().collect();
    let other_pump_b: Vec<f64> = other.pump_b().collect();
    let diff_pump_b = difference(&pump_b, &other_pump_b);

    let pump_f: Vec<f64> = profile.pump_f().collect();
    let other_pump_f: Vec<f64> = other.pump_f().collect();
    let diff_pump_f = difference(&pump_f, &other_pump_f);

    let mut plot = Plotter::new();
    plot.plot(&z, &diff_sgnl_b).label("backward sgnl diff");
    plot.plot(&z, &diff_sgnl_f).label("forward sgnl diff");
    plot.plot(&z, &diff_pump_b).label("backward pump diff");
    plot.plot(&z, &diff_pump_f).label("forward pump diff");
    plot.xlabel("z");
    plot.ylabel("Field amplitude diff");
    plot.title(title);
    plot.show()
}
