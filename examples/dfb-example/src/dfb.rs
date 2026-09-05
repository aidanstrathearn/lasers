use crate::controls::{fibre_param_sliders, grating_sliders, pump_sliders, solver_sliders};
use crate::threshold_plot::{ThresholdRange, threshold_sliders};
use laser_solver::grating::PiShift;
use laser_solver::lase::{
    Fibre, FibreGeometry, FieldMode, FieldProfile, Pump, ResolvedFibre, TwoLevelCrossSections,
    TwoLevelDopant,
};
use laser_solver::maths::rootfind::BisectionConfig;
use laser_solver::two_mode::PicardConfig;
use myplotlib::{AppDefinition, Plotter, Points, SliderGrid, ViewOption};

pub(crate) struct LaserParameters {
    pub(crate) fibre: Fibre<TwoLevelDopant, PiShift>,
    pub(crate) pump_mode: FieldMode,
    pub(crate) sgnl_mode: FieldMode,
    pub(crate) pump_interaction: TwoLevelCrossSections,
    pub(crate) signal_interaction: TwoLevelCrossSections,
    pub(crate) pump: Pump,
    pub(crate) steps: usize,
    pub(crate) config: BisectionConfig,
    pub(crate) picard_config: PicardConfig,
    pub(crate) threshold_range: ThresholdRange,
}

impl LaserParameters {
    pub(crate) fn resolved_fibre(&self) -> ResolvedFibre<'_, TwoLevelDopant, PiShift> {
        self.fibre.resolve_with_interactions(
            self.pump_mode,
            self.pump_interaction,
            self.sgnl_mode,
            self.signal_interaction,
        )
    }
}

impl Default for LaserParameters {
    fn default() -> Self {
        let fibre = Fibre {
            geometry: FibreGeometry {
                core_radius: 4e-6,
                numerical_aperture: 0.1,
                length: 5.0,
            },
            dopant: TwoLevelDopant {
                density: 0.50,
                lifetime: 1.0,
            },
            grating: PiShift {
                kappa_left: 0.6,
                kappa_right: 0.6,
                pi_shift_position: 0.5,
            },
        };
        let pump_mode = FieldMode::new(970e-9);
        let sgnl_mode = FieldMode::new(1060e-9);
        let pump_interaction = TwoLevelCrossSections::new(1.0, 0.0);
        let signal_interaction = TwoLevelCrossSections::new(0.0, 1.0);
        let resolved_fibre = fibre.resolve_with_interactions(
            pump_mode,
            pump_interaction,
            sgnl_mode,
            signal_interaction,
        );
        let pump_total = resolved_fibre.pump_power(10.0);
        let threshold_range = ThresholdRange::new_watts(
            resolved_fibre.pump_power(1e-6),
            resolved_fibre.pump_power(10.0),
            20,
        );

        Self {
            fibre,
            pump_mode,
            sgnl_mode,
            pump_interaction,
            signal_interaction,
            pump: Pump {
                total: pump_total,
                balance: 1.0,
            },
            steps: 100,
            config: BisectionConfig::default(),
            picard_config: PicardConfig {
                max_iterations: 5_000,
                relative_tolerance: 1e-6,
                absolute_tolerance: 1e-10,
                relaxation: 1.0,
            },
            threshold_range,
        }
    }
}

fn power_points(z: impl Iterator<Item = f64>, power_watts: impl Iterator<Item = f64>) -> Points {
    z.zip(power_watts)
        .map(|(z, power_watts)| [z, 1_000.0 * power_watts])
        .collect()
}

pub(crate) fn field_profile_plot(profile: &FieldProfile) -> Plotter {
    let mut plot = Plotter::new();
    plot.add_points(power_points(profile.z(), profile.signal_forward_power()))
        .label("Forward signal");
    plot.add_points(power_points(profile.z(), profile.signal_backward_power()))
        .label("Backward signal");
    plot.add_points(power_points(profile.z(), profile.pump_forward_power()))
        .label("Forward pump");
    plot.add_points(power_points(profile.z(), profile.pump_backward_power()))
        .label("Backward pump");
    plot.xlabel("Position (m)");
    plot.ylabel("Power (mW)");
    plot
}



impl LaserParameters {
    fn standard_controls(&mut self) -> SliderGrid<'_> {
        self.make_controls(false)
    }

    fn threshold_controls(&mut self) -> SliderGrid<'_> {
        self.make_controls(true)
    }

    fn make_controls(&mut self, include_threshold: bool) -> SliderGrid<'_> {
        let Self {
            fibre,
            pump_interaction,
            signal_interaction,
            pump,
            steps,
            config,
            threshold_range,
            ..
        } = self;
        let Fibre {
            geometry,
            dopant,
            grating,
        } = fibre;

        let mut groups = vec![
            fibre_param_sliders(geometry, dopant, pump_interaction, signal_interaction),
            grating_sliders(grating),
            pump_sliders(pump),
            solver_sliders(config, steps),
        ];
        if include_threshold {
            groups.push(threshold_sliders(threshold_range));
        }

        SliderGrid::new(4, groups)
    }
}



const DFB_VIEWS: &[ViewOption<LaserParameters>] = &[
    ViewOption::new(
        "Profile",
        LaserParameters::profile_plot,
        LaserParameters::standard_controls,
    ),
    ViewOption::new(
        "Populations",
        LaserParameters::pops_plot,
        LaserParameters::standard_controls,
    ),
    ViewOption::new(
        "Kappa",
        LaserParameters::kappa_plot,
        LaserParameters::standard_controls,
    ),
    ViewOption::new(
        "Threshold",
        LaserParameters::threshold_plot,
        LaserParameters::threshold_controls,
    ),
    ViewOption::new(
        "Pi position",
        LaserParameters::pi_pos_plot,
        LaserParameters::standard_controls,
    ),
    ViewOption::new(
        "Pi threshold",
        LaserParameters::pi_pos_threshold_plot,
        LaserParameters::threshold_controls,
    ),
];

const DFB_APP: AppDefinition<LaserParameters> =
    AppDefinition::new("Laser Solver", "plot-canvas", DFB_VIEWS);

#[cfg(not(target_arch = "wasm32"))]
pub fn run_native() -> myplotlib::NativeResult {
    myplotlib::run_native(DFB_APP)
}

#[cfg(target_arch = "wasm32")]
pub fn run_web() -> myplotlib::WebResult {
    myplotlib::run_web(DFB_APP)
}
