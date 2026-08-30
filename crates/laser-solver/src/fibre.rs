use crate::dopant::{DopantModel, TwoLevelDopant};
use crate::grating::{GratingModel, NoGrating};

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
const SPEED_OF_LIGHT_MS: f64 = 299_792_458.0;

fn numerical_aperture(n_core: f64, n_cladding: f64) -> f64 {
    (n_core * n_core - n_cladding * n_cladding).sqrt()
}

fn v_number(numerical_aperture: f64, core_radius: f64, wavelength: f64) -> f64 {
    TWO_PI * core_radius * numerical_aperture / wavelength
}

fn dimensionless_marcuse_radius(v_number: f64) -> f64 {
    // this is an approximation for the ratio mode_radius / core_radius
    // D. Marcuse, “Loss analysis of single-mode fiber splices,” Bell System Technical Journal, vol. 56, no. 5, pp. 703–718, 1977.
    0.65 + 1.619 / v_number.powf(1.5) + 2.879 / v_number.powi(6)
}

#[derive(Copy, Clone)]
pub struct FieldMode {
    wavelength: f64,
}

impl FieldMode {
    pub const fn new(wavelength: f64) -> Self {
        Self { wavelength }
    }

    pub fn wavelength(self) -> f64 {
        self.wavelength
    }
}

impl Default for FieldMode {
    fn default() -> Self {
        Self::new(1.0)
    }
}

pub struct ActiveMode<D: DopantModel> {
    pub(crate) mode: FieldMode,
    pub(crate) overlap: f64,
    pub(crate) interaction: D::Interaction,
}

impl<D: DopantModel> ActiveMode<D> {
    pub(crate) fn new<G: GratingModel>(
        fibre: &Fibre<D, G>,
        mode: FieldMode,
        interaction: D::Interaction,
    ) -> Self {
        Self {
            mode,
            overlap: fibre.geometry.mode_overlap(mode),
            interaction,
        }
    }
}

impl<D> Clone for ActiveMode<D>
where
    D: DopantModel,
    D::Interaction: Clone,
{
    fn clone(&self) -> Self {
        Self {
            mode: self.mode,
            overlap: self.overlap,
            interaction: self.interaction.clone(),
        }
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct BidirectionalAmplitude {
    pub forward: f64,
    pub backward: f64,
}

impl BidirectionalAmplitude {
    pub fn forward_flux(self) -> f64 {
        self.forward * self.forward
    }

    pub fn backward_flux(self) -> f64 {
        self.backward * self.backward
    }

    pub fn total_flux(self) -> f64 {
        self.forward_flux() + self.backward_flux()
    }
}

pub(crate) fn bidirectional_amplitudes(total_flux: f64, balance: f64) -> (f64, f64) {
    let forward_fraction = (balance + 1.0) * 0.5;
    (
        (forward_fraction * total_flux).sqrt(),
        ((1.0 - forward_fraction) * total_flux).sqrt(),
    )
}

#[derive(Clone)]
pub struct FibreGeometry {
    pub core_radius: f64,
    pub numerical_aperture: f64,
    pub length: f64,
}

impl FibreGeometry {
    fn single_mode_cutoff_wavelength(&self) -> f64 {
        TWO_PI * self.core_radius * self.numerical_aperture / 2.405
    }

    fn is_single_spatial_mode(&self, mode: FieldMode) -> bool {
        mode.wavelength > self.single_mode_cutoff_wavelength()
    }

    fn v_number(&self, mode: FieldMode) -> f64 {
        TWO_PI * self.core_radius * self.numerical_aperture / mode.wavelength
    }

    pub(crate) fn mode_overlap(&self, mode: FieldMode) -> f64 {
        //let _ = mode;
        let v = self.v_number(mode);
        let mode_over_core = dimensionless_marcuse_radius(v);
        let gamma = 1.0 - f64::exp(-2.0 / (mode_over_core * mode_over_core));
        gamma //* 1e-10 + 1.0
        //1.0
    }
}

#[derive(Clone)]
pub struct Fibre<
    D: DopantModel = TwoLevelDopant,
    G: GratingModel = NoGrating,
> {
    pub geometry: FibreGeometry,
    pub dopant: D,
    pub grating: G,
}

impl<G: GratingModel + Default> Default for Fibre<TwoLevelDopant, G> {
    fn default() -> Self {
        Self {
            geometry: FibreGeometry {
                core_radius: 4e-6,
                numerical_aperture: 0.1,
                length: 5.0,
            },
            dopant: TwoLevelDopant {
                density: 1.0,
                lifetime: 1.0,
            },
            grating: G::default(),
        }
    }
}


impl BidirectionalAmplitude {
    pub fn coupled_step(self, gain: f64, kappa: f64, dz: f64) -> Self {
        let (a, b, c, d) = transfer(gain, kappa, dz);
        Self {
            forward: a * self.forward + b * self.backward,
            backward: c * self.forward + d * self.backward,
        }
    }

    pub fn uncoupled_step(self, gain: f64, dz: f64) -> Self {
        let factor = (0.5 * gain * dz).exp();
        Self {
            forward: self.forward * factor,
            backward: self.backward / factor,
        }
    }
}

pub fn transfer(gain: f64, kappa: f64, dz: f64) -> (f64, f64, f64, f64) {
    let g_dz = 0.5 * gain * dz;
    let k_dz = kappa * dz;
    let x = (g_dz * g_dz + k_dz * k_dz).sqrt();

    let cosh = x.cosh();
    let sinch = if x > 1e-30 { x.sinh() / x } else { 1.0_f64 };
    (
        cosh + g_dz * sinch,
        k_dz * sinch,
        k_dz * sinch,
        cosh - g_dz * sinch,
    )
}
