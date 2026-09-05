mod controls;
mod dfb;
mod kappa_plot;
mod pi_pos_plot;
mod pi_pos_threshold_plot;
mod pop_plot;
mod profile_plot;
mod threshold_plot;

#[cfg(not(target_arch = "wasm32"))]
pub use crate::dfb::run_native;
#[cfg(target_arch = "wasm32")]
pub use crate::dfb::run_web;
#[cfg(not(target_arch = "wasm32"))]
pub use myplotlib::NativeResult;
