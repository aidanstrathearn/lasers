mod controls;
mod dfb;
mod kappa_plot;
mod pi_pos_plot;
mod pi_pos_threshold_plot;
mod pop_plot;
mod profile_plot;
mod threshold_plot;

#[cfg(not(target_arch = "wasm32"))]
pub fn run_native() -> myplotlib::NativeResult {
    myplotlib::run_native(dfb::definition())
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen::prelude::wasm_bindgen(js_name = mountApp)]
pub async fn mount_app(
    canvas: web_sys::HtmlCanvasElement,
) -> Result<myplotlib::WebHandle, wasm_bindgen::JsValue> {
    myplotlib::mount_web(canvas, dfb::definition()).await
}
