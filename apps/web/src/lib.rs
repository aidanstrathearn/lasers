#[cfg(target_arch = "wasm32")]
use laser_app::LaserApp;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::{JsCast as _, prelude::*};

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn start() -> Result<(), JsValue> {
    let window = web_sys::window().ok_or_else(|| JsValue::from_str("missing window"))?;
    let document = window
        .document()
        .ok_or_else(|| JsValue::from_str("missing document"))?;
    let canvas = document
        .get_element_by_id("plot-canvas")
        .ok_or_else(|| JsValue::from_str("missing #plot-canvas"))?
        .dyn_into::<web_sys::HtmlCanvasElement>()?;
    let mobile_layout = window
        .match_media("(pointer: coarse)")?
        .is_some_and(|media_query| media_query.matches());

    // This demo has one runner for the lifetime of the browser page.
    let runner = Box::leak(Box::new(eframe::WebRunner::new()));
    wasm_bindgen_futures::spawn_local(async move {
        runner
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(move |creation_context| {
                    Ok(Box::new(LaserApp::new_with_mobile_layout(
                        creation_context,
                        mobile_layout,
                    )))
                }),
            )
            .await
            .expect("failed to start eframe");
    });

    Ok(())
}
