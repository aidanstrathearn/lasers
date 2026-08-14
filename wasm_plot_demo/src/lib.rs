use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use wasm_bindgen::{JsCast as _, prelude::*};

#[derive(Default)]
struct PlotDemo;

impl eframe::App for PlotDemo {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("eframe + egui_plot + WebAssembly");

            let points = PlotPoints::from_iter((0..=200).map(|i| {
                let x = f64::from(i) * 0.05;
                [x, x.sin()]
            }));

            Plot::new("sine-wave")
                .x_axis_label("x")
                .y_axis_label("sin(x)")
                .show(ui, |plot_ui| {
                    plot_ui.line(Line::new("sine", points));
                });
        });
    }
}

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

    // This demo has one runner for the lifetime of the browser page.
    let runner = Box::leak(Box::new(eframe::WebRunner::new()));
    wasm_bindgen_futures::spawn_local(async move {
        runner
            .start(
                canvas,
                eframe::WebOptions::default(),
                Box::new(|creation_context| {
                    creation_context
                        .egui_ctx
                        .set_visuals(egui::Visuals::light());
                    Ok(Box::new(PlotDemo))
                }),
            )
            .await
            .expect("failed to start eframe");
    });

    Ok(())
}
