#[cfg(not(target_arch = "wasm32"))]
fn main() -> myplotlib::NativeResult {
    dfb_example::run_native()
}

#[cfg(target_arch = "wasm32")]
fn main() {}
