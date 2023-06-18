use num::Zero;
use num::complex::Complex;
use qdft::QDFT;

fn main() {
    let samplerate = 44100.0;
    let bandwidth = (50.0, samplerate / 2.0);
    let resolution = 24.0;
    let latency = 0.0;
    let window = Some((0.5, -0.5));

    let mut qdft = QDFT::new(
        samplerate,
        bandwidth,
        resolution,
        latency,
        window);

    let n = 1;
    let m = n * qdft.size();

    let mut samples = vec![f32::zero(); n];
    let mut dfts = vec![Complex::<f64>::zero(); m];

    qdft.qdft(&samples, &mut dfts);
    qdft.iqdft(&dfts, &mut samples);
}
