use num::Zero;
use num::complex::Complex;
use qdft::QDFT;

fn main() {
    let mut qdft = QDFT::new(44100.0, (100.0, 20000.0));

    let n = 1;
    let mut samples = vec![f32::zero(); n];
    let mut dfts = vec![Complex::<f64>::zero(); n];

    qdft.qdft(&samples, &mut dfts);
    qdft.iqdft(&dfts, &mut samples);
}
