use qdft::QDFT;

use num::Zero;
use std::time::Instant;

#[allow(non_camel_case_types)]
type c64 = num::complex::Complex<f64>;

fn main() {
    println!("RUST;\tQDFT;\tIQDFT");

    let samplerate = 44100.0;
    let bandwidth = (50.0, samplerate / 2.0);
    let resolution = 24.0;
    let quality = 0.0;
    let latency = 0.0;
    let window = Some((0.5, -0.5));

    let t0 = Instant::now();
    let mut qdft = QDFT::<f64, f64>::new(
        samplerate,
        bandwidth,
        resolution,
        quality,
        latency,
        window);
    let e0 = t0.elapsed().as_micros();

    println!("0;\t{e0};\t{e0}");

    let n = 1 * samplerate as usize;
    let m = qdft.size();

    let mut x = vec![f64::zero(); n];
    let mut y = vec![c64::zero(); n * m];

    let runs = 10;

    for run in 1 .. runs + 1 {
        x.fill(f64::zero());
        y.fill(c64::zero());

        let t1 = Instant::now();
        qdft.qdft(&x, &mut y);
        let e1 = t1.elapsed().as_micros();

        let t2 = Instant::now();
        qdft.iqdft(&y, &mut x);
        let e2 = t2.elapsed().as_micros();

        println!("{run};\t{e1};\t{e2}");
    }
}
