use qdft::QDFT;

use ndarray::Array1 as Vector;
use ndarray::Array2 as Matrix;
use ndarray_npy::write_npy;

#[allow(non_camel_case_types)]
type c64 = num::complex::Complex<f64>;

fn phase(t: &Vector<f64>, f: &Vector<f64>) -> Vector<f64>
{
    let pi = std::f64::consts::PI;
    let mut cs = Vector::<f64>::zeros(t.len());

    for i in 1 .. t.len() {
        let dt = t[i] - t[i - 1];
        let dp = 2.0 * pi * f[i] * dt;
        cs[i] = dp + cs[i - 1];
    }

  return cs;
}

fn main() {
    let samplerate = 44100.0;
    let bandwidth = (50.0, samplerate / 2.0);
    let resolution = 24.0;
    let quality = 0.0;
    let latency = 0.0;
    let window = Some((0.5, -0.5));

    let mut qdft = QDFT::new(
        samplerate,
        bandwidth,
        resolution,
        quality,
        latency,
        window);

    let n = 1 * samplerate as usize;
    let m = qdft.size();

    let t = Vector::<f64>::range(0.0, n as f64, 1.0) / samplerate;
    let f = Vector::<f64>::ones(n) * 1000.0;

    let x = phase(&t, &f).mapv_into(|v|v.sin());
    let mut y = Matrix::<c64>::zeros((n, m));

    qdft.qdft(x.as_slice().unwrap(), y.as_slice_mut().unwrap());

    let npy = std::path::Path::new("examples").join("analysis.npy");
    let error = format!("Unable to create {}!", npy.to_str().unwrap());

    write_npy(npy, &y).expect(error.as_str());
}
