#![allow(unused)]

use num::Float;
use num::Zero;
use num::complex::Complex;
use std::collections::VecDeque;

pub struct QDFT<T, F>
    where T: Float,
          F: Float {
    samplerate: f64,
    bandwidth: (f64, f64),
    resolution: f64,
    latency: f64,
    quality: f64,
    size: usize,
    window: Option<(f64, f64)>,
    frequencies: Vec<f64>,
    periods: Vec<usize>,
    offsets: Vec<usize>,
    weights: Vec<F>,
    fiddles: Vec<Complex<F>>,
    twiddles: Vec<Complex<F>>,
    inputs: VecDeque<T>,
    outputs: Vec<Complex<F>>,
}

pub type QDFT32 = QDFT<f32, f64>;
pub type QDFT64 = QDFT<f64, f64>;

impl<T, F> QDFT<T, F>
    where T: Float,
          F: Float {
    pub fn new(samplerate: f64,
               bandwidth: (f64, f64),
               resolution: f64,
               latency: f64,
               window: Option<(f64, f64)>,
    ) -> Self {
        let quality = f64::powf(f64::powf(2.0, 1.0 / resolution) - 1.0, -1.0);
        let size = f64::ceil(resolution * f64::log2(bandwidth.1 / bandwidth.0)) as usize;

        let mut frequencies = vec![f64::zero(); size];
        let mut periods = vec![usize::zero(); size];
        let mut offsets = vec![usize::zero(); size];
        let mut weights = vec![F::zero(); size];

        for i in 0..size {
          let frequency = bandwidth.0 * f64::powf(2.0, (i as f64) / resolution);
          frequencies[i] = frequency;

          let period = f64::ceil(quality * samplerate / frequency);
          periods[i] = period as usize;

          let offset = f64::ceil(((periods[0] as f64) - period)
                     * f64::clamp(latency * 0.5 + 0.5, 0.0, 1.0));
          offsets[i] = offset as usize;

          let weight = F::one() / F::from(period).unwrap();
          weights[i] = weight;
        }

        let mut fiddles = vec![Complex::<F>::zero(); 3];
        let mut twiddles = vec![Complex::<F>::zero(); size * 3];

        for k in [-1, 0, 1] {
            let pi = std::f64::consts::PI; // acos(-1)

            let fiddle = Complex::<F>::from_polar(F::one(), F::from(
                -2.0 * pi * (quality + (k as f64))).unwrap());
            fiddles[(k + 1) as usize] = fiddle;

            let mut i = 0;
            let mut j = 1;

            while i < size {
                let period = periods[i] as f64;

                let twiddle = Complex::<F>::from_polar(F::one(), F::from(
                    2.0 * pi * (quality + (k as f64)) / period).unwrap());
                twiddles[(j + k) as usize] = twiddle;

                i += 1;
                j += 3;
            }
        }

        let inputs = VecDeque::from(vec![T::zero(); periods[0] + 1]);
        let outputs = vec![Complex::<F>::zero(); size * 3];

        QDFT {
            samplerate,
            bandwidth,
            resolution,
            latency,
            quality,
            size,
            window,
            frequencies,
            periods,
            offsets,
            weights,
            fiddles,
            twiddles,
            inputs,
            outputs,
        }
    }

    pub fn qdft_scalar(&mut self, sample: &T, dft: &mut Complex::<F>) {
        *dft = Complex::<F>::zero();
    }

    pub fn iqdft_scalar(&mut self, dft: &Complex::<F>, sample: &mut T) {
        *sample = T::zero();
    }

    #[inline]
    pub fn qdft_vector(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        for i in 0..samples.len() {
            self.qdft_scalar(&samples[i], &mut dfts[i]);
        }
    }

    #[inline]
    pub fn iqdft_vector(&mut self, dfts: &[Complex::<F>], samples: &mut [T]) {
        for i in 0..samples.len() {
            self.iqdft_scalar(&dfts[i], &mut samples[i]);
        }
    }

    #[inline]
    pub fn qdft(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        self.qdft_vector(samples, dfts);
    }

    #[inline]
    pub fn iqdft(&mut self, dfts: &[Complex::<F>], samples: &mut [T]) {
        self.iqdft_vector(dfts, samples);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sanity_checks() {
        let samplerate = 44100.0;
        let bandwidth = (50.0, samplerate / 2.0);
        let resolution = 24.0;
        let latency = 0.0;
        let window = Some((0.5, -0.5));

        let qdft = QDFT::<f32, f64>::new(
            samplerate,
            bandwidth,
            resolution,
            latency,
            window);

        let mut i = 0;
        let mut j = 1;

        while i < qdft.size {
            assert!((qdft.offsets[i]) < qdft.inputs.len());
            assert!((qdft.offsets[i] + qdft.periods[i]) < qdft.inputs.len());

            assert!((-1 + 1) < qdft.fiddles.len() as isize);
            assert!(( 0 + 1) < qdft.fiddles.len() as isize);
            assert!(( 1 + 1) < qdft.fiddles.len() as isize);

            assert!((-1 + j) < qdft.twiddles.len() as isize);
            assert!(( 0 + j) < qdft.twiddles.len() as isize);
            assert!(( 1 + j) < qdft.twiddles.len() as isize);

            assert!((-1 + j) < qdft.outputs.len() as isize);
            assert!(( 0 + j) < qdft.outputs.len() as isize);
            assert!(( 1 + j) < qdft.outputs.len() as isize);

            i += 1;
            j += 3;
        }

        assert_eq!(i, qdft.size);
        assert_eq!(j, (qdft.size as isize) * 3 + 1);
    }
}
