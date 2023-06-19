#![allow(unused)]

use num::{Float, Zero};
use num::complex::Complex;
use std::collections::VecDeque;

pub trait CastFrom<T> {
    fn cast(value: T) -> Self;
}

macro_rules! impl_cast_from_x_for_y {
    ($x:ty, $y:ty) => {
        impl CastFrom<$x> for $y {
            #[inline]
            fn cast(value: $x) -> $y {
                value as $y
            }
        }
    }
}

impl_cast_from_x_for_y!(f32, f32);
impl_cast_from_x_for_y!(f32, f64);
impl_cast_from_x_for_y!(f64, f32);
impl_cast_from_x_for_y!(f64, f64);

pub struct QDFT<T, F>
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {
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
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {
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

        for i in 0 .. size {
          let frequency = bandwidth.0 * f64::powf(2.0, (i as f64) / resolution);
          frequencies[i] = frequency;

          let period = f64::ceil(quality * samplerate / frequency);
          periods[i] = period as usize;

          let offset = f64::ceil(((periods[0] as f64) - period)
                     * f64::clamp(latency * 0.5 + 0.5, 0.0, 1.0));
          offsets[i] = offset as usize;

          let weight = F::one() / F::cast(period);
          weights[i] = weight;
        }

        let mut fiddles = vec![Complex::<F>::zero(); 3];
        let mut twiddles = vec![Complex::<F>::zero(); size * 3];

        for k in [-1, 0, 1] {
            let pi = std::f64::consts::PI; // acos(-1)

            let fiddle = Complex::<F>::from_polar(
                F::one(),
                F::cast(-2.0 * pi * (quality + (k as f64))));
            fiddles[(k + 1) as usize] = fiddle;

            let mut i = 0;
            let mut j = 1;

            while i < size {
                let period = periods[i] as f64;

                let twiddle = Complex::<F>::from_polar(
                    F::one(),
                    F::cast(2.0 * pi * (quality + (k as f64)) / period));
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

    pub fn size(&self) -> usize { self.size }

    pub fn qdft_scalar(&mut self, sample: &T, dft: &mut [Complex::<F>]) {
        assert_eq!(dft.len(), self.size);

        let inputs = &mut self.inputs;
        let outputs = &mut self.outputs;

        inputs.pop_front();
        inputs.push_back(*sample);

        if self.window.is_some() {
            let w = self.window.unwrap();
            let a = F::cast(w.0);
            let b = F::cast(w.1 / 2.0);

            let mut i = 0;
            let mut j = 1;

            while i < self.size {
                let period = self.periods[i];
                let offset = self.offsets[i];
                let weight = self.weights[i];

                let fiddles = &self.fiddles;
                let twiddles = &self.twiddles;

                let left = F::cast(inputs[offset + period]);
                let right = F::cast(inputs[offset]);

                let deltas = (
                    (fiddles[0] * left - right) * weight,
                    (fiddles[1] * left - right) * weight,
                    (fiddles[2] * left - right) * weight,
                );

                let k = (
                    (-1 + j) as usize,
                    ( 0 + j) as usize,
                    ( 1 + j) as usize,
                );

                outputs[k.0] = twiddles[k.0] * (outputs[k.0] + deltas.0);
                outputs[k.1] = twiddles[k.1] * (outputs[k.1] + deltas.1);
                outputs[k.2] = twiddles[k.2] * (outputs[k.2] + deltas.2);

                dft[i] = outputs[k.1] * a + (outputs[k.0] + outputs[k.2]) * b;

                i += 1;
                j += 3;
            }
        }
        else {
            let mut i = 0;
            let mut j = 1;

            while i < self.size {
                let period = self.periods[i];
                let offset = self.offsets[i];
                let weight = self.weights[i];

                let fiddle = self.fiddles[1];
                let twiddle = self.twiddles[j];

                let left = F::cast(inputs[offset + period]);
                let right = F::cast(inputs[offset]);

                let delta = (fiddle * left - right) * weight;

                outputs[j] = twiddle * (outputs[j] + delta);

                dft[i] = outputs[j];

                i += 1;
                j += 3;
            }
        }
    }

    pub fn iqdft_scalar(&mut self, dft: &[Complex::<F>], sample: &mut T) {
        assert_eq!(dft.len(), self.size);

        let mut result = F::zero();

        let mut i = 0;
        let mut j = 1;

        while i < self.size {
            let twiddle = &self.twiddles[j];

            result = result + (dft[i] * twiddle).re;

            i += 1;
            j += 3;
        }

        *sample = T::cast(result);
    }

    #[inline]
    pub fn qdft_vector(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        assert_eq!(dfts.len(), samples.len() * self.size);
        for i in 0 .. samples.len() {
            let j = i * self.size .. (i + 1) * self.size;
            self.qdft_scalar(&samples[i], &mut dfts[j]);
        }
    }

    #[inline]
    pub fn iqdft_vector(&mut self, dfts: &[Complex::<F>], samples: &mut [T]) {
        assert_eq!(dfts.len(), samples.len() * self.size);
        for i in 0 .. samples.len() {
            let j = i * self.size .. (i + 1) * self.size;
            self.iqdft_scalar(&dfts[j], &mut samples[i]);
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
        assert_eq!(j, (qdft.size * 3 + 1) as isize);
    }
}
