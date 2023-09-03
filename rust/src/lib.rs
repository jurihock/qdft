use num::{Float,Zero};
use num::complex::Complex;
use std::collections::VecDeque;

/// Cast between primitive types such as `f32` and `f64`.
pub trait CastFrom<T> {
    /// Casts the specified source value
    /// to the destination primitive type
    /// by using the `as` expression.
    fn cast(value: T) -> Self;
}

/// Implements the [`CastFrom`] trait for the specified
/// source `x` and destination `y` types.
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

// Since only f32 and f64 types are intended,
// the cast trait must only be implemented
// for these two types:
impl_cast_from_x_for_y!(f32, f32);
impl_cast_from_x_for_y!(f32, f64);
impl_cast_from_x_for_y!(f64, f32);
impl_cast_from_x_for_y!(f64, f64);

/// Constant-Q Sliding Discrete Fourier Transform (QDFT).
///
/// # Arguments
/// - `T` - Time domain data type.
/// - `F` - Frequency domain data type.
///
/// Both `T` and `F` must be `f32` or `f64` respectively.
pub struct QDFT<T, F>
    where T: Float + CastFrom<F>,
          F: Float + CastFrom<T> + CastFrom<f64> {
    samplerate: f64,
    bandwidth: (f64, f64),
    resolution: f64,
    quality: f64,
    latency: f64,
    window: Option<(f64, f64)>,
    size: usize,
    frequencies: Vec<f64>,
    qualities: Vec<f64>,
    latencies: Vec<f64>,
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

    /// Returns a new QDFT plan instance for the specified parameters.
    ///
    /// # Arguments
    /// - `samplerate` - Sample rate in hertz.
    /// - `bandwidth`  - Lowest and highest frequency in hertz to be resolved.
    /// - `resolution` - Octave resolution, e.g. number of DFT bins per octave.
    /// - `quality`    - Bandwidth offset for determining filter lengths.
    /// - `latency`    - Analysis latency adjustment between -1 and +1.
    /// - `window`     - Cosine family window coeffs, e.g. (+0.5,-0.5) in case of hann window.
    pub fn new(samplerate: f64,
               bandwidth: (f64, f64),
               resolution: f64,
               quality: f64,
               latency: f64,
               window: Option<(f64, f64)>,
    ) -> Self {
        let size = f64::ceil(resolution * f64::log2(bandwidth.1 / bandwidth.0)) as usize;

        let alpha = f64::powf(2.0, 1.0 / resolution) - 1.0;
        let beta = if quality < 0.0 { alpha * 24.7 / 0.108 } else { quality };

        let mut frequencies = vec![f64::zero(); size];
        let mut qualities = vec![f64::zero(); size];
        let mut latencies = vec![f64::zero(); size];
        let mut periods = vec![usize::zero(); size];
        let mut offsets = vec![usize::zero(); size];
        let mut weights = vec![F::zero(); size];

        for i in 0 .. size {
          let frequency = bandwidth.0 * f64::powf(2.0, (i as f64) / resolution);
          frequencies[i] = frequency;

          let quality = frequency / (alpha * frequency + beta);
          qualities[i] = quality;

          let period = f64::ceil(quality * samplerate / frequency);
          periods[i] = period as usize;

          let offset = f64::ceil(((periods[0] as f64) - period)
                     * f64::clamp(latency * 0.5 + 0.5, 0.0, 1.0));
          offsets[i] = offset as usize;

          let latency = (periods[0] as f64 - offset) / samplerate;
          latencies[i] = latency;

          let weight = F::one() / F::cast(period);
          weights[i] = weight;
        }

        let mut fiddles = vec![Complex::<F>::zero(); size * 3];
        let mut twiddles = vec![Complex::<F>::zero(); size * 3];

        for k in [-1, 0, 1] {
            let pi = std::f64::consts::PI; // acos(-1)

            let mut i = 0;
            let mut j = 1;

            while i < size {
                let quality = qualities[i];
                let period = periods[i] as f64;

                let fiddle = Complex::<F>::from_polar(
                    F::one(),
                    F::cast(-2.0 * pi * (quality + (k as f64))));
                fiddles[(j + k) as usize] = fiddle;

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
            quality,
            latency,
            window,
            size,
            frequencies,
            qualities,
            latencies,
            periods,
            offsets,
            weights,
            fiddles,
            twiddles,
            inputs,
            outputs,
        }
    }

    /// Returns the number of DFT bins derived from the `bandwidth` and `resolution`.
    pub fn size(&self) -> usize { self.size }

    /// Returns the sample rate in hertz.
    pub fn samplerate(&self) -> f64 { self.samplerate }

    /// Returns the lowest and highest frequency in hertz to be resolved.
    pub fn bandwidth(&self) -> (f64, f64) { self.bandwidth }

    /// Returns the octave resolution, e.g. number of DFT bins per octave.
    pub fn resolution(&self) -> f64 { self.resolution }

    /// Returns the bandwidth offset for determining filter lengths.
    pub fn quality(&self) -> f64 { self.quality }

    /// Returns the analysis latency factor.
    pub fn latency(&self) -> f64 { self.latency }

    /// Returns frequency values in hertz of the individual DFT bins.
    pub fn frequencies(&self) -> &[f64] { self.frequencies.as_slice() }

    /// Returns Q-factors (relative frequency resolution) of the individual DFT bins.
    pub fn qualities(&self) -> &[f64] { self.qualities.as_slice() }

    /// Returns latency values in seconds of the individual DFT bins.
    pub fn latencies(&self) -> &[f64] { self.latencies.as_slice() }

    /// Returns the cosine family window coeffs, e.g. (+0.5,-0.5) in case of hann window.
    pub fn window(&self) -> Option<(f64, f64)> { self.window }

    /// Estimate the DFT vector for the given sample.
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

                let k = (
                    (-1 + j) as usize,
                    ( 0 + j) as usize,
                    ( 1 + j) as usize,
                );

                let deltas = (
                    (fiddles[k.0] * left - right) * weight,
                    (fiddles[k.1] * left - right) * weight,
                    (fiddles[k.2] * left - right) * weight,
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

                let fiddle = self.fiddles[j];
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

    /// Synthesize the sample from the given DFT vector.
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

    /// Estimate the DFT matrix for the given sample array.
    #[inline]
    pub fn qdft_vector(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        assert_eq!(dfts.len(), samples.len() * self.size);
        for i in 0 .. samples.len() {
            let j = i * self.size .. (i + 1) * self.size;
            self.qdft_scalar(&samples[i], &mut dfts[j]);
        }
    }

    /// Synthesize the sample array from the given DFT matrix.
    #[inline]
    pub fn iqdft_vector(&mut self, dfts: &[Complex::<F>], samples: &mut [T]) {
        assert_eq!(dfts.len(), samples.len() * self.size);
        for i in 0 .. samples.len() {
            let j = i * self.size .. (i + 1) * self.size;
            self.iqdft_scalar(&dfts[j], &mut samples[i]);
        }
    }

    /// Estimate the DFT matrix for the given sample array.
    /// This is a shortcut for the function [`qdft_vector`].
    #[inline]
    pub fn qdft(&mut self, samples: &[T], dfts: &mut [Complex::<F>]) {
        self.qdft_vector(samples, dfts);
    }

    /// Synthesize the sample array from the given DFT matrix.
    /// This is a shortcut for the function [`iqdft_vector`].
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
        let quality = 0.0;
        let latency = 0.0;
        let window = Some((0.5, -0.5));

        let qdft = QDFT::<f32, f64>::new(
            samplerate,
            bandwidth,
            resolution,
            quality,
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
