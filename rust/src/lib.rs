#![allow(unused)]

use num::Float;
use num::Zero;
use num::complex::Complex;
use std::marker::PhantomData;

pub struct QDFT<T, F>
    where T: Float,
          F: Float {
    samplerate: f64,
    bandwidth: (f64, f64),
    t: PhantomData<T>,
    f: PhantomData<F>,
}

pub type QDFT32 = QDFT<f32, f64>;
pub type QDFT64 = QDFT<f64, f64>;

impl<T, F> QDFT<T, F>
    where T: Float,
          F: Float {

    pub fn new(samplerate: f64,
               bandwidth: (f64, f64)
    ) -> Self {
        QDFT {
            samplerate: samplerate,
            bandwidth: bandwidth,
            t: PhantomData,
            f: PhantomData,
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
