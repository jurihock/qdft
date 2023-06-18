#[cfg(test)]
mod tests {

    use qdft::QDFT;

    use num::{One, Zero};

    #[allow(non_camel_case_types)]
    type c64 = num::complex::Complex<f64>;

    #[test]
    fn test_one_in_zero_out() {
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
        let m = qdft.size();

        let mut samples = vec![f32::zero(); n];
        let mut dfts = vec![c64::zero(); n * m];

        dfts[0] = c64::one();
        assert_ne!(dfts[0], c64::zero());

        qdft.qdft(&samples, &mut dfts);
        assert_eq!(dfts[0], c64::zero());

        samples[0] = f32::one();
        assert_ne!(samples[0], f32::zero());

        qdft.iqdft(&dfts, &mut samples);
        assert_eq!(samples[0], f32::zero());
    }
}
