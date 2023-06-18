#[cfg(test)]
mod tests {

    use num::One;
    use num::Zero;
    use num::complex::Complex;
    use qdft::QDFT;

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
        let mut samples = vec![f32::zero(); n];
        let mut dfts = vec![Complex::<f64>::zero(); n];

        dfts[0] = Complex::one();
        assert_ne!(dfts[0], Complex::zero());

        qdft.qdft(&samples, &mut dfts);
        assert_eq!(dfts[0], Complex::zero());

        samples[0] = f32::one();
        assert_ne!(samples[0], f32::zero());

        qdft.iqdft(&dfts, &mut samples);
        assert_eq!(samples[0], f32::zero());
    }
}
