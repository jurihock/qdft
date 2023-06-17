#[cfg(test)]
mod tests {

    use num::One;
    use num::Zero;
    use num::complex::Complex;
    use qdft::QDFT;

    #[test]
    fn test() {
        let mut qdft = QDFT::new(44100.0, (100.0, 20000.0));

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
