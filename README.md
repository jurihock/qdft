# Constant-Q Sliding DFT in C++, Rust, and Python

![language](https://img.shields.io/badge/languages-C%2B%2B%20Rust%20Python-blue)
![license](https://img.shields.io/github/license/jurihock/sdft?color=green)
![pypi](https://img.shields.io/pypi/v/qdft?color=gold)
![creates](https://img.shields.io/crates/v/qdft?color=gold)

Forward and inverse Constant-Q Sliding DFT (QDFT) according to [[1]](#1) with following features:

- Arbitrary octave resolution ([quarter tone](https://en.wikipedia.org/wiki/Quarter_tone) by default)
- Built-in parameterizable cosine family window (Hann by default)
- Customizable time and frequency domain data type in C++
- Endless single or multiple sample processing at once
- Optional quality control parameter to smoothly reduce low frequency bandwidth and improve the time resolution
- Optional analysis latency control parameter
- Real-time analysis and synthesis capability

The Constant-Q Sliding Discrete Fourier Transform (QDFT) is a recursive approach to compute the Fourier transform sample by sample. This is an efficient implementation without the FFT calculus. Just define an arbitrary frequency range and octave resolution to obtain the corresponding DFT estimate. In contrast to the linear [SDFT](https://github.com/jurihock/sdft), frequency bins of the QDFT are logarithmically spaced. Thus, both high and low frequencies are resolved with the same quality, which is particularly useful for audio analysis. Based on the QDFT, a chromagram feature with detailed instantaneous frequency estimation is planned for the future release.

## WIP

- [x] Readme
- [ ] Docstrings
- [x] PyPI package [qdft](https://pypi.org/project/qdft)
- [x] Rust package [qdft](https://crates.io/crates/qdft)
- [ ] Sliding [chromagram](https://en.wikipedia.org/wiki/Chroma_feature) as a bonus (a draft is already included in the Python package)

## Basic usage

### C++

```c++
#include <qdft/qdft.h> // see also cpp folder

double sr = 44100;                             // sample rate in hertz
std::pair<double, double> bw = { 50, sr / 2 }; // lowest and highest frequency in hertz to be resolved
double r = 24;                                 // octave resolution, e.g. number of DFT bins per octave

QDFT<float, double> qdft(sr, bw, r); // create qdft plan for custom time and frequency domain data types

size_t n = ...;         // number of samples
size_t m = qdft.size(); // number of dft bins

float* x = ...; // analysis samples of shape (n)
float* y = ...; // synthesis samples of shape (n)

std::complex<double>* dft = ...; // dft matrix of shape (n, m)

qdft.qdft(n, x, dft);  // extract dft matrix from input samples
qdft.iqdft(n, dft, y); // synthesize output samples from dft matrix
```

The time domain data type defaults to `float` and the frequency domain data type to `double`.

### Rust

```rust
use qdft::QDFT; // see also rust folder

// just a shortcut for our complex number type
#[allow(non_camel_case_types)]
type c64 = num::complex::Complex<f64>;

// zero number trait, e.g. c64::zero()
use num::Zero;

let samplerate = 44100.0;                 // sample rate in hertz
let bandwidth = (50.0, samplerate / 2.0); // lowest and highest frequency in hertz to be resolved
let resolution = 24.0;                    // octave resolution, e.g. number of DFT bins per octave
let latency = 0.0;                        // analysis latency adjustment between -1 and +1
let window = Some((0.5, -0.5));           // hann window coeffs

// create qdft plan for custom time and frequency domain data types
let mut qdft = QDFT::<f32, f64>::new(
    samplerate,
    bandwidth,
    resolution,
    latency,
    window);

let n = ...;         // number of samples
let m = qdft.size(); // number of dft bins

let mut x = vec![f32::zero(); n]; // analysis samples of shape (n)
let mut y = vec![f32::zero(); n]; // synthesis samples of shape (n)

let mut dft = vec![c64::zero(); n * m]; // dft matrix of shape (n, m)

qdft.qdft(&x, &mut dft);  // extract dft matrix from input samples
qdft.iqdft(&dft, &mut y); // synthesize output samples from dft matrix
```

Alternatively use [ndarray](https://github.com/rust-ndarray/ndarray) instead of a flat array to allocate the DFT matrix, as shown in the `analysis.rs` example.

### Python

```python
from qdft import QDFT # see also python folder

sr = 44100        # sample rate in hertz
bw = (50, sr / 2) # lowest and highest frequency in hertz to be resolved
r = 24            # octave resolution, e.g. number of DFT bins per octave

qdft = QDFT(sr, bw, r) # create qdft plan

n = ...       # number of samples
m = qdft.size # number of dft bins (if need to know in advance)

x = ... # analysis samples of shape (n)

dft = qdft.qdft(x)  # extract dft matrix of shape (n, m) from input samples
y = qdft.iqdft(dft) # synthesize output samples from dft matrix
```

Feel free to obtain current version from [PyPI](https://pypi.org/project/qdft) by executing `pip install qdft`.

## Examples

| QDFT | Chroma12 |
| :--: | :------: |
| ![SDFT](https://github.com/jurihock/qdft/raw/main/python/examples/face.png) | ![STFT](https://github.com/jurihock/qdft/raw/main/python/examples/cmajor.png) |
| [face.py](https://github.com/jurihock/qdft/blob/main/python/examples/face.py) | [cmajor.py](https://github.com/jurihock/qdft/blob/main/python/examples/cmajor.py) |
| ![SDFT](https://github.com/jurihock/qdft/raw/main/python/examples/face.wav.png) | ![STFT](https://github.com/jurihock/qdft/raw/main/python/examples/cmajor.wav.png) |

## See also

If you're interested in Sliding DFT with *linear* frequency resolution, don't forget to browse my [jurihock/sdft](https://github.com/jurihock/sdft) project!

## References

1. <span id="1">Russell Bradford et al. (2008). Sliding with a Constant Q. International Conference on Digital Audio Effects. https://www.dafx.de/paper-archive/2008/papers/dafx08_63.pdf</span>

2. <span id="2">Russell Bradford et al. (2005). Sliding is Smoother Than Jumping. International Computer Music Conference Proceedings. http://hdl.handle.net/2027/spo.bbp2372.2005.086</span>

3. <span id="3">Eric Jacobsen and Peter Kootsookos (2007). Fast, Accurate Frequency Estimators. IEEE Signal Processing Magazine. https://ieeexplore.ieee.org/document/4205098</span>

## License

[github.com/jurihock/qdft](https://github.com/jurihock/qdft) is licensed under the terms of the MIT license.
For details please refer to the accompanying [LICENSE](https://github.com/jurihock/qdft/raw/main/LICENSE) file distributed with it.
