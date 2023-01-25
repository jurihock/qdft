# Constant-Q Sliding DFT in C++ and Python (QDFT)

![language](https://img.shields.io/badge/languages-C%2B%2B%20Python-blue)
![license](https://img.shields.io/github/license/jurihock/sdft?color=green)
![pypi](https://img.shields.io/pypi/v/qdft?color=gold)

Forward and inverse Constant-Q Sliding DFT according to [[1]](#1) with following features:

- Arbitrary octave resolution ([quarter tone](https://en.wikipedia.org/wiki/Quarter_tone) by default)
- Built-in parameterizable cosine family window (Hann by default)
- Customizable time and frequency domain data type in C++
- Endless single or multiple sample processing at once
- Optional analysis latency control parameter
- Real-time analysis and synthesis capability

The Constant-Q Sliding Discrete Fourier Transform (QDFT) is a recursive approach to compute the Fourier transform sample by sample. This is an efficient implementation without the FFT calculus. Just define an arbitrary frequency range and octave resolution to obtain the corresponding DFT estimate. In contrast to the linear [SDFT](https://github.com/jurihock/sdft), frequency bins of the QDFT are logarithmically spaced. Thus, both high and low frequencies are resolved with the same quality, which is particularly useful for audio analysis. Based on the QDFT, a chromagram feature with detailed instantaneous frequency estimation is planned for the future release.

## WIP

- [x] Readme
- [ ] Docstrings
- [x] PyPI package [qdft](https://pypi.org/project/qdft)
- [ ] Sliding [chromagram](https://en.wikipedia.org/wiki/Chroma_feature) as a bonus (a draft is already included in the Python package)

## Basic usage

### C++

```c++
#include <qdft/qdft.h> // see also src/cpp folder

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

### Python

```python
from qdft import QDFT # see also src/python folder

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
| ![SDFT](https://github.com/jurihock/qdft/raw/main/examples/face.png) | ![STFT](https://github.com/jurihock/qdft/raw/main/examples/cmajor.png) |
| [face.py](https://github.com/jurihock/qdft/blob/main/examples/face.py) | [cmajor.py](https://github.com/jurihock/qdft/blob/main/examples/cmajor.py) |
| ![SDFT](https://github.com/jurihock/qdft/raw/main/examples/face.wav.png) | ![STFT](https://github.com/jurihock/qdft/raw/main/examples/cmajor.wav.png) |

## See also

If you're interested in Sliding DFT with *linear* frequency resolution, don't forget to browse my [jurihock/sdft](https://github.com/jurihock/sdft) project!

## References

1. <span id="1">Russell Bradford et al. (2008). Sliding with a Constant Q. International Conference on Digital Audio Effects. https://www.dafx.de/paper-archive/2008/papers/dafx08_63.pdf</span>

2. <span id="2">Russell Bradford et al. (2005). Sliding is Smoother Than Jumping. International Computer Music Conference Proceedings. http://hdl.handle.net/2027/spo.bbp2372.2005.086</span>

3. <span id="3">Eric Jacobsen and Peter Kootsookos (2007). Fast, Accurate Frequency Estimators. IEEE Signal Processing Magazine. https://ieeexplore.ieee.org/document/4205098</span>

## License

[github.com/jurihock/qdft](https://github.com/jurihock/qdft) is licensed under the terms of the MIT license.
For details please refer to the accompanying [LICENSE](https://github.com/jurihock/qdft/raw/main/LICENSE) file distributed with it.
