# Constant-Q Sliding DFT in C++ and Python (QDFT)

![language](https://img.shields.io/badge/languages-C%2B%2B%20Python-blue)
![license](https://img.shields.io/github/license/jurihock/sdft?color=green)

Forward and inverse Constant-Q Sliding DFT according to [[1]](#1) with following features:

- Arbitrary octave resolution ([quarter tone](https://en.wikipedia.org/wiki/Quarter_tone) by default)
- Built-in parameterizable cosine family window (Hann by default)
- Customizable time and frequency domain data type in C++
- Endless single or multiple sample processing at once
- Optional analysis latency control parameter
- Real-time analysis and synthesis capability

## WIP

- [ ] Readme
- [ ] Docstrings
- [ ] PyPI package
- [ ] Sliding [chromagram](https://en.wikipedia.org/wiki/Chroma_feature) as a bonus

## See also

If you're interested in Sliding DFT with *linear* frequency resolution, don't forget to browse my [jurihock/sdft](https://github.com/jurihock/sdft) project!

## References

1. <span id="1">Russell Bradford et al. (2008). Sliding with a Constant Q. International Conference on Digital Audio Effects. https://www.dafx.de/paper-archive/2008/papers/dafx08_63.pdf</span>

2. <span id="2">Russell Bradford et al. (2005). Sliding is Smoother Than Jumping. International Computer Music Conference Proceedings. http://hdl.handle.net/2027/spo.bbp2372.2005.086</span>

3. <span id="1">Krzysztof Duda (2010). Accurate, Guaranteed Stable, Sliding Discrete Fourier Transform. IEEE Signal Processing Magazine. https://ieeexplore.ieee.org/document/5563098</span>

## License

[github.com/jurihock/qdft](https://github.com/jurihock/qdft) is licensed under the terms of the MIT license.
For details please refer to the accompanying [LICENSE](https://github.com/jurihock/qdft/raw/main/LICENSE) file distributed with it.
