"""
Copyright (c) 2023 Juergen Hock

SPDX-License-Identifier: MIT

Constant-Q Sliding DFT implementation according to [1] and [2].

[1] Russell Bradford and Richard Dobson and John ffitch
    Sliding with a Constant Q
    International Conference on Digital Audio Effects (2008)
    https://www.dafx.de/paper-archive/2008/papers/dafx08_63.pdf

[2] Benjamin Blankertz
    The Constant Q Transform
    https://doc.ml.tu-berlin.de/bbci/material/publications/Bla_constQ.pdf

Source: https://github.com/jurihock/qdft
"""


import numba
import numpy


class QDFT:
    """
    Constant-Q Sliding Discrete Fourier Transform (QDFT).
    """

    def __init__(self, samplerate, bandwidth, resolution=24, latency=0, window=(+0.5,-0.5)):
        """
        Create a new QDFT plan.

        Parameters
        ----------
        samplerate : int
            Sample rate in hertz.
        bandwidth : tuple(float, float)
            Lowest and highest frequency in hertz to be resolved.
        resolution : int, optional
            Octave resolution, e.g. number of bins per octave.
        latency : float, optional
            Analysis latency factor between -1 and +1.
        window : tuple(float, float), optional
            Cosine family window coeffs, e.g. (+0.5,-0.5) in case of hann window.
        """

        kernels = numpy.array([0, +1, -1] if window is not None else [0])

        quality = (2 ** (1 / resolution) - 1) ** -1
        size = numpy.ceil(resolution * numpy.log2(bandwidth[1] / bandwidth[0])).astype(int)
        frequencies = bandwidth[0] * numpy.power(2, numpy.arange(size) / resolution)

        periods = numpy.ceil(quality * samplerate / frequencies).astype(int)
        offsets = numpy.ceil((periods[0] - periods) * numpy.clip(latency * 0.5 + 0.5, 0, 1)).astype(int)
        weights = 1 / periods

        fiddles = numpy.exp(-2j * numpy.pi * (quality + kernels))
        twiddles = numpy.exp(+2j * numpy.pi * (quality + kernels[:, None]) / periods)

        inputs = numpy.zeros(periods[0], dtype=float)
        outputs = numpy.zeros((kernels.size, size), dtype=complex)

        self.samplerate = samplerate
        self.bandwidth = bandwidth
        self.resolution = resolution
        self.latency = latency
        self.window = window
        self.quality = quality
        self.size = size
        self.frequencies = frequencies
        self.periods = periods
        self.offsets = offsets
        self.weights = weights
        self.fiddles = fiddles
        self.twiddles = twiddles
        self.inputs = inputs
        self.outputs = outputs
        self.kernels = kernels

        dfts = numpy.empty((kernels.size, 0, size), dtype=complex)
        QDFT.accumulate(0, dfts, twiddles, outputs)

    def qdft(self, samples):
        """
        Estimate the DFT matrix for the given sample array.

        Parameters
        ----------
        samples : ndarray, list, float
            Array of samples.

        Returns
        -------
        dfts : ndarray
            DFT matrix of shape (samples,frequencies).
        """

        samples = numpy.atleast_1d(samples).astype(float)
        assert samples.ndim == 1
        assert samples.size >= 1

        inputs = numpy.concatenate((self.inputs, samples))
        self.inputs = inputs[samples.size:]

        outputs = self.outputs

        periods = self.periods
        offsets = self.offsets + numpy.arange(samples.size)[:, None]
        weights = self.weights

        fiddles = self.fiddles
        twiddles = self.twiddles

        window = self.window
        kernels = self.kernels

        dfts = (fiddles[:, None, None] * inputs[offsets + periods] - inputs[offsets]) * weights

        QDFT.accumulate(samples.size, dfts, twiddles, outputs)

        numpy.copyto(outputs, dfts[:, -1])

        if window is not None:

            a, b = window[0], window[1] / 2

            return a * dfts[0] + b * (dfts[-1] + dfts[+1])

        return dfts[0]

    def iqdft(self, dfts):
        """
        Synthesize the sample array from the given DFT matrix.

        Parameters
        ----------
        dfts : ndarray
            DFT matrix of shape (samples,frequencies).

        Returns
        -------
        samples : ndarray
            Array of samples.
        """

        dfts = numpy.atleast_2d(dfts).astype(complex)
        assert dfts.ndim == 2
        assert dfts.size >= 1

        samples = numpy.real(dfts * self.twiddles[0])

        return numpy.sum(samples, axis=1)

    @numba.njit()
    def accumulate(n, dfts, twiddles, identity):

        if not n: return

        dfts[:, 0] = twiddles * (identity + dfts[:, 0])

        for i in range(1, n):

            dfts[:, i] = twiddles * (dfts[:, i - 1] + dfts[:, i])
