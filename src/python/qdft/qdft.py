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
            Octave resolution, e.g. number of DFT bins per octave.
        latency : float, optional
            Analysis latency adjustment between -1 and +1.
        window : tuple(float, float), optional
            Cosine family window coeffs, e.g. (+0.5,-0.5) in case of hann window.
        """

        kernels = numpy.array([0, +1, -1]) \
                  if window is not None else \
                  numpy.array([0])

        window = numpy.array([window[0], window[+1]/2, window[-1]/2]) \
                 if window is not None else \
                 numpy.array([1])

        quality = (2 ** (1 / resolution) - 1) ** -1
        size = numpy.ceil(resolution * numpy.log2(bandwidth[1] / bandwidth[0])).astype(int)
        frequencies = bandwidth[0] * numpy.power(2, numpy.arange(size) / resolution)

        periods = numpy.ceil(quality * samplerate / frequencies).astype(int)
        offsets = numpy.ceil((periods[0] - periods) * numpy.clip(latency * 0.5 + 0.5, 0, 1)).astype(int)
        weights = 1 / periods

        fiddles = numpy.exp(-2j * numpy.pi * (quality + kernels))
        twiddles = numpy.exp(+2j * numpy.pi * (quality + kernels[:, None]) / periods)

        inputs = numpy.zeros(periods[0], dtype=float)
        outputs = numpy.zeros((size, kernels.size), dtype=complex)

        dfts = numpy.zeros((0, size, kernels.size), dtype=complex)
        QDFT.transform(dfts, inputs, outputs, periods, offsets, weights, fiddles, twiddles, window)

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
        outputs = self.outputs
        periods = self.periods
        offsets = self.offsets
        weights = self.weights
        fiddles = self.fiddles
        twiddles = self.twiddles
        window = self.window

        numpy.copyto(self.inputs, inputs[samples.size:])

        dfts = numpy.zeros((samples.size, self.size, self.window.size), complex)

        return QDFT.transform(dfts, inputs, outputs, periods, offsets, weights, fiddles, twiddles, window)

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
    def transform(dfts, inputs, outputs, periods, offsets, weights, fiddles, twiddles, window):

        if not dfts.size: return

        dfts[-1] = outputs

        if window.size == 3:

            for i in range(dfts.shape[0]):

                for j in range(dfts.shape[1]):

                    left = inputs[offsets[j] + periods[j] + i]
                    right = inputs[offsets[j] + i]

                    delta0 = (fiddles[0] * left - right) * weights[j]
                    delta1 = (fiddles[1] * left - right) * weights[j]
                    delta2 = (fiddles[2] * left - right) * weights[j]

                    dfts[i, j, 0] = twiddles[0, j] * (dfts[i - 1, j, 0] + delta0)
                    dfts[i, j, 1] = twiddles[1, j] * (dfts[i - 1, j, 1] + delta1)
                    dfts[i, j, 2] = twiddles[2, j] * (dfts[i - 1, j, 2] + delta2)

        else:

            for i in range(dfts.shape[0]):

                for j in range(dfts.shape[1]):

                    left = inputs[offsets[j] + periods[j] + i]
                    right = inputs[offsets[j] + i]

                    delta0 = (fiddles[0] * left - right) * weights[j]

                    dfts[i, j, 0] = twiddles[0, j] * (dfts[i - 1, j, 0] + delta0)

        outputs = dfts[-1]

        dfts *= window

        return dfts.sum(axis=-1)
