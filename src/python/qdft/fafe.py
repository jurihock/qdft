"""
Copyright (c) 2023 Juergen Hock

SPDX-License-Identifier: MIT

Fast, Accurate Frequency Estimator according to [1].

[1] Eric Jacobsen and Peter Kootsookos
    Fast, Accurate Frequency Estimators
    IEEE Signal Processing Magazine (2007)
    https://ieeexplore.ieee.org/document/4205098

Source: https://github.com/jurihock/qdft
"""


import numpy


class FAFE:

    def __init__(self, mode=None):

        self.mode = mode

    def __call__(self, dfts):

        dfts = numpy.atleast_2d(dfts)
        assert dfts.ndim == 2

        l = numpy.roll(dfts, +1, axis=-1)
        m = dfts
        r = numpy.roll(dfts, -1, axis=-1)

        if self.mode is None:

            with numpy.errstate(all='ignore'):
                drifts = -numpy.real((r - l) / (2 * m - r - l))

        elif str(self.mode).lower() == 'p':

            p = 1.36  # TODO: hann?

            l = numpy.abs(l)
            m = numpy.abs(m)
            r = numpy.abs(r)

            with numpy.errstate(all='ignore'):
                drifts = p * (r - l) / (m + r + l)

        elif str(self.mode).lower() == 'q':

            q = 0.55  # TODO: hann?

            with numpy.errstate(all='ignore'):
                drifts = -numpy.real(q * (r - l) / (2 * m + r + l))

        else:

            drifts = numpy.zeros(dfts.shape)

        drifts[...,  0] = 0
        drifts[..., -1] = 0

        return drifts


class QFAFE:

    def __init__(self, qdft):

        self.qdft = qdft
        self.fafe = FAFE('p' if qdft.window is not None else None)

    def hz(self, dfts):

        oldfreqs = self.qdft.frequencies

        oldbins = numpy.arange(oldfreqs.size)
        newbins = oldbins + self.fafe(dfts)

        # TODO: is approximation possible? https://en.wikipedia.org/wiki/Cent_(music)
        newfreqs = self.qdft.bandwidth[0] * numpy.power(2, newbins / self.qdft.resolution)

        # TODO: does interp make sense?
        # newfreqs = numpy.interp(newbins, oldbins, oldfreqs)

        return newfreqs

    def cent(self, dfts):

        newfreqs = self.hz(dfts)
        oldfreqs = self.qdft.frequencies

        return 1200 * numpy.log2(newfreqs / oldfreqs)
