"""
Copyright (c) 2023 Juergen Hock

SPDX-License-Identifier: MIT

Source: https://github.com/jurihock/qdft
"""


import numpy

from .qdft import QDFT
from .scale import Scale


class Chroma:

    def __init__(self, samplerate, bandwidth=('A0', 'C#8'), concertpitch=440, decibel=True):

        scale = Scale(concertpitch)

        fmin = scale.frequency(bandwidth[0])
        fmax = scale.frequency(bandwidth[1])

        tmin = scale.semitone(fmin, 'abs')
        tmax = scale.semitone(fmax, 'abs')

        semitones = numpy.array(range(tmin, tmax))
        frequencies = numpy.array([scale.frequency(semitone) for semitone in semitones])
        notes = numpy.array([scale.note(semitone) for semitone in semitones])
        octaves = numpy.array([scale.octave(frequency) for frequency in frequencies])

        qdft = QDFT(samplerate, (fmin, fmax), 24)
        assert numpy.allclose(qdft.frequencies[::2], frequencies)

        size = qdft.size // 2

        self.samplerate = samplerate
        self.bandwidth = bandwidth
        self.concertpitch = concertpitch
        self.decibel = decibel
        self.semitones = semitones
        self.frequencies = frequencies
        self.notes = notes
        self.octaves = octaves
        self.size = size
        self.scale = scale
        self.qdft = qdft

    def analyze(self, dfts, mode=None):

        l = numpy.roll(dfts, +1, axis=-1)
        m = dfts
        r = numpy.roll(dfts, -1, axis=-1)

        if str(mode).lower() == 'p':

            p = 1.36

            l = numpy.abs(l)
            m = numpy.abs(m)
            r = numpy.abs(r)

            with numpy.errstate(all='ignore'):
                errors = p * (r - l) / (m + r + l)

        elif str(mode).lower() == 'q':

            q = 0.55

            with numpy.errstate(all='ignore'):
                errors = -numpy.real(q * (r - l) / (2 * m + r + l))

        else:

            with numpy.errstate(all='ignore'):
                errors = -numpy.real((r - l) / (2 * m - r - l))

        errors[...,  0] = 0
        errors[..., -1] = 0

        return errors

    def chroma(self, samples):

        dfts = self.qdft.qdft(samples)

        magnis = numpy.abs(dfts)
        errors = self.analyze(dfts)

        if self.decibel:

            with numpy.errstate(all='ignore'):
                magnis = 20 * numpy.log10(magnis)

        chromagram = magnis + 1j * errors

        chromagram = chromagram[..., ::2]
        assert chromagram.shape[-1] == self.frequencies.shape[-1]

        return chromagram
