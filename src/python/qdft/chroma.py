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
        frequencies = numpy.array([scale.frequency(t) for t in semitones])
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

    def chroma(self, samples):

        dfts = self.qdft.qdft(samples)

        dfts = numpy.abs(dfts)

        if self.decibel:

            with numpy.errstate(all='ignore'):
                dfts = 20 * numpy.log10(dfts)

        dfts = dfts[..., ::2]
        assert dfts.shape[-1] == self.frequencies.shape[-1]

        return dfts
