"""
Copyright (c) 2023 Juergen Hock

SPDX-License-Identifier: MIT

Source: https://github.com/jurihock/qdft
"""


import numpy

from .fafe import QFAFE
from .qdft import QDFT
from .scale import Scale


class Chroma:

    def __init__(self, samplerate, concertpitch=440, bandwidth=('A0', 'C#8'), latency=0, decibel=True, feature=None):

        scale = Scale(concertpitch)

        fmin = scale.frequency(bandwidth[0])
        fmax = scale.frequency(bandwidth[1])

        tmin = scale.semitone(fmin, 'abs')
        tmax = scale.semitone(fmax, 'abs')

        semitones = numpy.array(range(tmin, tmax))
        frequencies = numpy.array([scale.frequency(semitone) for semitone in semitones])
        notes = numpy.array([scale.note(semitone) for semitone in semitones])
        octaves = numpy.array([scale.octave(frequency) for frequency in frequencies])

        qdft = QDFT(samplerate, (fmin, fmax), 24, latency)
        assert numpy.allclose(qdft.frequencies[::2], frequencies)

        size = qdft.size // 2

        self.samplerate = samplerate
        self.concertpitch = concertpitch
        self.bandwidth = bandwidth
        self.latency = latency
        self.decibel = decibel
        self.feature = feature
        self.semitones = semitones
        self.frequencies = frequencies
        self.notes = notes
        self.octaves = octaves
        self.size = size
        self.scale = scale
        self.qdft = qdft

    def chroma(self, samples):

        dfts = self.qdft.qdft(samples)

        magnitudes = numpy.abs(dfts)
        features = None

        if self.decibel:

            with numpy.errstate(all='ignore'):
                magnitudes = 20 * numpy.log10(magnitudes)

        if str(self.feature).lower() in 'phase':

            features = numpy.angle(dfts)

        if str(self.feature).lower() in 'hz':

            fafe = QFAFE(self.qdft)
            features = fafe.hz(dfts)

        if str(self.feature).lower() in 'cent':

            fafe = QFAFE(self.qdft)
            features = fafe.cent(dfts)

        chromagram = (magnitudes + 1j * features) \
                     if features is not None \
                     else magnitudes

        chromagram = chromagram[..., ::2]
        assert chromagram.shape[-1] == self.frequencies.shape[-1]

        return chromagram


class Chroma12:

    def __init__(self, samplerate, concertpitch=440, bandwidth=('A0', 'C#8'), latency=0):

        superchroma = Chroma(samplerate, concertpitch, bandwidth, latency, decibel=False, feature=None)

        notes = superchroma.scale.scale
        size = len(notes)

        index = superchroma.scale.note(superchroma.notes[0])
        shift = -((size - index) % size)

        self.samplerate = samplerate
        self.concertpitch = concertpitch
        self.bandwidth = bandwidth
        self.latency = latency
        self.notes = notes
        self.size = size
        self.shift = shift
        self.superchroma = superchroma

    def chroma(self, samples):

        n = self.size

        src = self.superchroma.chroma(samples)
        dst = numpy.empty((len(src), n))

        for i in range(n):

            dst[:, i] = src[:, i::n].sum(axis=-1)

        dst = numpy.roll(dst, self.shift, axis=-1)

        return dst
