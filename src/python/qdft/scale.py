"""
Copyright (c) 2023 Juergen Hock

SPDX-License-Identifier: MIT

Source: https://github.com/jurihock/qdft
"""


import numpy
import re


class Scale:

    def __init__(self, concertpitch=440):

        self.concertpitch = concertpitch
        self.c0 = concertpitch * 2 ** (-(9 + 4*12) / 12)
        self.scale = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    def note(self, semitone):

        assert isinstance(semitone, (numpy.integer, int, str))

        if isinstance(semitone, (numpy.integer, int)):
            return self.scale[semitone % 12]

        if isinstance(semitone, str):
            semitone = re.sub(r'\d+', '', semitone)
            return self.scale.index(semitone.upper())

    def frequency(self, semitone, octave=None):

        assert isinstance(semitone, (numpy.integer, numpy.floating, int, float, str))
        assert isinstance(octave, (numpy.integer, numpy.floating, int, float, type(None)))

        if isinstance(semitone, str):
            numbers = ''.join(re.findall(r'\d+', semitone))
            semitone = self.note(semitone)
            octave = int(numbers) if numbers else octave

        if octave is not None:
            semitone = semitone % 12

        return 2 ** (semitone / 12 + (octave or 0)) * self.c0

    def semitone(self, frequency, mode='relative'):

        assert isinstance(frequency, (numpy.integer, numpy.floating, int, float))
        assert isinstance(mode, str)

        if mode.lower() in 'relative':
            return round(12 * numpy.log2(frequency / self.c0)) % 12

        if mode.lower() in 'absolute':
            return round(12 * numpy.log2(frequency / self.c0))

        raise ValueError(f'Invalid mode "{mode}", must be "relative" or "absolute"!')

    def octave(self, frequency):

        assert isinstance(frequency, (numpy.integer, numpy.floating, int, float))

        return round(12 * numpy.log2(frequency / self.c0)) // 12
