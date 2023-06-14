import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

from qdft import Chroma, Scale

import matplotlib.pyplot as plot
import numpy as np


def phase(t, f):

    dt = np.diff(t, prepend=0)
    dp = 2 * np.pi * f * dt
    return np.cumsum(dp)


def main():

    # 1) build test signal containing two sinusoids

    sr = 44100  # sample rate in hertz
    cp = 440    # concert pitch in hertz

    s = 1  # test signal duration in seconds
    t = np.arange(s * sr) / sr  # discrete timeline

    scale = Scale(cp)

    xgains = [1, 1]
    xtones = ['A', 'B']
    xoctaves = [2] * 2
    xfreqs = [scale.frequency(tone, octave) for tone, octave in zip(xtones, xoctaves)]

    xcents = 0  # number of cents to shift second sinusoid frequency
    xfreqs += [xfreqs[-1] * 2 ** (xcents / 1200)]  # append shifted sinusoid

    x1 = xgains[0] * np.sin(phase(t, xfreqs[0]))    # build first sinusoid
    x2 = xgains[-1] * np.sin(phase(t, xfreqs[-1]))  # build second shifted sinusoid

    x = (x1 + x2) / 2  # mix both sinusoids to the test signal


    # 2) analyze previously built test signal

    chroma = Chroma(sr, cp, feature='cent')

    chromagram = chroma.chroma(x)[-1]  # pick last chromagram entry

    qmagnis = np.real(chromagram)  # log magnitudes
    qdrifts = np.imag(chromagram)  # frequency drifts in cents
    qfreqs = chroma.frequencies    # chromagram bin frequencies


    # 3) print expected vs. estimated sinusoid frequency drifts

    xbins = np.interp(xfreqs[:-1], qfreqs, np.arange(qfreqs.size)).astype(int)
    qcents = qdrifts[xbins]
    xcents = [0, xcents]

    print('Expected drift in cents', '\t', xcents)
    print('Estimated drift in cents', '\t', qcents)


    # 4) plot results

    figure, axes1 = plot.subplots()
    axes2 = axes1.twinx()

    a, = axes1.plot(qfreqs, qmagnis, label='chroma bin magnitudes')
    b, = axes2.plot(qfreqs, qdrifts, color='gray', label='frequency drift estimates')

    c = axes1.axvline(x=xfreqs[0], color='magenta', linestyle='dashed', label='sinusoid bin frequency 1')
    d = axes1.axvline(x=xfreqs[1], color='magenta', linestyle='dotted', label='sinusoid bin frequency 2')
    e = axes1.axvline(x=xfreqs[2], color='gray', linestyle='dotted', label='shifted sinusoid frequency 2')

    axes1.legend(handles=(a, b, c, d, e))

    axes1.set_ylabel('chroma bin magnitude / dB')
    axes2.set_ylabel('chroma bin drift / cent')

    axes1.set_xlabel('chroma bin frequency / Hz')
    axes2.set_xscale('log')

    plot.show()


if __name__ == '__main__':

    main()
