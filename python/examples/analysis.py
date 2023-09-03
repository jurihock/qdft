import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

from qdft import QDFT

import matplotlib.patches as patches
import matplotlib.pyplot as plot
import numpy as np


def phase(t, f):

    dt = np.diff(t, prepend=0)
    dp = 2 * np.pi * f * dt
    return np.cumsum(dp)


def main():

    # 1) generate input signal

    sr = 44100  # sample rate in Hz
    n = 1 * sr  # number of samples

    t = np.arange(n) / sr  # timestamps in seconds

    f = 1000  # single frequency
    # f = np.linspace(0, 1000, n)  # linear chirp
    # f = np.sin(np.pi * t * t[-1]) * 1000 + 1000  # frequency wave

    x = np.sin(phase(t, f))  # sample vector of shape (n)


    # 2) estimate output dft

    qdft = QDFT(sr, (50, sr/2))  # create qdft plan

    dft = qdft.qdft(x)  # get dft matrix of shape (n, m), where m=qdft.size


    # 3) plot spectrogram

    with np.errstate(all='ignore'):
        db = 20 * np.log10(np.abs(dft))

    roi = (0, n / sr, 0, 1)

    args = dict(extent=roi,
                origin='lower',
                aspect='auto',
                cmap='inferno',
                interpolation='nearest')

    plot.imshow(db.T, **args)
    cbar = plot.colorbar()

    plot.xlabel('s')
    plot.ylabel('Hz')
    cbar.set_label('dB')

    plot.clim(-120, 0)

    ticks = np.linspace(0, 1, qdft.size, endpoint=True)
    freqs = qdft.frequencies
    lates = list(zip(qdft.latencies, ticks))

    plot.gca().yaxis.set_major_formatter(lambda tick, _: np.interp(tick, ticks, freqs).astype(int))
    plot.gca().add_patch(patches.Polygon(lates[::+1] + lates[::-1], linewidth=1, linestyle='--', color='w', alpha=0.5))

    plot.show()


if __name__ == '__main__':

    main()
