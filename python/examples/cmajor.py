import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

from qdft import Chroma12

import matplotlib.pyplot as plot
import matplotlib.ticker as ticker
import numpy as np
import numpy.lib.stride_tricks as tricks
import wavfile as wav


def main():

    x, sr, bits = wav.read('cmajor.wav', fmt='float')

    x = np.array(x).mean(axis=-1)

    chroma = Chroma12(sr)

    chromagram = np.empty((0, chroma.size))

    chunks = tricks.sliding_window_view(x, sr)[::sr]

    for i, chunk in enumerate(chunks):

        percent = int(100 * (i + 1) / len(chunks))

        print(f'{percent}%')

        chromagram = np.vstack((chromagram, chroma.chroma(chunk)))

    args = dict(extent=(0, chromagram.shape[0] / sr, 0, chromagram.shape[1]),
                origin='lower',
                aspect='auto',
                cmap='inferno',
                interpolation='nearest')

    plot.imshow(chromagram.T, **args)
    cbar = plot.colorbar()

    plot.xlabel('s')
    cbar.set_label('CP')

    plot.gca().yaxis.set_major_locator(ticker.IndexLocator(1, 0))
    plot.gca().yaxis.set_major_formatter(lambda i, _: chroma.notes[int(i) % chroma.size])

    plot.tight_layout()

    plot.show()

    # plot.plot(np.arange(x.size) / sr, x)
    # plot.xlabel('s')
    # plot.ylim(-0.35, +0.35)
    # plot.tight_layout()
    # plot.show()


if __name__ == '__main__':

    main()
