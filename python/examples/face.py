import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

from qdft import QDFT

import matplotlib.pyplot as plot
import matplotlib.ticker as ticker
import numpy as np
import numpy.lib.stride_tricks as tricks
import wavfile as wav


def main():

    x, sr, bits = wav.read('face.wav', fmt='float')

    x = np.array(x).mean(axis=-1)

    qdft = QDFT(sr, bandwidth=(100, min(20e3, sr/2)), resolution=50, latency=1)

    spectrogram = np.empty((0, qdft.size))

    chunks = tricks.sliding_window_view(x, sr)[::sr]

    for i, chunk in enumerate(chunks):

        percent = int(100 * (i + 1) / len(chunks))

        print(f'{percent}%')

        dfts = qdft.qdft(chunk)

        with np.errstate(all='ignore'):
            db = 20 * np.log10(np.abs(dfts))

        spectrogram = np.vstack((spectrogram, db))

    args = dict(extent=(0, spectrogram.shape[0] / sr, 0, spectrogram.shape[1]),
                origin='lower',
                aspect='auto',
                cmap='inferno',
                interpolation='nearest')

    plot.imshow(spectrogram.T, **args)
    cbar = plot.colorbar()

    plot.xlabel('s')
    plot.ylabel('Hz')
    cbar.set_label('dB')

    plot.clim(-120, 0)

    plot.gca().yaxis.set_major_locator(ticker.IndexLocator(qdft.size//10, 0))
    plot.gca().yaxis.set_major_formatter(lambda i, _: round(qdft.frequencies[int(i) % qdft.size]))

    plot.tight_layout()

    plot.show()

    # plot.plot(np.arange(x.size) / sr, x)
    # plot.xlabel('s')
    # plot.ylim(-1.1, +1.1)
    # plot.tight_layout()
    # plot.show()


if __name__ == '__main__':

    main()
