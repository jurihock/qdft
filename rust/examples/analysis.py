import matplotlib.pyplot as plot
import numpy as np
import os

npy = os.path.join("examples", "analysis.npy")
dft = np.load(npy)

sr = 44100
n = len(dft)

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

# freqs = qdft.frequencies
# ticks = np.linspace(0, 1, freqs.size, endpoint=True)
# plot.gca().yaxis.set_major_formatter(lambda tick, _: np.interp(tick, ticks, freqs).astype(int))

plot.show()
