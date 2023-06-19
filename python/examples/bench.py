import os, sys
src = os.path.join(os.path.dirname(__file__), '..', 'src')
sys.path.insert(0, src)

from qdft import QDFT
from timeit import default_timer as timer
import numpy as np

samplerate = 44100
bandwidth = (50, samplerate / 2)
resolution = 24
latency = 0
window = (+0.5, -0.5)

ta0 = timer()
qdft = QDFT(samplerate, bandwidth, resolution, latency, window)
tb0 = timer()
e0 = int((tb0 - ta0) * 1e+6)

print(f"PREP\tPYTHON\t{e0} us")

n = 1 * samplerate

x = np.ndarray(n, dtype=float)

runs = 10

for run in range(1, runs + 1):

    print(f"RUN\t{run}/{runs}")

    x[:] = 0

    ta1 = timer()
    y = qdft.qdft(x)
    tb1 = timer()
    e1 = int((tb1 - ta1) * 1e+6)

    ta2 = timer()
    x = qdft.iqdft(y)
    tb2 = timer()
    e2 = int((tb2 - ta2) * 1e+6)

    print(f"\tQDFT\t{e1} us")
    print(f"\tIQDFT\t{e2} us")
