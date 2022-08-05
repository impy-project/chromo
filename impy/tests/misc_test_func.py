"""Some common functions for tests"""


def plot_hist(xedges, ws, axes=None, facecolor=None, **kwargs):
    """
    Plots histogram data in ROOT style.

    Parameters
    ----------
    xedge: lower bin boundaries + upper boundary of last bin
    w: content of the bins
    """

    if axes is None:
        from matplotlib import pyplot as plt

        axes = plt.gca()

    import numpy as np

    m = len(ws)
    n = 2 * m + 2

    xs = np.zeros(n)
    ys = np.zeros(n)

    xs[0] = xedges[0]
    xs[-1] = xedges[-1]

    for i in xrange(m):
        xs[1 + 2 * i] = xedges[i]
        ys[1 + 2 * i] = ws[i]
        xs[1 + 2 * i + 1] = xedges[i + 1]
        ys[1 + 2 * i + 1] = ws[i]

    if not facecolor is None:
        return axes.fill(xs, ys, facecolor=facecolor, **kwargs)
    else:
        return axes.plot(xs, ys, **kwargs)
