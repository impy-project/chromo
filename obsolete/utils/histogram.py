import numpy as np

class Histogram():
    def __init__(self, title, x_bins):
        self.title = title

        self.x_bins = x_bins
        self.n_bins = len(x_bins) - 1
        self.widths = self.x_bins[1:] - self.x_bins[:-1]

        self.values = np.array(np.zeros(self.n_bins), dtype='d')
        self.n_entries = 0

        # Buffer 10^5 entries before running np.histogram
        self.s_buf = int(1e5)
        self.buffer = np.zeros(self.s_buf, dtype='d')
        self.last_idx = 0
        self.n_events_filled = 0


    def fill_event(self, event_series):
        s_entr = event_series.size
        if self.last_idx + s_entr < self.s_buf:
            self.buffer[self.last_idx:self.last_idx + s_entr] = event_series
            self.last_idx = self.last_idx + s_entr
        else:
            self.values += np.histogram(self.buffer[:self.last_idx],
                                        self.x_bins)[0]
            self.last_idx = 0

        self.n_entries += 1.

    def flush(self):
        if self.last_idx > 0:
            try:
                self.values += np.histogram(self.buffer[:self.last_idx],
                                        self.x_bins)[0]
            except:
                pass
        self.last_idx = 0

    def get(self):
        self.flush()
        return self.values / self.widths / float(self.n_events_filled)

    def __iadd__(self, other):
        if not bool(np.alltrue(self.x_bins == other.x_bins)):
            raise Exception("Histogram():: Trying to add histograms with different" +
                            " binning.")

        self.values += other.values
        self.n_entries += other.n_entries
        self.n_events_filled += other.n_events_filled
        return self

    def __getstate__(self):
        try:
            self.flush()
            del self.buffer
        except TypeError:
            pass
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state
        self.buffer = None


class EnergySpectrum(Histogram):
    def __init__(self, part_id, pdg_id, E_lab,
        grid_var, **kwargs):
        Histogram.__init__(self, **kwargs)

        self.particle_id = part_id
        self.pdg_id = pdg_id
        self.E_lab = E_lab
        if grid_var == 'E':
            self.scale = 1.
        elif grid_var == 'x':
            self.scale = self.E_lab

    def fill_event(self, event):
        event_series = event.E[np.where(
            event.p_ids == self.particle_id)]/self.scale
        Histogram.fill_event(self, event_series)


class EnergySpectrumBoost(EnergySpectrum):
    def __init__(self, lorentz_factors, **kwargs):
        self.gamma, self.beta_gamma = lorentz_factors
        EnergySpectrum.__init__(self, **kwargs)

    def fill_event(self, event):
        sel = np.where(event.p_ids == self.particle_id)
        event_series = (self.gamma * event.E[sel] +
            self.beta_gamma * event.pz[sel])/self.scale
        Histogram.fill_event(self, event_series)