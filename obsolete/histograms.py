'''
Created on 11.04.2013

@author: afedynitch
'''

import numpy as np

class Histogram(object):

    def __init__(self, unique_name, title, x_range, seldef, filldef,
                 n_bins=100, scale=None, log_scale=False):
        self.unique_name = unique_name
        self.title = title
        self.subset_name = unique_name
        self.x_range = x_range
        self.n_bins = n_bins
        if not scale:
            self.scale = '1.'
        else:
            self.scale = scale
        self.log_scale = log_scale
        self.n_events_filled = 0
        self.bufmax = int(1e4)
        self.buffil = 0
        self.buffer = np.zeros(self.bufmax)
        self.seldef = seldef
        self.filldef = filldef
        self._init_histogram()

    def __ne__(self, other):
        if (np.all(self.x_bins != other.x_bins)):
            return True
        else:
            return False

    def __eq__(self, other):
        return not (self != other)

    def __gt__(self, other):
        if (self.n_events_filled > other.n_events_filled):
            return True
        else:
            return False

    def __lt__(self, other):
        if (self.n_events_filled < other.n_events_filled):
            return True
        else:
            return False

    def _init_histogram(self):
        if self.x_range[1] - self.x_range[0] <= 0.0:
            raise Exception(
                "init_histograms(): You have specified a zero or negative bin range.")
        if not self.log_scale:
            self.x_bins, self.bin_width = np.linspace(self.x_range[0],
                                                      self.x_range[1],
                                                      self.n_bins + 1,
                                                      retstep=True)
            self.x_center_values = np.linspace(self.x_range[0] +
                                               0.5 * self.bin_width,
                                               self.x_range[1] -
                                               0.5 * self.bin_width, self.n_bins)
        else:
            self.x_bins = np.logspace(np.log10(self.x_range[0]),
                                      np.log10(self.x_range[1]),
                                      self.n_bins + 1)
            self.x_center_values = np.log10(0.5 *
                                            (10 ** self.x_bins[1:] +
                                             10 ** self.x_bins[:-1]))

        self.x_weights = self.x_bins[1:] - self.x_bins[:-1]
        self.run_hist_values = np.array(np.zeros(self.n_bins), dtype='d')
        self.run_sq_hist_values = np.array(np.zeros(self.n_bins), dtype='d')

    def __getstate__(self):
        result = self.__dict__.copy()
        try:
            del result['sel']
            del result['filler']
        except:
            pass
        try:
            del self.buffer
        except:
            pass
        try:
            del self.wbuffer
        except:
            pass

        return result

    def fill_event(self, event):
        self.n_events_filled += 1.
        
        try:
            event_series = self.filler(event, self.sel(event))
        except ValueError:
            print "Something weird happend."
            print event.eta, self.sel(event)
            print self.filler(event, self.sel(event))   
            self.n_events_filled -= 1
            return
        except:
            self.sel = eval(self.seldef)
            self.filler = eval(self.filldef)
            event_series = self.filler(event, self.sel(event))
        try:
            esize = event_series.size
        except AttributeError:
            esize = 1

        if self.buffil + esize < self.bufmax:
            self.buffer[self.buffil:self.buffil + esize] = event_series
            self.buffil += esize
        else:
            event_hist_values = np.histogram(self.buffer[:self.buffil],
                                             self.x_bins)[0]
            self.run_hist_values += event_hist_values
            self.run_sq_hist_values += event_hist_values ** 2
            self.buffer = np.zeros(self.bufmax)
            self.buffil = 0
            self.buffer[self.buffil:self.buffil + esize] = event_series
            self.buffil += esize

    def finalize_run(self, sigma_gen):
        if self.buffil > 0:
            event_hist_values = np.histogram(self.buffer[:self.buffil],
                                             self.x_bins)[0]
            self.run_hist_values += event_hist_values
            self.run_sq_hist_values += event_hist_values ** 2
        if self.n_events_filled < 1:
            print "Warning, trying to scale empty histogram", self.unique_name
            return
        variance = (self.run_sq_hist_values - self.run_hist_values ** 2 /
                    self.n_events_filled) / (self.n_events_filled - 1)
        self.run_hist_values = self.run_hist_values / self.x_weights / \
            self.n_events_filled
        self.errors = np.sqrt(variance) / self.x_weights / \
            self.n_events_filled
        self.sigma_gen = sigma_gen
        self.final_scale()

    def final_scale(self):
        self.scale = str(self.scale)
        if self.scale == 'norm':
            self.scale = 1. / np.sum(self.run_hist_values)
        elif self.scale == 'sigma':
            self.scale = self.sigma_gen
        else:
            self.scale = eval(self.scale)
        self.run_hist_values *= self.scale
        self.errors *= self.scale


class WeightedHistogram(Histogram):

    def __init__(self, *args, **kwargs):
        Histogram.__init__(self, *args, **kwargs)
        self.wbuffer = np.zeros(self.bufmax)
        self.n_vals_per_bin = np.zeros_like(self.run_hist_values)

    def fill_event(self, event):
        self.n_events_filled += 1.
        try:
            event_series, weights = self.filler(event, self.sel(event))
        except:
            self.sel = eval(self.seldef)
            self.filler = eval(self.filldef)
            event_series, weights = self.filler(event, self.sel(event))

        try:
            esize = event_series.size
        except AttributeError:
            esize = 1
        if np.any(np.isnan(weights)):
            return

        if self.buffil + esize < self.bufmax:
            self.buffer[self.buffil:self.buffil + esize] = event_series
            self.wbuffer[self.buffil:self.buffil + esize] = weights
            self.buffil += esize
        else:
            event_hist_values = np.histogram(self.buffer[:self.buffil],
                                             self.x_bins,
                                             weights=self.wbuffer[:self.buffil])[0]
            self.n_vals_per_bin += np.histogram(self.buffer[:self.buffil],
                                                self.x_bins)[0]
            self.run_hist_values += event_hist_values
            self.run_sq_hist_values += event_hist_values ** 2
            self.buffer = np.zeros(self.bufmax)
            self.wbuffer = np.zeros(self.bufmax)
            self.buffil = 0
            self.buffer[self.buffil:self.buffil + esize] = event_series
            self.wbuffer[self.buffil:self.buffil + esize] = weights
            self.buffil += esize

    def finalize_run(self, sigma_gen):
        if self.buffil > 0:
            event_hist_values = np.histogram(self.buffer[:self.buffil],
                                             self.x_bins,
                                             weights=self.wbuffer[:self.buffil])[0]
            self.n_vals_per_bin += np.histogram(self.buffer[:self.buffil],
                                                self.x_bins)[0]
            self.run_hist_values += event_hist_values
            self.run_sq_hist_values += event_hist_values ** 2

        if not (self.n_events_filled > 1):
            print "Warning, trying to scale empty histogram", self.unique_name
            return
        variance = (self.run_sq_hist_values - self.run_hist_values ** 2 /
                    self.n_events_filled) / (self.n_events_filled - 1)
        self.run_hist_values = self.run_hist_values / self.x_weights / \
            self.n_events_filled
        self.errors = np.sqrt(variance) / self.x_weights / \
            self.n_events_filled
        self.sigma_gen = sigma_gen
        self.final_scale()