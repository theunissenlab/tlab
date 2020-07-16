'''
Implementation of multitaper coherence estimation per Hsu et. al.

"Quantifying variability in neural responses and its application
for the validation of model predictions" Anne Hsu, Alexander Borst,
and FrŽderic E Theunissen
Network: Computation in Neural Systems 15 (2004) 91-109
'''

from numpy import (arange, diff, empty, log2, max, mean, nonzero, ones_like,
                   sqrt, sum, std, zeros_like)
from scipy.stats import t
from matplotlib.pyplot import gca
from nitime import descriptors as desc
from nitime.timeseries import Epochs, TimeArray, TimeInterface
from nitime.analysis import BaseAnalyzer, MTCoherenceAnalyzer
from .rate import Raster, psth, to_float

class SpikeCoherenceAnalyzer(BaseAnalyzer):
    
    '''
    Base class for combined windowed-multitaper coherence analyzers.
    '''
    
    def __init__(self, input, bandwidth=None, alpha=.05,
                 win_length=1024, overlap=None, pad=True):
        '''
        input:        a tlab.src.spike.rate.Raster object
        bandwidth:    the bandwidth for the multitaper coherence, in units
                        specified by nitime.analysis.MTCoherenceAnalyzer
                        Units are a fraction of the nyquist rate: if bw = 2,
                        the bandwidth is half the nyquist rate or 1/4 of the
                        sampling rate.
        alpha:        probability for confidence interval calculation
        '''
        
        BaseAnalyzer.__init__(self,input)
        
        self.bandwidth = bandwidth
        self.alpha = alpha
        self.pad = pad
        
        if not isinstance(win_length,TimeInterface):
            win_length = TimeArray(win_length/self.input.sampling_rate)
        self.win_length = win_length
        
        if overlap is None:
            overlap = win_length/2.
        if not isinstance(overlap,TimeInterface):
            overlap = TimeArray(overlap/self.input.sampling_rate)
        self.overlap = overlap
            
    @property
    def coherence(self):
        '''
        Coherence values, limited to areas where lower confidence bound
        is positive.
        '''
        return self._coherence
    
    @property
    def coherence_upper_bound(self):
        return self._coherence_upper
    
    @property
    def coherence_lower_bound(self):
        return self._coherence_lower * self.good_bands
    
    @property
    def info(self):
        return coherence_info(self.coherence, self.max_good_freq,
                              self.freqs)
    
    @property
    def info_upper_bound(self):
        return coherence_info(self.coherence_upper_bound, self.max_good_freq,
                              self.freqs)
    
    @property
    def info_lower_bound(self):
        return coherence_info(self.coherence_lower_bound, self.max_good_freq,
                              self.freqs)
    
    @property
    def confidence_interval(self):
        '''
        Confidence interval, limited to areas where lower confidence bound
        is positive.
        '''
        return self._confidence_interval
    
    @desc.setattr_on_read
    def good_bands(self):
        positive = self._coherence_lower > 0
        first_neg = nonzero(-positive)[0][0]
        positive[first_neg:] = False
        return positive
    
    @desc.setattr_on_read
    def max_good_freq(self):
        return self.freqs[self.good_bands][-1]
    
    def plot(self,axes=None,max_freq=100):
        
        axes = axes or gca()
        
        # Frequency range to plot
        plot_freqs = self.freqs <= max_freq
        freqs = self.freqs[plot_freqs]
        
        axes.plot(freqs,self.coherence[plot_freqs],'k-',
                  freqs,self.coherence_upper_bound[plot_freqs],'r--',
                  freqs,self.coherence_lower_bound[plot_freqs],'r--')

class AutoCoherenceAnalyzer(SpikeCoherenceAnalyzer):
    '''
    Analyzer implementing computation of coherence between trials
    of multi-trial spike data
    '''
    
    @property
    def _coherence(self):
        return self.coherence_halves_mean
    
    @property
    def _coherence_upper(self):
        return self._coherence + self.confidence_interval
    
    @property
    def _coherence_lower(self):
        return self._coherence - self._confidence_interval

    @property
    def _confidence_interval(self):
        t_int = t._ppf(1-self.alpha/2., self.n_windows-2)
        return t_int * self.coherence_halves_var(ddof=1)
    
    @property
    def coherence_halves_mean(self):
        return mean(self._coherence_halves,1)
    
    def coherence_halves_var(self,**kwargs):
        return std(self._coherence_halves,1,**kwargs)
    
    @property
    def freqs(self):
        return self._coherence_halves_raw[2]
    
    @property
    def n_windows(self):
        return self._coherence_halves.shape[1]
    
    @desc.setattr_on_read
    def _coherence_halves(self):
        '''
        Compute coherence using Hsu et. al. equation 8
        n.b. the coherences returned by nitime's coherence analyzers are the
        squared coherence
        '''
        m = self.input.n_trials / 2.
        lhs = .5 * (-m + m * sqrt(1 / self._coherence_halves_raw[0]))  # Left hand side of Hsu et. al. equation 8
        return 1 / (1 + lhs)
    
    @desc.setattr_on_read
    def _coherence_halves_raw(self):
        return coherence_windowed(self.input.two_halves(),
                                  win_length=self.win_length,
                                  overlap=self.overlap,
                                  pad=self.pad)
        
class CrossCoherenceAnalyzer(SpikeCoherenceAnalyzer):
    '''
    Analyzer implementing computation of coherence between multi-trial
    spike data and some other signal -- either multi-trial spike data,
    or an absolute signal.
    '''
    
    def __init__(self, input, signal, **kwargs):
        SpikeCoherenceAnalyzer(self,input,**kwargs)
        self.signal = signal

class PredictionCoherenceAnalyzer(SpikeCoherenceAnalyzer):
    '''
    Full analyzer implementing both auto-coherence and cross-coherence
    '''
    
    pass

def coherence(raster,prediction):
    '''
    Compute coherence between two 
    '''
    pass

def coherence_info(coherence, max_freq, frequencies):
    '''
    Compute normal mutual information from coherence.
    
    coherence:    array of coherence values.
    mask:         boolean array indicating which 
    '''
    df = mean(diff(frequencies))
    mask = (frequencies > 0) & (frequencies <= max_freq)
    in_range = max([1-coherence[mask],zeros_like(mask)],axis=0)
    return -df * sum(log2(in_range))

def auto_coherence(raster,method='halves',**kwargs):
    '''
    Estimate coherence between a single spike trial in raster
    and the true mean rate.
    
    switch 'method' determines the method to use in calculation:
        'halves':       use the approximation of Hsu et. al equation (8)
        'jackknife':    jackknife over the trials (not implemented)
        'bootstrap':    use a bootstrap estimate (not implemented)
    '''
    
    if method == 'halves':
        return auto_coherence_halves(raster,**kwargs)
    elif method == 'jackknife':
        raise ValueError, 'jackknife auto coherence not implemented yet'
    elif method == 'bootstrap':
        raise ValueError, 'bootstrap auto coherence not implemented yet'
    else:
        raise ValueError, 'method "%s" not recognized' % method

def auto_coherence_halves(raster,**kwargs):
    '''
    Estimate coherence between a single spike trial in raster
    and the true mean rate using the approximation
    of Hsu et. al equation (8)
    '''

    halves = raster.two_halves()
    m = raster.n_trials/2.
    coh_obj = MTCoherenceAnalyzer(raster.two_halves(),**kwargs)
    lhs = .5 * (-m + m * sqrt(1 / coh_obj.coherence ** 2))  # Left hand side of Hsu et. al. equation 8
    gamma_ar = 1 / sqrt(1 + lhs)
    return gamma_ar, coh.confidence_interval, coh_obj.frequencies

def coherence_windowed(time_series,win_length=1024,overlap=None,pad=False,**kwargs):
    
    # Overlap defaults to win_length/2
    if overlap is None:
        overlap = win_length/2
    
    # Convert so that overlap is TimeArray, n_overlap is in sample points
    if isinstance(overlap,TimeInterface):
        n_overlap = to_float(overlap) * time_series.sampling_rate
    else:
        n_overlap = overlap
        overlap = TimeArray(float(n_overlap)/time_series.sampling_rate)
:        
    # Convert win_length to TimeArray
    if not isinstance(win_length,TimeInterface):
        win_length = TimeArray(float(win_length)/time_series.sampling_rate)
        
    # Compute number of windows
    n_time_points = time_series.data.shape[-1]
    (n_windows,rem) = divmod(n_time_points - n_overlap, n_overlap)
    
    # Adjust all parameters for padding or truncation
    n_windows += pad
    n_time_points = n_time_points - rem + pad * n_overlap
    
    # If 'pad' is true, pad the series with zeros out to one extra window.
    # Otherwise, truncate the last n_overlap points.
    if pad:
        time_series = time_series.padded(n_time_points)
    else:
        time_series = time_series.truncated(n_time_points)
        
    # Make a list of data windows
    epochs = Epochs(arange(n_windows) * overlap, duration = win_length)
    
    # Initialize outputs
    coh = MTCoherenceAnalyzer(time_series.during(epochs[0]))
    n_frequencies = len(coh.frequencies)
    all_coherences = empty((n_frequencies,n_windows))
    all_confidences = empty((n_frequencies,n_windows))
    
    # Loop over windows
    for idx,ep in enumerate(epochs[:-1]):
        slice = time_series.during(ep)
        coh.set_input(slice)
        all_coherences[:,idx] = coh.coherence[0,1].squeeze()
        all_confidences[:,idx] = coh.confidence_interval[0,1].squeeze()
        
    return all_coherences, all_confidences, coh.frequencies