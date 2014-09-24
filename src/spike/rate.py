# This Python file uses the following encoding: utf-8
from copy import copy
from numpy import (asarray, arange, bool_, concatenate, convolve, dot,
                   histogram, hstack, mean, NaN, repeat, vstack, zeros)
from matplotlib.pyplot import (imshow, cm, xlim, plot, draw_if_interactive,
                               gca, setp)
from matplotlib.ticker import MultipleLocator
from math import ceil
from nitime.timeseries import TimeSeries
from ..stats import cc as sequence_cc
from ..signal import window as swindow
from ..keyword_arg import parse_keywords
from .timestamps import Timestamps, Trial, Multi

class Rate(TimeSeries):
    
    def __init__(self,spike_data,**kwargs):
        
        sampling_rate = kwargs.pop('sampling_rate',1000)
        starttime = float(kwargs.setdefault('starttime',0))
        compute = kwargs.pop('compute',True)
        
        # Parse multi-trial spike time data as input
        if isinstance(spike_data,Timestamps):
            multi = spike_data
            data = multi.to_rate(sampling_rate,**kwargs)
        
        # Attempt to convert input to array
        else:
            multi = None
            data = asarray(spike_data)
            
        TimeSeries.__init__(self,data,sampling_rate=sampling_rate,t0=starttime)
        self.multi = multi
        
    def __array__(self):
        return self.data
    
    def truncated(self,n_points):
        '''
        Return a copy whose data is truncated to have only n_points timepoints
        '''
        points = arange(len(self.time))
        data = self.data.compress(points<n_points,axis=-1)
        return type(self)(data,sampling_rate=self.sampling_rate,t0=self.t0)
        
    def padded(self,n_points):
        '''
        Return a copy whose data is zero-padded to have n_points timepoints
        '''
        # Create padding array of zeros
        n_pad = n_points - len(self.time)
        shape = list(self.data.shape)
        shape[-1] = n_pad
        pad = zeros(shape,dtype=self.data.dtype)
        
        # Append zeros and return array
        data = concatenate((self.data,pad),axis=-1)
        return type(self)(data,sampling_rate=self.sampling_rate,t0=self.t0)
         
    def _last_spike(self):
        if self.multi is None:
            return self.data.shape[-1] / float(self.sampling_rate)
        else:
            return self.multi.max_time
    last_spike = property(_last_spike)
        
    def _first_spike(self):
        if self.multi is None:
            return 0
        else:
            return self.multi.min_time
    first_spike = property(_first_spike)
    
    def _get_sampling_rate(self):
        return self.sampling_rate
        
    def _set_sampling_rate(self,rate):
        self.sampling_rate = rate
    t_rate = property(_get_sampling_rate,_set_sampling_rate)
    
    def _starttime(self):
        return self.t0
    starttime = property(_starttime)
    
    def _starttime_float(self):
        return to_float(self.t0)
    starttime_float = property(_starttime_float)
    
    def _endtime(self):
        return self.t0 + self.duration
    endtime = property(_endtime)
    
    def _endtime_float(self):
        return to_float(self.endtime)
    endtime_float = property(_endtime_float)

    def _n_time_points(self):
        return self.data.shape[-1]
    n_time_points = property(_n_time_points)
    
    def _duration_float(self):
        return to_float(self.duration)
    float_duration = property(_duration_float)

class Raster(Rate):
    
    def __init__(self,spike_data,**kwargs):
        '''
        __init__(self,multi,compute=True,trial_len=None,sample_rate=1000,
                 starttime=0,endtime=trial_len)
        '''
        Rate.__init__(self,spike_data,**kwargs)
        
        if kwargs.pop('compute',True):
            if self.data is None:
                self._compute()
        else:
            self.n_bins = None
            
    def _compute(self):
        if self.endtime is None:
            self.endtime = self.last_spike
        else:
            self.endtime = float(self.endtime)
            
        if self.starttime is None:
            self.starttime = min(self.first_spike,0)
                
        self.len = self.endtime - self.starttime
        
        self.n_bins = ceil(self.len*self.t_rate)
        
        if self.multi is None:
            
            self.data = self.data.reshape()
        else:
            self.data = asarray([histogram(asarray(trial.times,dtype='float'),
                                           self.n_bins,[self.starttime,
                                                        self.endtime])[0].squeeze()\
                                 for trial in self.multi.trials],
                                 dtype='float')
    
    def psth(self,win_type=None,**kwargs):
        if win_type is None:
            window = None
        else:
            win_length = round(kwargs.pop('smooth_time',.031) * self.sampling_rate)
            window = swindow(win_type,win_length,normalize=True)
        data = mean(self.data,axis=0)
        return PSTH(data,window,sampling_rate=self.sampling_rate,
                    starttime=self.starttime,endtime=self.endtime)
    
    def two_halves(self,**kwargs):
        trial_sets = bool_(arange(self.n_trials) % 2) # Separate even and odd trials
        half_1 = mean(self.data[-trial_sets,:],axis=0)
        half_2 = mean(self.data[trial_sets,:],axis=0)
        return Raster(vstack([half_1,half_2]),
                      sampling_rate=self.sampling_rate,
                      starttime=self.starttime,
                      endtime=self.endtime)
        
    def plot(self,axes=None,xlabel=True,ylabel=True,frameon=True,xticks=None,
             yticks=None,**kwargs):
        if axes==None:
            axes = gca()
        blend = kwargs.pop('blend',1)
        maxspikes = kwargs.pop('maxspikes',None)
        
        extents = [self.starttime_float,self.endtime_float,
                   self.n_trials+.5,.5]
        im = axes.imshow(repeat(self.data,blend,axis=0),
                         aspect='auto', cmap=cm.bone_r,
                         interpolation='bilinear',extent=extents,
                         origin='upper')
        axes.yaxis.set_major_locator(MultipleLocator(1))
        
        if maxspikes is not None:
            im.set_clim([0,maxspikes])
            
        if xlabel:
            axes.set_xlabel('Time [s]')
        if ylabel:
            axes.set_ylabel('Trial')
        
        axes.set_frame_on(frameon)
        
        if xticks is not None:
            axes.set_xticks(xticks)

        if yticks is not None:
            axes.set_yticks(yticks)
            
        draw_if_interactive()
        return axes
    
    def _n_trials(self):
        return self.data.shape[-2]
    n_trials = property(_n_trials)
    
class PSTH(Rate):
    
    def __init__(self,spike_data,window=None,**kwargs):
        '''
        __init__(self,multi,window=None,compute=True,trial_len=None,
                 sampling_rate=1000,
                 starttime=0,endtime=(trial_len or max(multi)))
        '''
        Rate.__init__(self,spike_data,**kwargs)
        if self.data.ndim > 1:
            self.data = mean(self.data,axis=0)
        
        self.window = window
        
        if kwargs.pop('compute',True):
            self._smooth()
    
    def _smooth(self):
        if self.window is not None:
            self.data = convolve(self.data,self.window.normalize(),'same')
        
    def plot(self,axes=None,max_rate=None,**kwargs):
        if axes==None:
            axes = gca()
        times = to_float(self.time)
        axes.hold(False)
        
        # Get axis keywords
        axes_kwargs, kwargs = parse_keywords(kwargs,['axis'])
        
        '''
        Allow xlabel ylabel, xticks, yticks to be the same as axes_xlabel etc.
        axes_ etc. take precedence.
        '''
        xlabel = kwargs.pop('xlabel',True)
        xlabel = 'Time [s]' if xlabel is True \
                            else '' if xlabel is False \
                            else xlabel
        axes_kwargs.setdefault('xlabel',xlabel)
        
        ylabel = kwargs.pop('ylabel',True)
        ylabel = 'Rate [spikes/s]' if ylabel is True \
                                   else '' if ylabel is False \
                                   else ylabel
        axes_kwargs.setdefault('ylabel',ylabel)
        
        xticks = kwargs.pop('xticks',True)
        if xticks is not True:
            axes_kwargs['xticks'] = xticks
            
        yticks = kwargs.pop('yticks',True)
        if yticks is not True:
            axes_kwargs['yticks'] = yticks
        
        # Plot the PSTH line
        axes.plot(times,self.t_rate*self.data,**kwargs)
        
        # Set tight limits on axis
        axes_kwargs['xlim'] = [self.starttime_float,self.endtime_float]
        
        # Set limit on y axis
        if max_rate is not None:
            axes_kwargs['ylim'] = [0,max_rate]
        
        # Set all axis properties
        setp(axes,**axes_kwargs)
            
        draw_if_interactive()
        return axes

def raster(spike_times,sampling_rate=1000,**kwargs):
    kwargs.update({'sampling_rate':sampling_rate})
    return Raster(Multi(spike_times),**kwargs)

def psth(spike_times,win_type='hanning',**kwargs):
    '''
    psth(spike_times,window,sampling_rate=1000,trial_len=None,pre=0,
         post=(trial_len or max(spike_times)),
    '''
    
    # Handle both old-style 'sample_rate' and new 'sampling_rate'
    sample_rate = kwargs.pop('sample_rate',1000)
    sampling_rate = kwargs.setdefault('sampling_rate',sample_rate)
    
    # Make window if specified
    if win_type is None:
        window = None
    else:
        win_length = round(kwargs.pop('smooth_time',.031) * sampling_rate)
        window = swindow(win_type,win_length,normalize=True)
    
    return PSTH(Multi(spike_times),window,**kwargs)

def average_psth(*args):
    
    # Start by copying first psth object
    p = copy(args[0])
    
    # We need to make a copy of the 'multi', 'data', 'trials' subobjects
    p.multi = copy(p.multi)
    if p.multi is None:
        p.data = copy(p.data)
        
        # Take the mean of the raw rate data
        p.data = mean(vstack(a.data for a in args),0)
        
    else:
        p.multi.trials = copy(p.multi.trials)
        
        # Append all of the spiketimes together
        for psth_obj in args[1:]:
            p.multi.trials.extend(psth_obj.multi.trials)
            
    p.data = copy(p.data)
    
    # Decide how to compute the final data
    do_average = not any(a.data is None for a in args)
    
    if do_average:
        all_data = asarray([a.data for a in args])
        p.data = mean(all_data,0)
    else:
        p._compute()
        
    return p

def concatenate_rasters(rasters,truncate=False):
    
    # Get sampling rates and numbers of trials
    (rates,n_trials) = zip(*((r.sampling_rate,r.n_trials) for r in rasters))
    
    # Check that all rasters have the same sampling rate
    if len(set(rates)) > 1:
        raise ValueError, "Rasters must all have the same sampling rate"
    else:
        rate = rates[0]
    
    # Check that all rasters have the same number of trials
    if len(set(n_trials)) > 1:
        if truncate:
            n_trials = min(n_trials)
        else:
            raise ValueError, "If truncate=False, all rasters must have the same number of trials"
    else:
        n_trials = n_trials[0]

    # Create new raster object by concatenating the data in order
    data = hstack([r.data[:n_trials,:] for r in rasters])
    return Raster(data,sampling_rate=rate)
        
def cc(*args):
    return sequence_cc(*[a.data for a in args])

def correlate(b,a):
    '''
    Alternate correlation metric suggested by reviewers for masker paper:
    instead of dot(a,b)/sqrt(dot(a,a) * dot(b,b)), return dot(a,b) / dot(a,a)
    '''
    ma = mean(a.data)
    mb = mean(b.data)
    denom = dot(a.data-ma,a.data-ma)
    if denom != 0:
        return dot(a.data-ma,b.data-mb) / denom
    else:
        return NaN

def to_float(time_obj):
    return asarray(time_obj,'float')/time_obj._conversion_factor
