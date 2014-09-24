from numpy import arange, asarray, ceil, diff, histogram, hstack, vstack
from matplotlib.pyplot import hist
from nitime.timeseries import TimeSeries

class Timestamps(object):
    
    def __init__(self,sample_rate_Hz=1,starttime=0,endtime=None):
        self.sample_rate_Hz=sample_rate_Hz
        self.starttime = starttime
        self.endtime = endtime
    
    def _end(self):
        return self.endtime or self.max_time
    end = property(_end)
    
    def _duration(self):
        return self.end - self.starttime
    duration = property(_duration)
    
    def interspike_interval_histogram(self,resolution=.0001,max_time=.1):
        # resolution, max_time are in seconds
        bins = int(ceil((max_time*self.sample_rate_Hz/float(resolution))))
        return histogram(self.interspike_times(),bins,(0,max_time))
    
    def to_timeseries(self,sampling_rate,starttime=None,endtime=None):
        rate_kwargs = {'sampling_rate':sampling_rate,
                       'starttime':starttime,
                       'endtime':endtime}
        ts_kwargs = {'sampling_rate':sampling_rate,
                     't0':starttime}
        return TimeSeries(self.to_rate(**rate_kwargs),**ts_kwargs)
    
    def plot_isi_histogram(self,resolution=.0001,max_time=.1):
        n_bins = int(ceil((max_time*self.sample_rate_Hz/float(resolution))))
        bins = arange(0,n_bins)*resolution
        return hist(self.interspike_times(),bins)
    
class Trial(Timestamps):
    
    def __init__(self,times,**kwargs):
        Timestamps.__init__(self,**kwargs)
        self.times = asarray(times).flatten()
        
    def _nspikes(self):
        return len(self.times)
    nspikes = property(_nspikes)
    
    def _min_time(self):
        return min(hstack((self.times,0)))
    min_time = property(_min_time)
    
    def _max_time(self):
        return max(hstack((self.times,0)))
    max_time = property(_max_time)

    def interspike_times(self):
        return self.sample_rate_Hz * diff(self.times)
        
    def to_rate(self,sampling_rate,starttime=None,endtime=None):
        
        # Allow for external setting of start and end times
        if starttime is None:
            starttime = self.starttime
        
        if endtime is None:
            endtime = self.end
        
        duration = endtime - starttime
        
        # Compute number of bins for output rate
        n_bins = ceil(float(duration)*float(sampling_rate))
        
        
        return histogram(asarray(self.times,dtype='float'),
                         n_bins,
                         [float(starttime),float(endtime)])[0].squeeze()

class Multi(Timestamps):
    
    def __init__(self,trials,sample_rate_Hz=1,**kwargs):
        Timestamps.__init__(self,sample_rate_Hz=sample_rate_Hz)
        self.trials = [Trial(trial,sample_rate_Hz=sample_rate_Hz) \
                       for trial in trials]
            
    def _nspikes(self):
        return sum([t.nspikes for t in self.trials])
    nspikes = property(_nspikes)
    
    def _ntrials(self):
        return len(self.trials)
    ntrials = property(_ntrials)
            
    def _max_time(self):
        if self.ntrials > 0:
            return max([trial.max_time \
                        for trial in self.trials \
                        if trial.nspikes>0]\
                       + [None])
        else:
            return None
    max_time = property(_max_time)
    
    def _min_time(self):
        if self.ntrials > 0:    
            return min([trial.min_time \
                        for trial in self.trials \
                        if trial.nspikes>0])
        else:
            return None
    min_time = property(_min_time)
        
    def interspike_times(self):
        return concatenate(t.interspike_times() for t in self.trials)
    
    def to_rate(self,sampling_rate,**kwargs):
        return vstack([trial.to_rate(sampling_rate,**kwargs) for trial in self.trials])