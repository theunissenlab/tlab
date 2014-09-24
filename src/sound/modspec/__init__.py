from numpy import array, concatenate, log10, mean
from matplotlib.pyplot import gca
from ...keyword_arg import parse_keywords
from ..data import gaussian_spectrum, gaussian_window
  
class ModSpec(object):
    '''
    Modulation spectrum for a sound
    '''
    
    def __init__(self,data,temp_freqs,spec_freqs):
        self.data = data
        self.temp_freqs = temp_freqs
        self.spec_freqs = spec_freqs
    
    def plot(self,ax=None,temp_lims=[-100,100],spec_lims=[0,3],
             spec_units='cycles per kHz',**kwargs):
        
        # Default to plot to current axes
        if ax is None:
            ax = gca()
        
        # Handle frequency scale units
        if spec_units.lower() == 'cycles per khz':
            spec_scale = 1000.
        elif spec_units.lower() == 'cycles per hz':
            spec_scale = 1.
        else:
            raise ValueError, \
                "Spectral scale units '%s' not recognized" % spec_units
        spec_lims_cph = [l/spec_scale for l in spec_lims]
        
        # Get data and extents
        data = self.get_data(temp_lims,spec_lims_cph)
        extent=[l for s in [temp_lims,spec_lims] for l in s]
        
        # Handle log scale
        do_log = kwargs.pop('logscale',True)
        if do_log:
            data = log10(data)
        
        # Plot        
        ax.imshow(data,extent=extent,aspect='auto',origin='lower',**kwargs)
        return ax
    
    def get_data(self,temp_lims,spec_lims):
        
        temp = self._temp_mask(temp_lims)
        spec = self._spec_mask(spec_lims)
        
        return self.data[spec,:][:,temp]
    
    def _spec_mask(self,spec_lims):
        spec_freqs = array(self.spec_freqs)
        return array([spec_freqs >= spec_lims[0],
                      spec_freqs <= spec_lims[1]]).all(0)
                      
    def _temp_mask(self,temp_lims):
        temp_freqs = array(self.temp_freqs)
        return array([temp_freqs >= temp_lims[0],
                      temp_freqs <= temp_lims[1]]).all(0)
                      
    def get_spec_freqs(self,spec_lims):
        return self.spec_freqs[self._spec_mask(spec_lims)]
    
    def get_temp_freqs(self,temp_lims):
        return self.temp_freqs[self._temp_mask(temp_lims)]
    
def emps(sounds,window,temp_rate_Hz,**kwargs):
    
    spec_kwargs,kwargs = parse_keywords(kwargs,['spec'])
        
    samples = [gaussian_spectrum(s,**spec_kwargs)\
                   .sample_mps(window,temp_rate_Hz,**kwargs) \
               for s in sounds]
    mod_data = concatenate([s[0] for s in samples],0)
    tf = samples[0][1]
    sf = samples[0][2]
    return ModSpec(mean(mod_data,0),tf,sf)

def gaussian_emps(sounds,win_duration=0.2,**kwargs):
    
    # Construct gaussian window of specified duration
    spec_t_rate = kwargs.setdefault('spec_samplerate_Hz',1000)
    win_length = int(round(win_duration * spec_t_rate))
    window = gaussian_window(win_length)
    
    # Default to 50% overlap if temp_rate_Hz is not specified
    temp_rate_Hz = kwargs.pop('temp_rate_Hz',round(2./win_duration))
    
    return emps(sounds,window,temp_rate_Hz,**kwargs)
    
def gaussian_mps(spectrogram,win_duration=0.2,**kwargs):
    
    
    win_length = int(round(win_duration * spectrogram.t_rate))
    window = gaussian_window(win_length)
    
    # Default to 50% overlap
    temp_rate_Hz = kwargs.pop('temp_rate_Hz',round(2./win_duration))
    
    return ModSpec(*spectrogram.calc_mps(window,temp_rate_Hz,**kwargs))