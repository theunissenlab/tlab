from numpy import array, concatenate, interp, nonzero
from matplotlib.pyplot import gca, setp
from . import STRF
from ..stats.fit import fitgaussian

def eMTF(mtf_list):
    return n.average(n.concatenate(mtf_list,2),2)

class MTF(object):
    
    def __init__(self,data,temp_mods,spec_mods):
        self.data = data
        self.spec_mods = spec_mods
        self.temp_mods = temp_mods
        
    def plot(self,axes=None,**kwargs):
        
        # Determine scale factor for spectral units
        spec_units = kwargs.pop('spec_units','cyc_per_kHz')
        
        if spec_units.lower() == 'cyc_per_khz':
            sf = 1000
        else:
            sf = 1
            
        if axes is None:
            axes = gca()
        
        tlarg = {'temp_lim':kwargs.pop('temp_lim',100)}
        
        im= axes.imshow(self._center_data(**tlarg),
                        origin='lower',
                        aspect='auto',
                        extent=[self._center_temps(**tlarg)[0],
                                self._center_temps(**tlarg)[-1],
                                self.spec_mods_pos[0]*sf,
                                self.spec_mods_pos[-1]*sf])
        setp(axes,**kwargs)
        return im
    
    def _spec_mod_is_positive(self):
        return (self.spec_mods>=0).flatten()
    
    def _temp_mod_is_mid(self,temp_lim=100):
        return (abs(self.temp_mods)<=temp_lim).flatten()
    
    def _upper_freqs(self):
        return self.spec_mods[self._spec_mod_is_positive()]
    
    def _center_temps(self,**kwargs):
        return self.temp_mods[self._temp_mod_is_mid(**kwargs)]
    
    def _upper_data(self):
        return self.data[self._spec_mod_is_positive(),:]
    
    def _center_data(self,**kwargs):
        return self.upper_data[:,self._temp_mod_is_mid(**kwargs)]
    
    def get_data(self,temp_lims,spec_lims):
        
        temp = self._temp_mask(temp_lims)
        spec = self._spec_mask(spec_lims)
        
        return self.data[spec,:][:,temp]
    
    def _spec_mask(self,spec_lims):
        spec_mods = array(self.spec_mods)
        return array([spec_mods >= spec_lims[0],
                      spec_mods <= spec_lims[1]]).all(0)
                      
    def _temp_mask(self,temp_lims):
        temp_mods = array(self.temp_mods)
        return array([temp_mods >= temp_lims[0],
                      temp_mods <= temp_lims[1]]).all(0)
                      
    def get_spec_mods(self,spec_lims):
        return self.spec_mods[self._spec_mask(spec_lims)]
    
    def get_temp_mods(self,temp_lims):
        return self.temp_mods[self._temp_mask(temp_lims)]
    
    def fit_gaussian(self):
        [p,x,y,width_x,width_y] = fitgaussian(self.upper_data)
        wt = interp(y,range(len(self.temp_mods)),self.temp_mods)
        wf = interp(x,range(len(self.spec_mods_pos)),self.spec_mods_pos)
        return wt,wf
    
    def best_freqs(self):
        m = self.upper_data.max()
        x,y = nonzero(self.upper_data == m)
        wt = self.temp_mods[y[0]]
        wf = self.spec_mods_pos[x[0]]
        return wt,wf
    
    upper_data = property(_upper_data)
    center_data = property(_center_data)
    spec_mods_pos = property(_upper_freqs)
    temp_mods_mid = property(_center_temps)
    
def mtf(strf):
    return MTF(*strf.calc_mptf())