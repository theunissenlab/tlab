from numpy import (average, diff, abs, max, ones, round, vstack, zeros, sum,
                   convolve, array)
from numpy.fft import fft2, fftfreq, fftshift
import matplotlib as mpl
from matplotlib import cm
from matplotlib.pyplot import gca, gcf, setp, Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from .SI import similarity_index as si
from ..sound import colormap
from ..stats.fit import fitgaussian

class STRF(object):
    
    def __init__(self,data,freqs,lags):
        self.data = data
        self.freqs = freqs
        self.lags = lags
        
    def fit_gaussian_2d(self):
        pass
    
    def _lag_is_causal(self):
        return (self.lags>=0).flatten()
    
    def _n_causal(self):
        return sum(self._lag_is_causal())
    
    def _n_bands(self):
        return len(self.freqs)
    
    def _omega_t(self):
        return average(diff(self.lags))
    
    def _omega_f(self):
        return average(diff(self.freqs))
    
    def _absmax(self):
        return max(abs(self.data))
    
    n_causal = property(_n_causal)
    n_bands = property(_n_bands)
    omega_t = property(_omega_t)
    omega_f = property(_omega_f)
    absmax = property(_absmax)
    
    def predict(self,spectrogram,causal=True):
        nbands = self.n_bands
        spec_nbands = len(spectrogram.freqs)
        if causal:
            filter = self.causal()
        else:
            filter = self.data

        if spec_nbands == nbands:
            band_preds = [convolve(spectrogram.data[jj],filter[jj],'same') for jj in range(nbands)]
            return sum(vstack(band_preds),0)
        else:
            raise ValueError, "Spectrogram and STRF have different numbers of bands (%d,%d)" % (spec_nbands,
                                                                                                nbands)
    
    def causal(self):
        return self.data[:,self._lag_is_causal()]
    
    def causal_lags(self):
        return self.lags[self._lag_is_causal()]
    
    def pad_causal(self):
        '''
        Pad data so that the zero band is in the data
        '''
        f_min = min(self.freqs)
        n_pad = round(f_min/self.omega_f)
        pad = zeros((n_pad,self.n_causal))
        return vstack([pad,self.causal()])
        
    def calc_mtf(self):
        '''
        mtf, omega_t, omega_f = calc_mtf(self)
        '''

        return (fft2(self.causal()),
                fftfreq(self.n_causal,self.omega_t),
                fftfreq(self.n_bands,self.omega_f))
        
    def calc_mptf(self):
        mtf,tf,sf = self.calc_mtf()
        return abs(fftshift(mtf*mtf.conj())),fftshift(tf),fftshift(sf)
        
    def plot(self,pos_only=False,max_lag=None,min_lag=None,axes=None,
             time_units='ms',freq_units='kHz',show_xlabel=True,
             show_ylabel=True,**kwargs):
        if axes is None:
            axes=gca()
        
        defaults = {'aspect':'auto',
                    'origin':'lower',
                    'cmap':cm.jet,
                    'vmin':-self.absmax,
                    'vmax':self.absmax}
        for k,v in defaults.items():
            kwargs.setdefault(k,v)
        
        if max_lag is None:
            low_lags = ones(self.lags.shape,dtype='bool').squeeze()
            max_lag = max(self.lags)
        else:
            low_lags = (abs(self.lags) <= max_lag).squeeze()
        
        if pos_only:
            min_lag = 0
        
        if min_lag is None:
            good_lags = low_lags
            min_lag = min(self.lags)
        else:
            high_lags = (self.lags >= min_lag).squeeze()
            good_lags = low_lags * high_lags

        data = self.data[:,good_lags]
        
        # Compute time axis scaling factor
        if time_units == 'ms':
            time_scale = 1000
        elif time_units == 's':
            time_scale = 1
        else:
            raise ValueError, 'Time unit "%s" not recognized' % time_units
        
        # Compute frequency axis scaling factor
        if freq_units.lower() == 'khz':
            freq_scale = .001
        elif freq_units.lower() == 'hz':
            freq_scale = 1
        else:
            raise ValueError, 'Frequency unit "%s" not recognized' % freq_units
        
        # Set extents for display
        extent = [min_lag * time_scale,
                  max_lag * time_scale,
                  self.freqs[0] * freq_scale,
                  self.freqs[-1] * freq_scale]
        
        # Add x and y axis labels
        if show_xlabel:
            axes.set_xlabel('Time [%s]' % time_units)
        if show_ylabel:
            axes.set_ylabel('Frequency [%s]' % freq_units)
        
        return axes.imshow(data,extent=extent,**kwargs)
    
    def plot_with_marginals(self,fig=None,**kwargs):
        
        strf_kwargs = kwargs.copy();
        
        # Compute time axis scaling factor
        time_units = kwargs.pop('time_units','ms')
        if time_units == 'ms':
            time_scale = 1000
        elif time_units == 's':
            time_scale = 1
        else:
            raise ValueError, 'Time unit "%s" not recognized' % time_units
        
        # Compute frequency axis scaling factor
        freq_units = kwargs.pop('freq_units','kHz')
        if freq_units.lower() == 'khz':
            freq_scale = .001
        elif freq_units.lower() == 'hz':
            freq_scale = 1
        else:
            raise ValueError, 'Frequency unit "%s" not recognized' % freq_units
        
        if fig is None:
            fig = gcf()
        
        # Make axes
        strf_ax = fig.add_axes([.1,.3,.6,.6])
        spec_ax = fig.add_axes([.75,.3,.15,.6],sharey=strf_ax)
        temp_ax = fig.add_axes([.1,.1,.6,.15],sharex=strf_ax)
        
        # Plot STRF
        self.plot(axes=strf_ax,pos_only=True,**strf_kwargs)
        xt = strf_ax.get_xticklabels()
        yt = strf_ax.get_yticklabels()
        setp(xt,'visible',False)
        setp(yt,'visible',False)
        strf_ax.xaxis.set_ticks_position('none')
        strf_ax.yaxis.set_ticks_position('none')
        
        # Plot marginal spectrum
        spec = self.causal().mean(1)
        maxspec = max(abs(spec))
        freqs = array(self.freqs) * freq_scale
        spec_ax.plot(spec,freqs)
        spec_ax.set_xlim([-1.1*maxspec,1.1*maxspec])
        spec_ax.set_xticks([0])
        spec_ax.set_ylim([min(freqs),max(freqs)])
        spec_ax.yaxis.set_ticks_position('right')
        
        # Plot time marginal
        temp = self.causal().mean(0)
        maxtemp = max(abs(temp))
        lags = array(self.causal_lags()) * time_scale
        temp_ax.plot(lags,temp)
        temp_ax.set_yticks([0])
        temp_ax.set_xlim([min(lags),max(lags)])
        temp_ax.set_ylim([-1.1*maxtemp,1.1*maxtemp])
        temp_ax.xaxis.set_ticks_position('bottom')

    def saveplot(self,filename,**kwargs):
        dpi = kwargs.pop('dpi',72)
        format = kwargs.pop('format','png')
        backend_canvas = kwargs.pop('backend_canvas',FigureCanvasAgg)
        
        ins_orig = mpl.rcParams['svg.image_noscale']
        mpl.rcParams['svg.image_noscale'] = True
        
        fig = Figure()
        canvas = backend_canvas(fig)
        ax = fig.add_subplot(111)
        self.plot(axes=ax,**kwargs)
        fig.savefig(filename,dpi=dpi,format=format)
        
        # Put svg image config back to original
        #mpl.rcParams['svg.embed_char_paths'] = ecp_orig
        mpl.rcParams['svg.image_noscale'] = ins_orig
        
        return filename