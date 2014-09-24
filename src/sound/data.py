from __future__ import absolute_import
import os
from copy import copy, deepcopy
from math import floor, ceil, pi
from numpy import (array, asarray, arange, diff, squeeze, zeros,
                   matrix, mean, diag, max, log10, hstack, round)
from numpy.fft import fft, fftshift, fftfreq, fft2
from scipy.io.wavfile import write as wav_write
from matplotlib.colors import ListedColormap
from matplotlib.pyplot import imshow, draw_if_interactive, gca
from . import colormap
from ..stats import cc
from ..signal import gaussian_window

__all__ = ['PCM','Spectrogram']

class TimeDomain(object):
    '''
    Time-domain waveforms
    '''
    
    def __init__(self,data,samplerate):
        self.data = squeeze(asarray(data))
        self.samplerate = samplerate
        
    def __len__(self):
        return len(self.data)
         
    def _timeslice(self,window,framerate,zeropad=True):
        
        # Hanlde window objects, arrays, or integer lengths
        if isinstance(window,int):
            n_points = window
            do_window = False
        elif isinstance(window,float):
            raise ValueError,'input argument "window" must be either\
                              an integer or an array type'
        else:
            n_points = len(window)
            do_window = True
            
        # Zeropad if needed        
        if zeropad:
            time_data = self._zeropad(floor(n_points/2.),
                                      ceil(n_points/2.))
        else:
            time_data = self.data
        
        # Compute parameters and initialize output
        increment = float(self.samplerate) / framerate
        n_frames = int(floor((len(time_data) - n_points) / increment) + 1)
        frame_starts = increment * arange(n_frames)
        output = zeros([n_points,n_frames])
        
        # Assign output
        for s in range(n_frames):
            start = increment * s
            output[:,s] = time_data[start:start+n_points]
            
        # Window if needed:
        if do_window:
            return asarray(matrix(diag(window)) * matrix(output))
        else:
            return output
        
    def _zeropad(self,npre,*args):
        if len(args) > 0: npost = args[0]
        else: npost = npre
        padded = zeros(len(self) + npre + npost)
        padded[npre:npre+len(self)] = self.data
        return padded

class PCM(TimeDomain):
    '''
    PCM representations, including standard audio files
    '''

    def __init__(self,data,samplerate,**kwargs):
        TimeDomain.__init__(self,data,samplerate)
        self.file = kwargs.pop('file',None)
    
    # Testing removal 1/10/11 RCM
    #def spectrogram(self,window,slice_rate,logscale=True):
    #    sliced = self._timeslice(window,slice_rate)
    #    if logscale:
    #        return Spectrogram(log10(abs(fft(sliced,len(window),0))),
    #                           slice_rate)
    #    else:
    #        return Spectrogram(abs(fft(sliced,len(window),0)),slice_rate)
        
    def write_wav(self,filename):
        wav_write(filename,round(float(self.samplerate)),self.data)

class TimeFrequency(object):
    '''
    Time-frequency representations
    '''
    freq_axis = 0
    time_axis = 1
    
    def __init__(self,sound,slice_rate,**kwargs):
        self.sound = sound
        self.t_rate = slice_rate
        
    def add_silence(self,**kwargs):
        pre = kwargs.pop('pre',0)
        post = kwargs.pop('post',pre)
        nbands = self.data.shape[0]
        
        if pre > 0:
            nsamples = pre * self.t_rate
            self.data = hstack((zeros([nbands,nsamples]),self.data))
            newtimes = arange(-nsamples,0,dtype='float')/self.t_rate
            self.times = hstack((newtimes,self.times))
            
        if post > 0:
            nsamples = post * self.t_rate
            self.data = hstack((self.data,zeros([nbands,nsamples])))
            newtimes = arange(1,nsamples+1,dtype='float')/self.t_rate\
                       + max(self.times)
            self.times = hstack((self.times,newtimes))
            
        return self

    def subtract_constant(self,band_constants):
        
        '''Subtract a constant offset from each frequency band'''
        
        if self.n_bands != len(band_constants):
            raise ValueError, "Input 'band_constants' is the wrong length"
        for jj,const in enumerate(band_constants):
            self.set_band(jj,self.get_band(jj) - const)
            
    def constant_subtracted(self,band_constants):
        
        '''
        Return a new spectrogram with the constants in band_constants
        subtracted from the frequency bands
        '''
        
        sub_spec = self._copy()
        sub_spec.subtract_constant(band_constants)
        return sub_spec
    
    def mean_subtracted(self):
        return self.constant_subtracted(self.band_means)
    
    def get_band(self,index):
        return self.data[index,:]

    def set_band(self,index,data):
        self.data[index,:] = data
    
    @property   
    def band_means(self):
        return mean(self.data,self.freq_axis)
    
    @property
    def n_bands(self):
        return self.data.shape[self.freq_axis]
    
    @property
    def n_timepoints(self):
        return self.data.shape[self.time_axis]
    
    def _timeslice(self,window,framerate,zeropad=True):
        '''
        Return array of spectrogram slices, possibly windowed
        '''
        
        # Handle window objects, arrays, or integer lengths
        if isinstance(window,int):
            n_points = window
            do_window = False
        elif isinstance(window,float):
            raise ValueError,'input argument "window" must be either\
                              an integer or an array type'
        else:
            n_points = len(window)
            do_window = True
            
        # Zeropad if needed
        if zeropad:
            tf_data = self._zeropad(floor(n_points/2.),
                                      ceil(n_points/2.))
        else:
            tf_data = self.data
        
        # Compute parameters
        increment = float(self.t_rate) / framerate
        n_frames = int(floor((self.n_timepoints - n_points) / increment) + 1)
        frame_starts = increment * arange(n_frames)
        
        # Compute output shape and initialize
        out_shape = [n_frames,0,0]
        out_shape[self.time_axis + 1] = n_points
        out_shape[self.freq_axis + 1] = self.n_bands
        output = zeros(out_shape)
        
        # Assign output
        slice_min = [0,0]
        slice_max = [0,0]
        slice_max[self.freq_axis] = self.n_bands
        out_min = [0,0]
        out_max = [0,0]
        out_max[self.freq_axis] = self.n_bands
        for s,start in enumerate(frame_starts):
            slice_min[self.time_axis] = start
            slice_max[self.time_axis] = start + n_points
            out_min[self.time_axis] = 0
            out_max[self.time_axis] = n_points
            output[s,out_min[0]:out_max[0],
                     out_min[1]:out_max[1]]\
                = tf_data[slice_min[0]:slice_max[0],
                          slice_min[1]:slice_max[1]]
            
        # Window if needed:
        if do_window:
            return output * window
        else:
            return output
        
    def _zeropad(self,npre,*args):
        if len(args) > 0: npost = args[0]
        else: npost = npre
        
        # Set output shape
        padded_shape = [0,0]
        padded_shape[self.freq_axis] = self.n_bands
        padded_shape[self.time_axis] = self.n_timepoints + npre + npost
        
        # Create empty matrix for output
        padded = zeros(padded_shape)
        
        # Get output range
        mins = [0,0]
        mins[self.time_axis] = npre
        maxs = [-1,-1]
        maxs[self.time_axis] = npre + self.n_timepoints
        maxs[self.freq_axis] = self.n_bands
        
        # Set output shape
        padded[mins[0]:maxs[0],mins[1]:maxs[1]] = self.data
        return padded

class Spectrogram(TimeFrequency):
    ''' 
    Spectrograms
    '''
    
    def __init__(self,sound,window,slice_rate,**kwargs):
        TimeFrequency.__init__(self,sound,slice_rate,**kwargs)
        self.window = window
        
        self.logscale = kwargs.pop('logscale',True)
        self.freqscale = kwargs.pop('freqscale','linear')        
        compute = kwargs.pop('compute',True)
        self.f_lims = kwargs.pop('freq_lims',[0,float(sound.samplerate)/2])
        self.db_noise = kwargs.pop('db_noise',40)
        
        if self.freqscale != 'linear':
            raise ValueError,'nonlinear frequency scale warpings\
                              are not supported yet'
        
        if compute:
            self._compute()
        else:
            self.data = None
            self.freqs = None
        
    def _compute(self):
        
        # Slice and compute FFT
        sliced = self.sound._timeslice(self.window,self.t_rate)
        self.times = arange(sliced.shape[1],dtype='float')/self.t_rate
        self.data = fftshift(fft(sliced,axis=0),axes=[0])
        self.data2 = self.data.copy()
        freqs = fftshift(fftfreq(len(self.window),
                                 1/float(self.sound.samplerate)))
        
        # Select correct frequencies
        keep_freqs = (freqs >= self.f_lims[0]) & (freqs <= self.f_lims[1])
        self.data = abs(self.data[keep_freqs,:])
        self.freqs = freqs[keep_freqs]
        
        if self.logscale:
            self.data = 20*log10(self.data)
            self.data = self.data - max(self.data) + self.db_noise 
            self.data[self.data<0] = 0
            
    def _copy(self):
        '''
        Sane 1-deep copy operation
        Avoids issues of copying data buffers; 'sound' is the only attribute
        left as a link to the original.
        '''
        if self.__class__ is Spectrogram:
            s = Spectrogram(self.sound,copy(self.window),copy(self.t_rate),
                            logscale=copy(self.logscale),
                            freqscale=copy(self.freqscale),
                            freq_lims=copy(self.f_lims),
                            db_noise=copy(self.db_noise),
                            compute=False)
            s.freqs = copy(self.freqs)
            s.data = copy(self.data)
            return s
        else:
            return copy(self)
            
    def plot(self,axes=None,cmap=colormap.fet,xlabel=True,ylabel=True,
             title=True,xticks=None,yticks=None,frameon=True,freq_tick_scale=1,
             freq_decimal_points=0,**kwargs):
        if axes == None:
            axes = gca()
        extents = [min(self.times),max(self.times),min(self.freqs),max(self.freqs)]
        axes.imshow(self.data,aspect='auto',origin='lower',
                    extent=extents,cmap=cmap,**kwargs)
        
        # Default x and y axis labels
        if xlabel is True:
            xlabel = 'Time [s]'
        elif xlabel is False:
            xlabel = ''
        axes.set_xlabel(xlabel)

        if ylabel is True:
            ylabel = 'Frequency [Hz]'
        elif ylabel is False:
            ylabel = ''
        axes.set_ylabel(ylabel)
        
        if title and (self.sound.file is not None):
            axes.set_title(os.path.split(self.sound.file.name)[1])
        
        axes.set_frame_on(frameon)
            
        if xticks is not None:
            axes.set_xticks(xticks)
            
        if yticks is not None:
            axes.set_yticks(yticks)
            freq_label_format = '%%.%df' % freq_decimal_points
            labels = ['%.0f' % (t/freq_tick_scale) for t in yticks]
            axes.set_yticklabels(labels)
        
        draw_if_interactive()
        return axes
    
    @property
    def band_spacing(self):
        if self.freqscale.lower() == 'linear':
            return mean(diff(self.freqs))
        else:
            return None
    
    def _sample_modspec(self,*args,**kwargs):
        '''
        Return samples of complex modulation spectrum
        '''
        sliced = self._timeslice(*args,**kwargs)
        winlength = sliced.shape[self.time_axis+1]
        tf = fftfreq(winlength,1./self.t_rate)
        sf = fftfreq(self.n_bands,self.band_spacing)
        ft = fft2(sliced)
        return ft, tf, sf
    
    def sample_mps(self,window,temp_rate_Hz=10,zeropad=True):
        '''
        Return samples of modulation power spectrum
        '''
        
        ft, tf, sf = self._sample_modspec(window,temp_rate_Hz,
                                          zeropad=zeropad)
        return fftshift(abs(ft)**2), fftshift(tf), fftshift(sf)
    
    def calc_mps(self,*args,**kwargs):
        '''
        Return estimate of modulation power spectrum
        '''
        
        mps, tf, sf = self.sample_mps(*args,**kwargs)
        return mean(mps,0), tf, sf
        
def gaussian_spectrum(sound,slice_rate=1000,win_length=256,**kwargs):
    band_width = kwargs.pop('band_width',None)
    if band_width is None:
        return Spectrogram(sound,gaussian_window(win_length),slice_rate,**kwargs)
    else:
        pass

def strfpak_resample(sound,slice_rate,resample_type):
    
    '''
    Resample a sound object using scikits.samplerate.resample.
    Sound is upsampled to nearest multiple of slice_rate;
    if sound.samplerate == slice_rate, no resampling is done.
    Sound is also converted to float.
    Currently, error is less than 10dB max, mostly very good.
    '''
    
    from scikits.samplerate import resample
    from copy import deepcopy
    
    # Calculate the resampling frequency for STRFPAK 5.3:
    # round up to the nearest multiple of 'slice_rate'
    input_freq = floor(sound.samplerate)
    output_freq = ceil(input_freq/slice_rate)*slice_rate
    
    # Copy so that we don't alter original object
    sound2 = deepcopy(sound)
    
    # Scale factor -- MATLAB makes .wav data be -1:1, scipy does integer
    sound2.data = sound2.data/float(2**15)
    
    # Resample if necessary
    if output_freq > input_freq:
        sound2.data = resample(sound2.data,output_freq/input_freq,
                               resample_type)
        sound2.samplerate = output_freq
    
        # Pad because matlab resampling works differently
        sound2.data = sound2._zeropad(0,2)
    
    return sound2
    
def strfpak_spectrum(sound,slice_rate=1000,bandwidth_hz=125,
                     resample_type='sinc_best',**kwargs):
    
    '''
    Compute spectrogram using the same parameters as STRFPAK 5.3
    As of 1/11/11, mean-squared-errors are as follows:
        2e-5 between computed spectrogram and matlab's version
        3e-7 between spectrogramp computed from matlab's resampling
            and matlab's version
        Thus some of the error is actually due to resampling.
        In all cases, maimum errors are 8dB on a 9dB signal, very near the
            noise floor.
    '''
    
    # STRFPAK default is 80dB range; default for Spectrogram is 40
    kwargs.setdefault('db_noise',80)
    
    # Resample
    sound2 = strfpak_resample(sound,slice_rate,resample_type)
    
    # Compute window length for requested bandwidth
    kwargs['win_length'] = floor(1/(2*pi*bandwidth_hz)*6*float(sound2.samplerate))
    
    # Zero-pad in accordance with STRFPAK
    #if win_length % 2:
    #    n_pre = (win_length+1)/2
    #    n_post = (win_length+1)/2-1
    #else:
    #    n_pre = win_length/2
    #    n_post = win_length/2
    #sound2.data = sound2._zeropad(n_pre,n_post)
    
    # Compute number of frames for STRFPAK
    increment = sound2.samplerate/slice_rate
    frame_count = floor(float(len(sound2))/increment) + 1
    
    # Compute spectrogram
    spec = gaussian_spectrum(sound2,slice_rate,**kwargs)
    
    # Trim to STRFPAK size
    spec.data = spec.data[:,:frame_count]
    return spec

def spec_si(*specs):
    return cc(*(s.data.flatten() for s in specs))

def time_si(*sounds):
    if sounds[1].samplerate == sounds[0].samplerate:
        if len(sounds[1]) == len(sounds[0]):
            return cc(*(s.data for s in sounds))
        else:
            return correlate(*(s.data for s in sounds))
    else:
        raise ValueError, "Sounds do not have the same sampling rate"
    
def si(*sounds,**kwargs):
    si_type = kwargs.pop('si_type','spec')
    if si_type == 'spec':
        return spec_si(*(gaussian_spectrum(s,**kwargs) for s in sounds))
    elif si_type == 'time':
        return time_si(*sounds)
    else:
        raise ValueError, "Argument 'si_type' must be one of 'time','spec'"
