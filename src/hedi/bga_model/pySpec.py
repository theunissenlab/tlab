from numpy import *
import pyBGA

sampleRate=44100.0
nstd=6
fband=125
f=1000.0
freq_min=200
freq_max=8000
twindow = nstd/(fband*2.0*pi)
winLength = round(twindow*sampleRate);
winLength = round(winLength/2)*2; 
increment = round(sampleRate/f); 



class pySpec(object):
	def __init__(self,wave_file):
		self.wave_file=wave_file[:]
		spec=pyBGA.pyGaussianSpectrum(wave_file,increment,winLength)
		self.whole_spec=spec[:]		
		fft,fc=shape(spec)
		t=arange(0,fc)/f
		fs=arange(0,fft)*sampleRate/fft
		a= sum(abs(spec),0)
		imin=ceil(freq_min/sampleRate*fft);
		imax=ceil(freq_max/sampleRate*fft);
		fs=fs[imin:imax]
		t=t[abs(a)>max(a)/10]
		self.spectrum=abs(spec[:,abs(a)>max(a)/10])
		self.fs=fs[:]
		self.t=t[:]
		self.imin=imin
		self.imax=imax