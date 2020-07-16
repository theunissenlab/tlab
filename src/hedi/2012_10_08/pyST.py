import sys
import warnings
warnings.filterwarnings("ignore")
from numpy import *
from numpy import fft
import scipy
from scipy import signal,stats
import sndfile
hamming = scipy.signal.hamming


class pyCall:
    def __init__(self,filename,divisor,fftwindow,dB,LowPass):
        self.fname=filename
	self.sfile = sndfile.open(filename,"r")
	self.info = self.sfile.get_info()
	ch = self.info.channels
        self.samplerate=self.info.samplerate
        self.frames=self.info.frames
        self.r=self.sfile.read_double(self.frames)
        self.sfile.close()
	self.rescale()
        self.fftwindow=fftwindow
        self.divisor=divisor
        self.dB=dB
        self.LowPass=LowPass
    def rescale(self):
	m = max(abs(self.r))
	self.r  = self.r/m
    def computespectrum(self):
	ch = self.info.channels
	samplerate = self.info.samplerate
        fftwindow=self.fftwindow
        LowPass=self.LowPass
        divisor=self.divisor
	Hz = min(fftwindow*LowPass/samplerate,fftwindow/2-1)        
	self.env = zeros((self.frames/divisor))
       	for k in range(self.frames/divisor):
	    a = self.r[(ch*k*divisor):(ch*k*divisor+ch*fftwindow):ch]
	    a = a -mean(a)
	    if a.shape !=(fftwindow,):
		a.resize((fftwindow,))
	    a = a*hamming(fftwindow)
	    y = abs(fft.fft(a,fftwindow))**2
	    y = y[:fftwindow/2]
	    self.env[k] = sum(y[Hz:len(y)])
    def save_wave(self,filename,start=0,end=-1):
        outinfo = sndfile.SF_INFO(samplerate=44100, 
                                  channels = 1, 
                                  format = (sndfile.SF_FORMAT_WAV|
					 sndfile.SF_FORMAT_PCM_16), 
                                  sections = 1, 
                                  seekable = 1)
        file = sndfile.open(filename,mode="w",info=outinfo)
	file.write_double(self.r[::self.info.channels])
	file.close()
    def getcall(self):
	self.call = []
        divisor=self.divisor
        fftwindos=self.fftwindow
        dB=self.dB
  	m = max(self.env)
	idmax = self.env.argmax()
	start = 0
	star1 = 0
	star2 = 0
	end = 0
	x = m
	i = idmax
	while x>m*dB and i<self.frames/divisor:
	    x = self.env[i]
	    i += 1
	end = i
	x = m
	i = idmax
	while x>m*dB and i >=0:
	    x = self.env[i]
	    i += -1
	start = i
	if start<0:
	    start =0
	    star1 = 1
	if end > self.frames/divisor:
	    end = self.frames/divisor
	    star2 = 1
	self.call = (start*divisor,end*divisor)
	return m,idmax*divisor,self.call,star1,star2
   
    def extract(self,dir):
        self.computespectrum()
        self.getcall();
        self.save_wave(dir+'/'+self.fname,self.call[0],self.call[1])

import sys


fname=sys.argv[1]
fftwindow=1024
divisor=10
dB=0.1
LowPass=500
c=pyCall(fname,divisor,fftwindow,dB,LowPass)
c.extract('../data/')
