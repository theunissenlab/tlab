from numpy import *
from pySpec import *
from scipy.interpolate import interp1d
import pyBGA
samplerate=44100.0
tmax=0.5
gamma=24000
c=340


def getSpectrumAlpha(alpha,beta,a,L,amp):
	size=floor(samplerate*tmax)
	data=zeros(size)
	pyBGA.pyBGAs(data,gamma,alpha,beta)
	S=pySpec(data)

	fft,fc=shape(S.whole_spec)
	fw=c/(4*L)*1000
	spec=S.whole_spec[S.imin:S.imax,fc/2]
##	print S.fs,fw
##	filter=amp/(1+a*exp(2*pi*S.fs*2*L/c/1000*1.j)) 
	filter= amp * exp(-(S.fs-fw)**2/(2*a**2))
	specf=spec*filter
	return abs(specf)##,abs(spec),abs(filter)


def getSpectrumCubicFilter(alpha,beta,x,y):
	# x and y must have the same size 
	size=floor(samplerate*tmax)
	data=zeros(size)
	pyBGA.pyBGAs(data,gamma,alpha,beta)
	S=pySpec(data)

	fft,fc=shape(S.whole_spec)
	#print x,y
	spec=S.whole_spec[S.imin:S.imax,fc/2]
	fil_func=interp1d(x,y)
	#print fil_func(S.fs)
	specf=spec*fil_func(S.fs)

	#print specf,shape(specf)
	return abs(specf)##,abs(spec),abs(filter)

if __name__ == '__main__':
	slope=-0.6025;
	intercept=-0.5979;	
	alpha=-0.41769
	a,b,c= getSpectrumAlpha(alpha,intercept+slope*alpha,0.8	,2000,19,2)
#	for x,y,z in zip(a,b,c):
#		print x,y,z