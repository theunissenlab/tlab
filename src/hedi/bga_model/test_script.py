from numpy import *
import pyBGA

slope=-0.6025;
intercept=-0.5979;	
alpha=-0.41769

samplerate=44100.0
tmax=0.5
size=floor(samplerate*tmax)
data=zeros(size)
beta=intercept + slope * alpha
gamma=24000

pyBGA.pyBGAs(data,gamma,alpha,beta)

for x in data:
	print x
