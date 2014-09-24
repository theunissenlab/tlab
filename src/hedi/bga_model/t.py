from numpy import *

import test

slope=-0.6025;
intercept=-0.5979;

samplerate=44100.0
tmax=0.5
size=floor(samplerate*tmax)
data=zeros(size)
alpha=-0.401
beta=intercept + slope * alpha
gamma=24000

test.smBGA(data,gamma,alpha,beta)


for inter in arange(-.65,-0.55,0.001):
	beta=inter + slope * alpha
	test.smBGA(data,gamma,alpha,beta)
	print inter,max(abs(data-mean(data[:size/2]))),argmax(abs(fft.fft(data[:size/2]-mean(data))))