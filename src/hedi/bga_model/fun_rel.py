from numpy import *
import pyBGA



fFile=open("stf.fun.dat")


data={}
freqs={500,750,1000} #Hz
fband=25

for f in freqs:
	data[f]=[]

for l in fFile:
	s=l.split()
	alpha=float(s[0])
	beta=float(s[1])
	f0=float(s[2])
	for freq in freqs:
		f1=freq-fband
		f2=freq+fband
		if f0>f1 and f0<f2:
			data[freq].append([alpha,beta])
fFile.close()

#slope=-0.6025;
#intercept=-0.5979;

#print data
samplerate=44100.0
gamma=24000
tmax=0.5
size=floor(samplerate*tmax)
da=zeros(size)

for freq in data.keys():
	fFile=open("%f-power.dat"%freq,"w")
	for p in data[freq]:
		alpha=p[0]
		beta=p[1]
		pyBGA.pyBGAs(da,gamma,alpha,beta)
		fFile.write("%lf %lf %lf\n" %(alpha,beta,max(abs(da-mean(da[:size/2])))))
	fFile.close()



