import warnings
warnings.filterwarnings("ignore")

from numpy import *
from pySpec import *
from pySynth import *
from scipy.interpolate import interp1d
import sndfile
from scipy.optimize import fmin

slope=-0.6025;
intercept=-0.5979;	
data_dir='../data/'
fit_dir='../output_txt/'

import sys

def chi2all(p,f,y,v=0):
	n=len(p)
	alpha=p[0]
	beta=intercept + slope*alpha
	c=zeros(2+(n-1)/2,dtype=double)
	d=zeros(2+(n-1)/2,dtype=double)
	c[1:-1]=p[1:(n-1)/2+1]
	d[1:-1]=p[((n-1)/2+1):]
	c[0]=0
	c[-1]=10000
	d[0]=0
	d[-1]=0
	
	s=getSpectrumCubicFilter(alpha,beta,c,d)
	dd= sum((y-s)**2)
	##print dd,alpha,beta,c,d
	if v>0:
		fi=open('out.out','w')
		fil_func=interp1d(c,d)
		e = fil_func(f)
		for ff,r,v,ss in zip(f,y,s,e):
			fi.write('%lg %lg %lg %lg\n'%(ff,r,v,ss))
		fi.close()
	return dd

def compute_fitall(f,tspec,a0):
	v = fmin(chi2all, a0,args=(f,tspec),disp=False,maxiter=100)
	return v

def load_file(fname):
	sfile = sndfile.open(fname,"r")
	info = sfile.get_info()
	d=sfile.readf_double(info.frames)
	sfile.close()
	spec_wav=pySpec(d)
	return spec_wav,d



class fit_procedure(object):
	def __init__(self,fname):
		self.filename=fname.split('/')[-1][:-4]
		self.target,self.wave=load_file(data_dir+self.filename+'.wav')
		self.spec=self.target.spectrum[self.target.imin:self.target.imax,:]
		self.fs=len(self.target.fs)
		self.tt=len(self.target.t)
		self.dumpname='dump-'+self.filename
	def set_guess(self,guess):
		self.initial_guess=guess[:]
		self.current_guess=guess[:]
	def init_result(self):
		self.result=zeros((self.tt,2+len(self.initial_guess)))		
	def fit(self,t):
		v= compute_fitall(self.target.fs,self.spec[:,t],self.current_guess)
		dv=chi2all(v,self.target.fs,self.spec[:,t],1)
		r=[t,dv]
		r.extend(v)
		st=''.join(['%.5f ']*len(r))
		print	st%tuple(r)
		return array(r)
	def refit(self,t,all_param):
		v= compute_fitall(self.target.fs,self.spec[:,t],all_param)
		dv=chi2all(v,self.target.fs,self.spec[:,t])
		alpha0=v[-1]
		if alpha0>-0.40:
			alpha0=-0.40
		v=v[:-1]
		r=[t,dv,alpha0,slope * alpha0 + intercept]
		r.extend(v)
		st=''.join(['%.2f ']*len(r))
		print	st%tuple(r)
		return array(r)
	def initial_fit_all(self,initial_guess):
		self.set_guess(initial_guess)
		self.init_result()
		tt=self.tt
		for t in range(tt/2,tt):
			vt=self.fit(t)
			self.current_guess=vt[2:]
			self.result[t,:]=vt
		self.current_guess=initial_guess
		for t in range(tt/2-1,-1,-1):
			vt=self.fit(t)
			self.current_guess=vt[2:]
			self.result[t,:]=vt	
		self.result[self.result[:,0].argsort(),]
	def save_result(self,outdir):		
		fFile=open(outdir + self.filename +'.fit',"w")
		st=''.join(['%lg ']*(len(self.result[0,:]))) + '\n'
		for x in self.result:
			fFile.write(st%tuple(x))
		fFile.close()
	def reparse(self,thr):
		bound=5
		w=where(self.result[:,1]>thr)
		for worst in w[0]:
			nbound=worst-bound 
			if nbound<0: nbound=0
			ubound=worst+bound
			if ubound>=len(self.result[:,1]):
				ubound=len(self.result[:,1])-1
			best_around=argmin(self.result[worst-bound:worst+bound,1])
			new_guess=self.result[best_around+worst-bound,:]
			new_param=concatenate((new_guess[4:],array([new_guess[2]])))
			new_r=self.refit(worst,new_param)
			print self.result[worst,:]
			if new_r[1]<self.result[worst,1]:
				self.result[worst,:]=new_r



if __name__ == '__main__':

	import cPickle
	import os.path
	fname=sys.argv[1]
	a0=float(sys.argv[2])
	guess=[a0,100,2000,4500,6000,0.01,0.05,1.5,0.3]
	fp=fit_procedure(fname)
	fp.initial_fit_all(guess)
	fp.save_result('../python_output/')