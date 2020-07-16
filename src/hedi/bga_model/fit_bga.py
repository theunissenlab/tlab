from numpy import *
from pySpec import *
from pySynth import *

import sndfile
from scipy.optimize import fmin

slope=-0.6025;
intercept=-0.5979;	
data_dir='../data/'
fit_dir='../output_txt/'

import sys

def chi2all(p,f,y,v=0):
	alpha=p[3]
	beta=intercept + slope*alpha

	a=p[0]
	L=p[1]
	amp=p[2]
	s=getSpectrumAlpha(alpha,beta,a,L,amp)
	d= sum((y-s)**2)
	if v:
		print alpha,a,L,amp,d
		fi=open('out.out','w')
		for r,v in zip(y,s):
			fi.write('%lg %lg\n'%(r,v))
		fi.close()
	return d

def compute_fitall(f,tspec,a0):
	v = fmin(chi2all, a0,args=(f,tspec),disp=False,maxiter=50)
	return v


def chi2(p,f,y,alpha,v=0):
	beta=intercept + slope*alpha
	a=p[0]
	L=p[1]
	amp=p[2]
	s=getSpectrumAlpha(alpha,beta,a,L,amp)
	d= sum((y-s)**2)
	if v:
		print alpha,a,L,amp,d
		fi=open('out.out','w')
		for r,v in zip(y,s):
			fi.write('%lg %lg\n'%(r,v))
		fi.close()
	return d

def compute_fit(f,tspec,alpha,a0):
	v = fmin(chi2, a0,args=(f,tspec,alpha),disp=False,maxiter=50)
	return v


def load_file(fname):
	sfile = sndfile.open(fname,"r")
	info = sfile.get_info()
	d=sfile.readf_double(info.frames)
	sfile.close()
	spec_wav=pySpec(d)
	return spec_wav,d



class fit_procedure(object):
	def __init__(self,sex,bird,call):
		self.sex=sex
		self.bird=bird
		self.call=call
		filename='%c'%sex+'bird%dcall%d'%(bird,call)
		self.target,self.wave=load_file(data_dir+filename+'.wav')
		self.fit_text=fromfile(fit_dir+filename,sep='\n',dtype=double)
		self.fit_text=reshape(self.fit_text,(-1,9))
		self.spec=self.target.spectrum[self.target.imin:self.target.imax,:]
		self.fs=len(self.target.fs)
		self.tt=len(self.target.t)
		self.filename=filename
		self.dumpname='dump-'+self.filename

		n=len(self.fit_text)
		if n<self.tt:
			row=self.fit_text[-1,:]
	##			print len(self.fit_text),self.tt
##			print row
			while len(self.fit_text)<self.tt:
				self.fit_text=vstack((self.fit_text,row))
		
	def set_guess(self,guess):
		self.initial_guess=guess[:]
		self.current_guess=guess[:]
	def init_result(self):
		self.result=zeros((self.tt,4+len(self.initial_guess)))
	def fit(self,t):
		alpha0=self.fit_text[t,1]
		v= compute_fit(self.target.fs,self.spec[:,t],alpha0,self.current_guess)
		dv=chi2(v,self.target.fs,self.spec[:,t],alpha0)
		r=[t,dv,alpha0,slope * alpha0 + intercept]
		r.extend(v)
		st=''.join(['%.2f ']*len(r))
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
			self.current_guess=vt[4:]
			self.result[t,:]=vt	
		for t in range(tt/2-1,-1,-1):
			vt=self.fit(t)
			self.current_guess=vt[4:]
			self.result[t,:]=vt	
		self.result[self.result[:,0].argsort(),]
	def save_result(self):
		fFile=open(fit_dir + 'fit.'+self.filename,"w")
		st=''.join(['%lg ']*(4+len(self.initial_guess))) + '\n'
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



import cPickle
import os.path
sex=sys.argv[1]
bird=int(sys.argv[2])
call=int(sys.argv[3])
b0=[1034,19.7,2]
fp=fit_procedure(sex,bird,call)
# fp.set_guess(b0)
# t=49
# v=fp.fit(t)
# chi2(v[4:],0,fp.spec[:,t],v[2],v=1)

if os.path.isfile(fp.dumpname):
	fp=cPickle.load(open(fp.dumpname,'rb'))
	fp.reparse(200)
	cPickle.dump(fp,open(fp.dumpname,'wb'))
	fp.save_result()
else:
	fp.initial_fit_all(b0)
	cPickle.dump(fp,open(fp.dumpname,'wb'))
	fp.save_result()

