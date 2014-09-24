from libc.stdlib cimport free
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

cdef extern from "numpy/arrayobject.h":
	ctypedef int intp
	ctypedef extern class numpy.ndarray [object PyArrayObject]:
		cdef char *data
		cdef int nd
		cdef intp *dimensions
		cdef intp *strides
		cdef int flags

cdef extern from "src/bga.h":
	void smBGAs(double * data, int size,double gamma, double alpha, double beta)
	void smBGAsl(double * data, int size,double gamma, double * alpha, double * talpha, int nap,double * beta, double *tbeta, int nbp,double * filter, double * freq, int nfp)
cdef extern from "src/spectrum.h":
	double * GaussianSpectrumS(double * input, int inputsize,int increment, int winLength,int * fc, int *fft)
	void InvertAndAddS(double * spec,int frameCount,int winLength,int increment,double * data)


def pyGaussianSpectrum(ndarray input,int increment,int winLength):
	cdef double *inp = <double *>input.data
	cdef int in_size = input.dimensions[0]
	cdef int fc
	cdef int fft
	cdef double * res=GaussianSpectrumS(inp, in_size,increment,winLength,&fc,&fft)
	cdef np.ndarray h = np.zeros([fft, fc], dtype=complex)
	for i in range(fc):
		for j in range(fft):
			h[j,i]=res[2*(j+fft*i)]+ res[2*(j+fft*i)+1]*1.j 
	free(res)
	return h

def pyInvertAndAdd(ndarray spec,int increment):
	[fs,fc]=np.size(spec)	
	frameCount=fc
	winLength=fs
	cdef np.ndarray fspec = np.zeros(2*fs*fc, dtype=np.double)
	for i in range(fc):
		for j in range(fs):
			fspec[2*(j+fs*i)]=np.real(spec[j,i])
			fspec[2*(j+fs*i)+1]=np.imag(spec[j,i])
	cdef double * sp=<double*>fspec.data
	cdef np.ndarray wav = np.zeros(frameCount*increment,dtype=complex)
	cdef double * data=<double *>wav.data
	InvertAndAddS(sp,frameCount,winLength,increment,data)


def pyBGAs(ndarray a, double gamma, double alpha, double beta):
	cdef double *p = <double *>a.data
	cdef int size = a.dimensions[0]
	smBGAs(p,size,gamma,alpha,beta)

def pyBGAsl(ndarray d, double gamma, ndarray alpha, ndarray talpha, ndarray beta,ndarray tbeta,ndarray filter, ndarray freq):
	cdef double *p = <double *>d.data
	cdef int size = d.dimensions[0]
	cdef double *a = <double *>alpha.data
	cdef double *ta = <double *>talpha.data
	cdef int sizea = alpha.dimensions[0]
	cdef int sizeta = talpha.dimensions[0]
	if (sizea!=sizeta):
		raise 'SizeError'
	cdef double *b = <double *>beta.data
	cdef double *tb = <double *>tbeta.data
	cdef int sizeb = beta.dimensions[0]
	cdef int sizetb = beta.dimensions[0]

	if (sizeb!=sizetb):
		raise 'SizeError'
	cdef double *f = <double *>filter.data
	cdef double *tf = <double *>freq.data
	cdef int sizef = filter.dimensions[0]
	cdef int sizetf = freq.dimensions[0]

	if (sizeb!=sizetb):
		raise 'SizeError'
	smBGAsl(p,size,gamma,a,ta,sizea,b,tb,sizeb,f,tf,sizef)
