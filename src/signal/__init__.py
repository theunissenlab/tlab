from numpy import array, asarray, arange, sum, ndarray, any
from scipy.stats import rv_continuous, norm
#from scipy.stats.distributions import rv_frozen, rv_generic
from scipy.stats.distributions import rv_frozen, rv_continuous
from inspect import isfunction

__all__ = ['Window','FuncWindow','DistWindow','window',
           'gaussian_window','wrap_to_integer']

class Window(object):
    '''
    Signal processing window
    '''
    
    def __init__(self,win_length,**kwargs):
        '''
        __init__(self,win_length,normalize=False)
        Currently, this only supports conversion of distributions
        Conversion from numpy windows (e.g. Hanning, Hamming) not yet supported
        '''
        self._len = win_length
        
        if kwargs.pop('normalize',False):
            self.normalize()
        else:
            self.sf = 1

    def __len__(self):
        return self._len
    
    def __array__(self):
        return self._raw()/self.sf
    
    def normalize(self):
        self.sf = sum(self._raw())
        return self
        
    def raw(self):
        self.sf = 1
    
class FuncWindow(Window):
    '''
    Window defined from a function like Numpy's function windows
    '''
    
    def __init__(self,func,*args,**kwargs):
        '''
        __init__(self,func,win_length,normalize=False)
        
        'func' can be either a string, or a function.
            If func is a string, it is evaluated in Numpy's base 
            namespace; thus entries like 'blackman' or 'hanning'
            will return the respective window functions.
            If func is a function, the function is presumed to
            work like one of Numpy's window functions. There is
            currently no provision to pass any additional arguments
            to the function.
        
        'win_length' is the number of points in the window.
        
        'normalize', if True, scales the coefficients so that the sum
            of the window is 1.
        '''
        
        if isinstance(func,str):
            import numpy as n
            # Evaluate string arguments in base Numpy namespace.
            # Hey, it's your funeral. Don't pass bad ideas here.
            try:
                self.func = getattr(n,func.lower())
            except AttributeError:
                raise ValueError,"Numpy has no function named '%s'" % (func)
        elif isfunction(func):
            self.func = func
        else:
            raise TypeError,"Input argument 'func' must be of type 'string' or 'function'"
        
        Window.__init__(self,*args,**kwargs)
    
    def _raw(self):
        return self.func(self._len)
    
    def __array__(self):
        return self._raw()/self.sf

class DistWindow(Window):
    '''
    Window defined from probablity distribution
    Defined to have its "location" parameter
    at the center of the window unless otherwise specified.
    '''
    
    def __init__(self,dist,*args,**kwargs):
        '''
        __init__(self,dist,win_length,loc,scale,**kwargs)
        '''
        # Set location parameter for the distribution so that loc=0
        # puts loc at the middle of the window
        Window.__init__(self,*args,**kwargs)
        
        if len(args) > 1:
            kwargs.update({'loc':args[1]})
        if len(args) > 2:
            kwargs.update({'scale':args[2]})
        
        kwargs.update({'loc': 0.5 * (self._len - 1) + kwargs.get('loc',0)})
        
        if isinstance(dist,rv_frozen):
            # Must specify the loc and scale parameters for the distribution
            # after binding
            raise
        
        self.dist = dist(*args,**kwargs)
        
    def _raw(self):
        return self.dist.pdf(arange(len(self)))
    
    def normalize(self):
        self.sf = sum(self._raw())
        
    def __array__(self):
        return self._raw()/self.sf
    
def window(win_type,*args,**kwargs):
    #if isinstance(win_type,rv_generic):
    if isinstance(win_type,rv_continuous):
        return DistWindow(win_type,*args,**kwargs)
    elif isinstance(win_type,str):
        return FuncWindow(win_type,*args,**kwargs)
    elif isfunction(win_type):
        return FuncWindow(win_type,*args,**kwargs)

def gaussian_window(win_length,n_sd=6):
    sd = float(win_length)/n_sd
    return DistWindow(norm,win_length,loc=0,scale=sd)

def wrap_to_integer(arr,dtype):
    nbytes = array([],dtype=dtype).itemsize
    intrange = 2 ** (nbytes * 8)
    if dtype[0] == 'u':
        minint = 0
        maxint = intrange - 1
    else:
        minint = -intrange/2
        maxint = -minint - 1
    
    wrap = 1
    while any(wrap):
        # Wrap max values
        wrap = arr > maxint
        arr -= intrange * wrap
    
    wrap = 1
    while any(wrap):
        # Wrap min values
        wrap = arr < minint
        arr += intrange * wrap
        
    return asarray(arr,dtype=dtype)