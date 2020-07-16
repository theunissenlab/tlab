from numpy import mean, isnan, diag, std, var, sqrt, corrcoef, sum, array
from scipy.stats import linregress

__all__ = ['zscore','dprime']

def zscore(input,*args):
    try:
        return mean(input,*args)/std(input,*args)
    except ZeroDivisionError:
        return 0
    
def dprime(input_plus,input_minus,**kwargs):
    '''
    Compute d-prime statistic from two sample distributions.
    '''
    kwargs.setdefault('ddof',1)
    return dprime_param(len(input_plus),mean(input_plus),std(input_plus,**kwargs),
                        len(input_minus),mean(input_minus),std(input_minus,**kwargs))
    
def dprime_param(n1,mean1,se1,n2,mean2,se2):
    '''
    Compute d-prime statistic from the means and variances
    of two distributions.
    '''
    
    denom = se1**2 + se2**2
    if denom == 0:
        return 0
    else:
        return (mean1-mean2)/sqrt(denom)

def cohen_d(input_plus,input_minus,**kwargs):
    '''
    Compute Cohen's d metric, which is similar to the d-prime
    Ported 2/25/11 by RCM from 
    http://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d
    '''
    
    kwargs.setdefault('ddof',1)
    
    n1 = len(input_plus)
    mean1 = mean(input_plus)
    se1 = std(input_plus,**kwargs)
    
    n2 = len(input_minus)
    mean2 = mean(input_minus)
    se2 = std(input_minus,**kwargs)
    
    return cohen_d_param(n1,mean1,se1,n2,mean2,se2)
    
def cohen_d_param(n1,mean1,se1,n2,mean2,se2):
    '''
    Compute Cohen's d metric from parameters
    Ported 2/25/11 by RCM from 
    http://en.wikipedia.org/wiki/Effect_size#Cohen.27s_d
    '''
    
    # Denominator
    s1n = (n1 - 1) * s1**2
    s2n = (n2 - 2) * s2**2
    s = sqrt((s1n + s2n) / (n1 + n2))
    
    return (mean1 - mean2) / s
    
def cc(*args):
    ccs = corrcoef([a for a in args])
    cc = ccs[0,1]
    if isnan(cc):
        if all(diag(ccs)==0):
            return 1
        else:
            return 0
    else:
        return cc

def linregress(*args,**kwargs):
    if kwargs.pop('intercept',True):
        return linregress(*args)
    else:
        x = array(args[0])
        y = array(args[1])
        return sum(x*y)/sum(x*x)
    
