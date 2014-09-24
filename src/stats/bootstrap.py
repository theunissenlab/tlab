from numpy import (asarray, mean, sqrt, std, isfinite)
from numpy.random import normal
from scipy.stats import linregress, ss, t
from matplotlib.pyplot import gca, gcf
from .resample import subsample

test_data = normal(size=8)
#test_data = range(4)

def jackknife_bias_correct(pairs,confidence=None,return_all=False,
                           nan_remove=True,return_raw=False):
    '''
    Return jackknife-bias-corrected estimate from estimate-nsamples pairs
    Pairs can be either a list of tuples, or a 2 x nestimates array.
    If 'confidence' is between 0 and 1, return the mean with lower and upper
    bounds at -/+ the confidence interval.
    If 'confidence' is None, return the mean and standard error.
    If 'return_all' is True, return the mean, standard error, number of points,
    and confidence interval size.
    '''
    
    data = asarray(pairs)
    
    if nan_remove:
        data = data[isfinite(data)[:,0],:]
    
    y = data[:,0]
    x = 1./data[:,1]
    n = len(x)
    
    # Compute linear regression and standard error of intercept
    (slope,intercept,r,p,slope_se) = linregress(x,y)
    intercept_se = slope_se * sqrt(ss(x)/n)
    
    # Return mean and SE if no value is specified:
    if confidence is None:
        if return_all:
            if return_raw:
                np = data[:,1]
                max_n = max(np)
                raw_mean = mean(data[np==max_n,0])
                return intercept, intercept_se, n, raw_mean
            else:
                return intercept, intercept_se, n
        else:
            return intercept, intercept_se
    
    # Otherwise return intercept with confidence
    else:
        t_int = t._ppf((1+confidence)/2,n-2)
        intercept_int = t_int * intercept_se
        
        if return_all:
            if return_raw:
                np = data[:,1]
                max_n = max(np)
                raw_mean = mean(data[np==max_n,0])
                return intercept, intercept_se, n, intercept_int, raw_mean
            else:
                return intercept, intercept_se, n, intercept_int
        else:
            return intercept, intercept - intercept_int, intercept + intercept_int

def jk_plot(pairs,**kwargs):
    '''
    Plots a scatterplot of the pairs with fit line and error bars
    '''
    
    # Plot into current axes if none specified
    axes = kwargs.pop('axes',None)
    if axes is None:
        axes = gca()
    
    # Split pairs into vectors of cc estimates and n samples for plotting
    array_pairs = asarray(pairs)
    ccs = array_pairs[:,0]
    n_samples = array_pairs[:,1]
    
    # Plot array pairs
    h1 = axes.plot(1/n_samples,ccs,'b+',**kwargs)
    
    # Comppute jackknifed estimate of the cc with confidence interval
    jk_cc,jk_cc_down, jk_cc_up = jackknife_bias_correct(pairs,confidence=.95)
    
    # Compute means and errorbars for fit line
    sample_sizes = sorted(set(n_samples),reverse=True)
    means = [jk_cc]
    errors = [jk_cc-jk_cc_down]
    for sample_size in sample_sizes:
        this_set = ccs[n_samples==sample_size]
        n_points = len(this_set)
        means.append(mean(this_set))
        errors.append(std(this_set)/sqrt(n_points))

    # Plot fit line with error bars
    fit_xdata = [0] + [1./n for n in sample_sizes]
    h2 = axes.errorbar(fit_xdata,means,errors,fmt='k-')
    
    # Draw errorbars at weight 2
    h2[0].set_lw(2)
    for h in h2[1]:
        h.set_markeredgewidth(2)
    h2[2][0].set_lw(2)
    
    # Set x axis limits to good display range (range + 5% on each end)
    x_pad = kwargs.pop('x_pad',.05)
    x_range = 1./min(n_samples)
    x_min = -x_pad * x_range
    x_max = (1 + x_pad) * x_range
    axes.set_xlim(x_min,x_max)
    
    # Label x axis
    axes.set_xlabel('Number of samples')
    axes.set_xticks(fit_xdata)
    axes.set_xticklabels(['Inf']+['%d' % n for n in sample_sizes])
    
    # Label y axis
    axes.set_ylabel('CC')
    
    # Set figure title
    axes.set_title('Jackknife-bias-corrected estimate of masker cc')
    
    return h1,h2
    
def jk_hist(pairs,**kwargs):
    '''
    Plots a set of histograms of the pairs by number of samples
    '''
    
    fig = kwargs.pop('fig',None)
    if fig is None:
        fig = gcf()

def test(func=mean):
    subs = subsample(test_data)
    pairs = [(func(s),len(s)) for s in subs]
    return pairs