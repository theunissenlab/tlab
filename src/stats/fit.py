'''
Copied wholesale from the SciPy cookbook 10/29/10
RCM
http://www.scipy.org/Cookbook/FittingData#head-11870c917b101bb0d4b34900a0da1b7deb613bf7
'''

from numpy import *
from scipy import optimize

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
    
def gabor(height, center_x, center_y, width_x, width_y, lambda_, theta, psi):
    return lambda x,y: height * \
                exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)\
                * cos(2 * pi * (x * cos(theta) + y * sin(theta))/lambda_ + psi)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def freqs(data):
    f_data = fft.fft2(data)
    mag = f_data * f_data.conj()
    f_max = amax(f_data)
    x_max,y_max = nonzero(f_data==f_max)
    lambda_ = sqrt(x_max**2 + y_max**2)
    theta = arctan(y_max/x_max)
    return lambda_, theta

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

def fitgabor(data):
    params = list(moments(data)) + list(freqs(data)) + [0]
    errorfunction = lambda p: ravel(gabor(*p)(*indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p