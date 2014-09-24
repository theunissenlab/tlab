'''
Welch t-test ported from
http://mail.scipy.org/pipermail/scipy-user/2007-August/013311.html
2/25/11
'''

import numpy as n
import scipy

def welchs_approximate_ttest(n1, mean1, sem1, \
                             n2, mean2, sem2, alpha):
    '''Welch''s approximate t-test for the difference of two means of
    heteroscedasctic populations.
    
    Implemented from Biometry, Sokal and Rohlf, 3rd ed., 1995, Box 13.4
    
    :Parameters:
        n1 : int
            number of variates in sample 1
        n2 : int
            number of variates in sample 2
        mean1 : float
            mean of sample 1
        mean2 : float
            mean of sample 2
        sem1 : float
            standard error of mean1
        sem2 : float
            standard error of mean2
        alpha : float
            desired level of significance of test
    
    :Returns:
        significant : bool
            True if means are significantly different, else False
        t_s_prime : float
            t_prime value for difference of means
        t_alpha_prime : float
            critical value of t_prime at given level of significance
    
    Copyright (c) 2007, Angus McMorland
    
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the University of Auckland, New Zealand nor
        the names of its contributors may be used to endorse or promote
        products derived from this software without specific prior written
        permission.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL,EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.'''
    
    svm1 = sem1**2 * n1
    svm2 = sem2**2 * n2
    t_s_prime = (mean1 - mean2)/n.sqrt(svm1/n1+svm2/n2)
    t_alpha_df1 = scipy.stats.t.ppf(1-alpha/2, n1 - 1)
    t_alpha_df2 = scipy.stats.t.ppf(1-alpha/2, n2 - 1)
    t_alpha_prime = (t_alpha_df1 * sem1**2 + t_alpha_df2 * sem2**2) / \
                    (sem1**2 + sem2**2)
    return abs(t_s_prime) > t_alpha_prime, t_s_prime, t_alpha_prime

def welch_fischer_stars(*args):
    '''
    Quickie method to return Fischer's significance stars from Welch
    ttest.
    ''    means p > 0.05
    '*'   means p < 0.05
    '**'  means p < 0.01
    '***' means p < 0.001
    '''
    args = list(args)
    args.append(.05)
    significant,t_sample,t_critical = welchs_approximate_ttest(*args)
    if significant:
        t_args = [args[0],args[2],args[3],args[5]]
        if t_alpha_df(t_sample,0.01,*t_args):
            if t_alpha_df(t_sample,0.001,*t_args):
                return '***'
            else:
                return '**'
        else:
            return '*'
    else:
        return ''
    
def t_alpha_df(t_sample,alpha,n1,sem1,n2,sem2):
    '''
    Utility function for welch_fischer_stars. Gets the cutoff at a given
    significance.
    '''
    
    t_alpha_df1 = scipy.stats.t.ppf(1-alpha/2, n1 - 1)
    t_alpha_df2 = scipy.stats.t.ppf(1-alpha/2, n2 - 1)
    t_alpha_prime = (t_alpha_df1 * sem1**2 + t_alpha_df2 * sem2**2) / \
                    (sem1**2 + sem2**2)
    return abs(t_sample) >= t_alpha_prime

def welch_t(x,y,alpha=0.05):
    n1 = len(x)
    n2 = len(y)
    mean1 = n.mean(x)
    mean2 = n.mean(y)
    sem1 = n.std(x,ddof=1)/n.sqrt(n1)
    sem2 = n.std(y,ddof=1)/n.sqrt(n2)
    return welchs_approximate_ttest(n1,mean1,sem1,n2,mean2,sem2,alpha)