from .resample import jackknife_subsample

def jaccknife_estimate(samples,func,n_delete=1):
    subsamples = jackknife_subsample(samples,n_delete)
    