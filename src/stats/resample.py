from numpy import (arange, asarray, empty, unpackbits, sum, bool_, hstack)

def sub_masks(n_points,n_samples):
    '''
    Return all possible choices of n_samples from n_points as boolean masks
    for indexing use.
    n_points must be 8 or less because we use numpy.unpackbits
    '''
    
    if (n_points < 1):
        raise ValueError,'n_points must be at least 1'
    
    # We use numpy.unpackbits for each octet
    elif n_points <= 8:
        
        # Create a column vector with all numbers from 0 to 2^n_points
        n_masks = 2**n_points
        num_masks = arange(n_masks,dtype='uint8').reshape(n_masks,1)
        
        # Transform column vector into binary masks with n_points bits
        bin_masks = unpackbits(num_masks,axis=1)[:,-n_points:]
        
        # Find the masks that have n_samples ones
        good_masks = sum(bin_masks,axis=1) == n_samples
        
        # Return these masks as numpy boolean arrays for indexing
        return bool_(bin_masks[good_masks,:])
    
    # For more than 8 bits, recurse
    else:
        #raise ValueError,'n_points cannot be more than 8'
        # Truncate one octet
        n_points_rem = n_points - 8
        
        # Compute how many samples we can take from the remainder
        min_samples_rem = max(n_samples - 8, 0)
        max_samples_rem = min(n_samples, n_points_rem)
        
        # Initialize output
        output = []
        
        # Loop over allowed valuse
        for n_samples_rem in range(min_samples_rem, max_samples_rem + 1):
            
            # Compute how many samples to take from the octet
            n_samples_oct = n_samples - n_samples_rem
            
            # Compute sample masks for the octet
            oct_masks = sub_masks(8,n_samples_oct)
            
            # Compute sample masks for the remainder
            rem_masks = sub_masks(n_points_rem,n_samples_rem)
            
            # Combine the sets    
            output.extend([hstack((rem,oct)) for rem in rem_masks for oct in oct_masks])
        return output

def jackknife_subsample(sample_data,n_delete,**kwargs):
    n_samples = len(sample_data) - n_delete
    return subsample_n(sample_data,n_samples,**kwargs)

def subsample_n(sample_data,n_samples,**kwargs):
    '''
    Compute all leave-n jackknife subsamples from the data
    '''
    n_points = len(sample_data)
    
    # Copy things into array -- kludgy but needs to be done
    sample_array = empty(len(sample_data),dtype='object')
    for idx,obj in enumerate(sample_data):
        sample_array[idx] = obj
        
    # Read out masks
    for mask in sub_masks(n_points,n_samples):
        #yield list(sample_array[mask])
        yield sample_array[mask]

def subsample(sample_data,**kwargs):
    '''
    Return a full set of subsamples from data
    '''
    
    # Total size of sample set
    n_points = len(sample_data)
    
    # Maximum number of samples to remove -- defaults to n_points-1,
    # computing all subsample sizes from 1 to n_points
    max_n_remove = kwargs.pop('max_n_remove',n_points-1)
    
    # Size of the smallest sample
    min_subsample_size = n_points - max_n_remove
    
    # Return list of all samples in range
    return list(sub for n_samples in range(min_subsample_size,n_points+1)\
                    for sub in subsample_n(sample_data,n_samples,**kwargs))