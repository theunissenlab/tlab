function A = make_toeplitz(stim, ndelays)

    stim_sz = size(stim);
    stim_len = stim_sz(1);
    num_freqs = stim_sz(2);
    
    A = zeros(stim_len, num_freqs*ndelays);
    %A_sz = size(A)
    
    for t = 1:stim_len
        z = zeros(ndelays, num_freqs);
        %z_sz = size(z)
       
        tstart = max(1, t - ndelays + 1);
        zlen = t - tstart + 1;
        zstart = ndelays - zlen + 1;
        %z2_sz = size(z(zstart:end, :))
        %stim_sz = size(stim(tstart:t, :))
        z(zstart:end, :) = stim(tstart:t, :);
        
        A(t, :) = reshape(z, 1, num_freqs*ndelays);        
    end
    
end