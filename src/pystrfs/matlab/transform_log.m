function wholeStim = transform_log(wholeStim)

    dbnoise = 80;
    refpow = max(wholeStim(:));    
    wholeStim = max(0, 20*log10(wholeStim/refpow) + dbnoise);
    
    