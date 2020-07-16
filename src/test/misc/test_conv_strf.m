function test_conv_strf
    
    nChannels = 2;
    timeLen = 100;
    strfLen = 10;
    delays = 0:(strfLen-1);

    stim1 = randn(timeLen, nChannels);
    stim2 = randn(timeLen, nChannels);
    strf = randn(nChannels, strfLen);
    
    gi1 = ones(1, timeLen);
    gi2 = ones(1, timeLen);
    
    resp1 = conv_strf(stim1, delays, strf, gi1);
    resp2 = conv_strf(stim2, delays, strf, gi2);
    mresp = [resp1 resp2];
    
    groupIndex = [gi1 gi2*2];
    stim = [stim1; stim2];
    resp = conv_strf(stim, delays, strf, groupIndex);
    
    respDiff = norm(resp - mresp)
    figure; hold on;
    plot(mresp, 'k-');
    plot(resp, 'r-');