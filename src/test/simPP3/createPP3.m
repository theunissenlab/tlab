function pp = createPP3(strf, sr, duration, binSize, stim)

    pp = struct;
    pp.strf = strf;
    pp.sr = sr;
    pp.duration = duration;
    pp.binSize = binSize;
    pp.sampleRate = round(1 / binSize);
    pp.stim = stim;    
    pp.nl = @(x) exp(x);
    pp.dnl = @(x) exp(x);
