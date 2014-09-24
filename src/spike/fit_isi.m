function fit_isi(srData, desc)

    allISIs = [];
    
    nDatasets = length(srData.datasets);
    
    for k = 1:nDatasets
     
        ds = srData.datasets{k};
        stimLengthMs = ds.stim.stimLength * 1e3;
        
        for m = 1:length(ds.resp.rawSpikeTimes)
       
            rawSpikes = ds.resp.rawSpikeTimes{m};
            
            rawSpikes(rawSpikes < 0) = [];
            rawSpikes(rawSpikes > stimLengthMs) = [];
            
            isis = diff([0 rawSpikes]);
            allISIs = [allISIs isis];
            
        end
        
    end
    
    gparams = gamfit(allISIs);
    ga = gparams(1);
    gb = gparams(2);
    
    emu = expfit(allISIs);
    
    [n, xout] = hist(allISIs, 200);
    normISIs = n / sum(n);
    
    gpdf = gampdf(xout, ga, gb);
    epdf = exppdf(xout, emu);
    
    figure; hold on;
    plot(xout, normISIs, 'k-');
    plot(xout, gpdf, 'r-');
    plot(xout, epdf, 'b-');
    legend('Real', 'Gamma', 'Exp');
    if nargin > 1
        title(desc);
    end
    