function fitStrf = testPP(strf)

    duration = 100;
    binSize = 0.001;
    tpl = 0:binSize:duration;

    %generate stimulus
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen);
	
    %generate fitting structure       
    pp = struct;
    %pp.params = [0.75 0.4]*2;
    pp.params = strf;
    nparams = length(strf);
    %pp.params = randn(1, nparams)*1e-1;
    pp.duration = duration;
    pp.binSize = binSize;
    pp.stim = stim;    
    pp.srFilt = [-0.3];
    pp.rateFunc = @(t, p) glmRate(t, p, pp);
    
    fprintf('Generating event times...\n');
    [eventTimes, rate] = simPP(pp);
    %{
    figure; hold on;
    plot(tpl, stim, 'b-');
    plot(eventTimes, zeros(size(eventTimes)), 'r.');
    axis tight;
    %}
    %fit the model with the data
    initialGuess = randn(1, nparams);    
    fitParams = fitPP(pp, eventTimes, initialGuess);
    
    actualParams = pp.params
    fitParams
    pdiff = norm(actualParams - fitParams)
    
    ppFit = pp;
    ppFit.params = fitParams;
    ppFit.rateFunc = @(t, p) glmRate(t, p, ppFit);
    
    fitStrf = fitParams;
    
    %generate psths from each model to compare          
    %{
    [etimes, realRate] = simPP(pp);
    [etimes, fitRate] = simPP(ppFit);
    
    realRate = realRate / max(realRate);
    fitRate = fitRate / max(fitRate);
    
    figure; hold on;
    plot(realRate, 'k-');
    plot(fitRate, 'r-');
    legend('Real', 'Fit');
    axis tight;
    %}
    