function testSim()

    duration = 5;
    binSize = 0.001;
    tpl = 0:binSize:duration;

    %generate stimulus
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen);
	
    %generate fitting structure       
    pp = struct;
    pp.params = [0.75 0.4]*1.75;
    pp.duration = duration;
    pp.binSize = binSize;
    pp.stim = stim;    
    pp.rateFunc = @(t, p) glmRate(t, p, pp);
    
    fprintf('Generating event times...\n');
    [eventTimes, rate] = simPP(pp);
    
    %generate a bunch of spike trains and see if they
    %match up to the rate
    numTrials = 1000;
    psth = zeros(1, round(duration/binSize)+1);
    for k = 1:numTrials
      
      [etimes, rate] = simPP(pp);
      eindx = round(etimes / binSize) + 1;
      psth(eindx) = psth(eindx) + 1;       
      
    end
      
    psth = psth / numTrials;    
    rate = rate / max(rate);
    
    figure; hold on;   
    plot(rate, 'k-');
    plot(psth, 'r-');
    axis tight;
    
    