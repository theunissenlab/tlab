function testSim3()

    duration = 100;
    binSize = 0.001;
    tpl = 0:binSize:duration;

    %generate stimulus
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen);
	
    %generate fitting structure           
    strf = [0.75 0.4]*2;
    sr = [-3];    
    pp = createPP(strf, sr, duration, binSize, stim);
    
    fprintf('Generating event times...\n');
    [eventTimes, rate, stimCurrent, srCurrent] = simPP2(pp);
    
    figure; hold on;
    
    subplot(3, 1, 1);
    plot(tpl, stim, 'k-');
    title('Stimulus');
    
    subplot(3, 1, 2); hold on;
    plot(tpl, stimCurrent, 'b-');
    plot(tpl, srCurrent, 'r-');
    plot(eventTimes, zeros(size(eventTimes)), 'g.');
    legend('stim', 'sr');
    title('Linear Currents');
    
    subplot(3, 1, 3); hold on;
    plot(tpl, rate, 'k-');
    plot(eventTimes, zeros(size(eventTimes)), 'r.');
    title('Rate/Event Times');
    
    %generate a bunch of spike trains and see if they
    %match up to the rate
    numTrials = 100;
    psth = zeros(1, round(duration/binSize)+1);
    for k = 1:numTrials      
      [etimes, rate, stimCurrent, srCurrent] = simPP2(pp);    
      eindx = round(etimes / binSize) + 1;
      psth(eindx) = psth(eindx) + 1;             
    end
      
    psth = psth / numTrials;    
    rate = rate / max(rate);
    
    pdiff = norm(psth - rate)
    
    figure; hold on;   
    plot(rate, 'k-', 'LineWidth', 2);
    plot(psth, 'r-');
    axis tight;
    title(sprintf('diff=%f', pdiff));
    