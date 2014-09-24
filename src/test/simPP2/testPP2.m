function fitStrf = testPP2(strf)

    duration = 500;
    binSize = 0.001;
    tpl = 0:binSize:duration;

    %generate stimulus
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen) + 1;
	
    %generate fitting structure           
    if nargin < 1
      strf = [0.75 0.4]*0;
      makePlots = 1;
    else
      makePlots = 0;
    end
    sr = [3 2 1 0];    
    pp = createPP(strf, sr, duration, binSize, stim);
    
    fprintf('Generating event times...\n');
    [eventTimes, rate, stimCurrent, srCurrent] = simPP2(pp);
    fprintf('# of event times: %d\n', length(eventTimes));
          
    if makePlots
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
    end  
  
    %fit the model with the data
    fprintf('Fitting...\n');
    initialGuess = randn(1, length(strf) + length(sr));    
    fitParams = fitPP2(pp, eventTimes, initialGuess);
    
    actualParams = [pp.strf pp.sr]
    fitParams
    pdiff = norm(actualParams - fitParams)
    
    ppFit = pp;
    ppFit.strf = fitParams;
    
    %generate psths from each model to compare     
    if makePlots
      [eventTimesReal, realRate, stimCurrentReal, srCurrentReal] = simPP2(pp);
      [eventTimesFit, fitRate, stimCurrentFit, srCurrentFit] = simPP2(ppFit);
    
      realRate = realRate / max(realRate);
      fitRate = fitRate / max(fitRate);
    
      figure; hold on;
      plot(realRate, 'k-', 'LineWidth', 2);
      plot(fitRate, 'r-');
      legend('Real', 'Fit');
      axis tight;
    end
    
    fitStrf = fitParams;
    
    
    