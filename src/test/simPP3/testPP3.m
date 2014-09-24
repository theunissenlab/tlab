function fitStrf = testPP3(strf)

    duration = 1000;
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
    sr = [3 2 1];    
    pp = createPP3(strf, sr, duration, binSize, stim);
    
    fprintf('Generating event times...\n');
    [eventTimes, rate, stimCurrent, srCurrent] = simPP3(pp);
          
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
    initialGuess = randn(1, length(sr))*1e-3;    
    %initialGuess = sr;
    fitParams = fitPP3(pp, eventTimes, initialGuess);
    
    actualParams = pp.sr
    fitParams
    pdiff = norm(actualParams - fitParams)
    
    ppFit = pp;
    ppFit.sr = fitParams;
    
    %generate psths from each model to compare     
    if makePlots
      [eventTimesReal, realRate, stimCurrentReal, srCurrentReal] = simPP3(pp);
      [eventTimesFit, fitRate, stimCurrentFit, srCurrentFit] = simPP3(ppFit);
    
      realRate = realRate / max(realRate);
      fitRate = fitRate / max(fitRate);
    
      figure; hold on;
      plot(realRate, 'k-', 'LineWidth', 2);
      plot(fitRate, 'r-');
      legend('Real', 'Fit');
      axis tight;
    end
    
    fitStrf = fitParams;
    