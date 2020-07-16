function errBndPP()
   
    strf = [0.75 0.4]*2;

    nTrials = 20;
    errs = zeros(nTrials, length(strf));
    for k = 1:nTrials      
      fitStrf = testPP(strf);
      errs(k, :) = strf - fitStrf;          
    end
    
    errs
    avgErr = sum(errs, 1) / nTrials
    squareError = sqrt(sum(errs .^ 2, 1) / nTrials)
    