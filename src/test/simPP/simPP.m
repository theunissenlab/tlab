function [eventTimes, rate] = simPP(pp)

    duration = pp.duration;
    binSize = pp.binSize;
    
    eventTimes = [];    
    t = 0:binSize:duration;

    rateSum = 0;
    nextEventTime = exprnd(1);
    rate = zeros(size(t));
    
    for k = 1:length(t)
        
        tk = t(k);
        rval = pp.rateFunc(tk, pp.params);
        rateSum = rateSum + rval*binSize;
        rate(k) = rval;
        
        if rateSum >= nextEventTime
            eventTimes = [eventTimes tk];
            rateSum = 0;
            nextEventTime = exprnd(1);
        end 
        
    end
    