function [eventTimes, rate, stimCurrent, srCurrent] = simPP2(pp)

    duration = pp.duration;
    binSize = pp.binSize;
    
    t = 0:binSize:duration;
    tindx = round(t * pp.sampleRate) + 1;
    
    eventTimes = [];    
    rate = zeros(size(t));    
    srCurrent = zeros(size(t));
    
    %generate stimulus current through convolution
    stimCurrent = conv(pp.strf, pp.stim);
    stimCurrent = stimCurrent(1:length(tindx));    
    linCurrent = stimCurrent;
    
    rateSum = 0;
    nextEventTime = exprnd(1);
    
    %generate event times using time-rescaling
    maxT = max(tindx);
    for k = 1:length(tindx)
        
        tk = tindx(k);
        
        %evalute rate at time tk
        rval = pp.nl(linCurrent(tk));
        rateSum = rateSum + rval*binSize;
        rate(tk) = rval;
        
        if rateSum >= nextEventTime
            %record event time, reset rate sum
            etime = (tk-1) / pp.sampleRate;
            eventTimes = [eventTimes etime];
            rateSum = 0;
            nextEventTime = exprnd(1);
            
            %add spike current to linear current
            if ~isempty(pp.sr)
                swts = pp.sr;
                stend = min(tk+length(swts), maxT);
                stindx = (tk+1):stend;
                linCurrent(stindx) = linCurrent(stindx) + swts(1:length(stindx));
                srCurrent(stindx)  = srCurrent(stindx) + swts(1:length(stindx));
            end
        end 
        
    end
    