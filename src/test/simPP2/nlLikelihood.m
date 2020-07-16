function [err, grad] = nlLikelihood(params, pp, eventTimes)
    
    xlen = length(pp.strf);
    srlen = length(pp.sr);
    x = params(1:xlen);
    swts = params(xlen+1:end);

    t = 0:pp.binSize:pp.duration;
    tindx = round(t * pp.sampleRate) + 1;
    eindx = round(eventTimes * pp.sampleRate) + 1;
    maxT = max(tindx);
    
    %% compute error
    
    %generate stimulus current through convolution
    stimCurrent = conv(x, pp.stim);
    stimCurrent = stimCurrent(1:length(tindx));  
    
    %generate spike current through convolution
    srCurrent = zeros(size(tindx));
    if ~isempty(swts)
      for k = 1:length(eventTimes)            
        ek = round(eventTimes(k) * pp.sampleRate) + 1;        
        eend = min(ek+length(swts), maxT);
        etindx = (ek+1):eend;
        srCurrent(etindx) = srCurrent(etindx) + swts(1:length(etindx));
      end
    end
    
    %generate linear current
    linCurrent = stimCurrent + srCurrent;
    
    %generate nonlinear rate
    rate = pp.nl(linCurrent);
    
    %compute terms of likelihood function
    l1 = sum(log(rate(eindx)));
    l2 = sum(rate)*pp.binSize;
    
    %compute error, negative of likelihood
    err = -(l1 - l2);

    %% compute gradient
    if nargout > 1       

        drate = pp.dnl(linCurrent);        
        idrate = drate .^ -1;        
	
        %form matrix of stimulus history vectors, and indicator
    	%vectors for gradient w/ respect to sr filter
        shist = zeros(length(x), length(tindx));
    	indmat = zeros(length(swts), length(tindx));
        for k = 1:length(tindx)
            tk = tindx(k);       
            sstart = max(1, tk-length(x)+1);
            send = tk;
            sindx = sstart:send;      
            shist(1:length(sindx), tk) = pp.stim(sindx);

            tdiff = tk - eindx(eindx < tk & eindx >= (tk-length(swts)));
            for m = 1:length(swts)
              indmat(m, tk) = sum(tdiff == m);
            end	    
        end
        shist = flipud(shist);
        fullMat = [shist; indmat];
        
        %compute both parts of likelihood gradient
        rhs = (idrate .* drate)';
        g1 = fullMat(:, eindx) * rhs(eindx);
        g2 = (fullMat * drate')*pp.binSize;
         
        %gradient of error is negative gradient of likelihood
        grad = -(g1 - g2);
	
    end
    

