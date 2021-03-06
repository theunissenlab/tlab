function [r, dr] = glmRate2(t, params, pp)

    x = params;
    binSize = pp.binSize;
    tindx = round(t / binSize) + 1;
    
    linVal = zeros(1, length(tindx));
    r = zeros(1, length(tindx));
    dr = zeros(length(x), length(r));
    
    for k = 1:length(tindx)
      
      tk = tindx(k);
            
      rstart = max(1, tk-length(x)+1);
      rend = tk;
      indx = rstart:rend;      
      stim = pp.stim(indx);
      
      filt = x(1:length(indx));
      linVal(tk) = stim * filt';       
      
      r(k) = exp(linVal);
      
      rhist = zeros(length(x), 1);
      %rhist(1:length(indx)) = fliplr(stim);
      rhist(1:length(indx)) = stim;
      dr(:, k) = r(k) * rhist;
      
    end
    
    
        
