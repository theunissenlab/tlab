function y = expDecay(b, x, ceilval)
    % Returns an exponential decay fucntion for fitting purposes
    % expParam sets reasonable limites
    
    [k, t0, tau] = expParam(b);    
    y = k*ceilval*(1.0-exp(-(x-t0)/tau));
    

end