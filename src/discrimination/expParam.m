function [k, t0, tau] = expParam(b)
    % Fixes limits on the parameters used in expDecay
    
    k = b(1);
    if k < 0.0
        k = 0.0;
    elseif k > 1.0
        k = 0.99;
    end
    
    t0 = b(2);
    if t0 < 0.0
        t0 = 0.0;
    end
    
    tau = b(3);
    if (tau < 0)
        tau = -tau;
    end
    if tau < 0.1
        tau = 0.1;
    end
    
    
end