function fitParams = fitPP(pp, eventTimes, initialGuess)

    objFunc = @(x) nlLikelihood(x, pp, eventTimes);
    options = optimset('GradObj', 'on', 'DerivativeCheck', 'off', 'Display', 'iter');
    fitParams = fminunc(objFunc, initialGuess, options);    
    
end

function [err, grad] = nlLikelihood(params, pp, eventTimes)
    
    t = 0:pp.binSize:pp.duration;
    
    [r, dr] = pp.rateFunc(t, params);
    [re, dre] = pp.rateFunc(eventTimes, params);

    %{
    rsz = size(r)
    drsz = size(dr)
    resz = size(re)
    dresz = size(dre)
    %}
    
    zindx = (re == 0);
    re(zindx) = 1e-100;
    dre(zindx) = 0;
        
    l1 = sum(log(re));
    l2 = sum(r)*pp.binSize;
    
    err = -(l1 - l2);

    if nargout > 1       
        
        reInv = re .^ -1;
        
        gl1Mat = zeros(size(dre));
        for k = 1:length(reInv)
            gl1Mat(:, k) = reInv(k) * dre(:, k);
        end        
        gl1 = sum(gl1Mat, 2);        
        gl2 = sum(dr, 2)*pp.binSize;
         
        grad = -(gl1 - gl2);
        grad;
    end
    
end
