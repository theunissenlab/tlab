function fitParams = fitPP2(pp, eventTimes, initialGuess)

    objFunc = @(x) nlLikelihood(x, pp, eventTimes);
    options = optimset('GradObj', 'on', 'DerivativeCheck', 'on', 'Display', 'iter');
    fitParams = fminunc(objFunc, initialGuess, options);    
