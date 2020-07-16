function fitParams = fitPP3(pp, eventTimes, initialGuess)

    objFunc = @(x) nlLikelihood3(x, pp, eventTimes);
    options = optimset('GradObj', 'on', 'DerivativeCheck', 'on', 'Display', 'iter');
    fitParams = fminunc(objFunc, initialGuess, options);    
