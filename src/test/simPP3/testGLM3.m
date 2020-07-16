function testGLM3()

    duration = 100;
    binSize = 0.001;

    %generate stimulus
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen);
	
    %generate fitting structure           
    strf = [0.75 0.4]*2;
    sr = [-0.5 -0.4 -0.3 -0.2 -0.1 0];    
    pp = createPP3(strf, sr, duration, binSize, stim);
    
    %generate event times
    [eventTimes, rate, stimCurrent, srCurrent] = simPP3(pp);
    
    x = [sr];
    
    %test gradent of likelihood function
    [err, grad] = nlLikelihood3(x, pp, eventTimes);

    deps = 1e-10;
    gradCd = zeros(size(grad));
    gradFd = zeros(size(grad));
    for m = 1:length(x)
        xFwd = x;
        xBack = x;	
        xFwd(m) = x(m) + deps;
        xBack(m) = x(m) - deps;

        errFwd = nlLikelihood3(xFwd, pp, eventTimes);
        errBack = nlLikelihood3(xBack, pp, eventTimes);

        gradCd(m) = (errFwd - errBack) / (2*deps);
	gradFd(m) = (errFwd - err)  / deps;
    end
    
    
    grad
    gradCd
    gradFd
    gdiff = norm(grad-gradFd);