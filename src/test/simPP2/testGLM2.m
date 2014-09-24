function testGLM2()

    duration = 100;
    binSize = 0.001;

    %generate stimulus
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen);
	
    %generate fitting structure           
    strf = [0.75 0.4]*2;
    sr = [-0.5 -0.2];    
    pp = createPP(strf, sr, duration, binSize, stim);
    
    %generate event times
    [eventTimes, rate, stimCurrent, srCurrent] = simPP2(pp);
    
    x = [strf, sr];
    
    %test gradent of likelihood function
    [err, grad] = nlLikelihood(x, pp, eventTimes);

    deps = 1e-10;
    gradCd = zeros(size(grad));
    gradFd = zeros(size(grad));
    for m = 1:length(x)
        xFwd = x;
        xBack = x;	
        xFwd(m) = x(m) + deps;
        xBack(m) = x(m) - deps;

        errFwd = nlLikelihood(xFwd, pp, eventTimes);
        errBack = nlLikelihood(xBack, pp, eventTimes);

        gradCd(m) = (errFwd - errBack) / (2*deps);
	gradFd(m) = (errFwd - err)  / deps;
    end
    
    
    grad
    gradCd
    gradFd
    gdiff = norm(grad-gradFd);