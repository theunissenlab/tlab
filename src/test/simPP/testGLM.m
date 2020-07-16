function testGLM()

    duration = 0.050;
    binSize = 0.001;
    
    stimLen = round(duration / binSize) + 1;
    stim = randn(1, stimLen);
    
    pp = struct;
    pp.params = [0.75 0.4]*5;
    pp.duration = duration;
    pp.binSize = binSize;
    pp.stim = stim;    
    pp.rateFunc = @(t, p) glmRate(t, p, pp);
         
    
    x = [0.5, 0.3];
    t = 0:binSize:duration;
    
    eventTimes = simPP(pp);
    
    for k = 1:length(t)
      
      tk = t(k);
      [r, dr] = glmRate(tk, x, pp);
      
      deps = 1e-12;
      drFd = zeros(size(dr));
      for m = 1:length(x)
	xFwd = x;
	xBack = x;	
	xFwd(m) = x(m) + deps;
	xBack(m) = x(m) - deps;
	
	[rFwd, drFwd] = glmRate(tk, xFwd, pp);
	[rBack, drBack] = glmRate(tk, xBack, pp);
	
	drFd(m) = (rFwd - rBack) / (2*deps);	
      end
      
      dr
      drFd
      dddiff = norm(dr - drFd)      
      
    end
    