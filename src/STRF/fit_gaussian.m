function [A, X0, Bw, fitvals, sse] = fit_gaussian(vals);
% fits a 1-D array to a Gabor Function. 

% find maximum time and frequency
maxvals = max(vals);
minvals = min(vals);

if ( abs(maxvals) > abs(minvals) )
    X0_start = find(vals==maxvals);
    A_start = maxvals;
    Xf = find(vals==minvals);
else
    X0_start = find(vals==minvals);
    A_start = minvals;
    Xf = find(vals==maxvals);
end;

xvals=1:length(vals);

% obtain a guess for bandwidth and frequency.    
Bw_start  = abs(X0_start-Xf)/sqrt(2*abs(log(abs(maxvals)/abs(minvals))));

params(1)= X0_start;
params(2)= A_start;
params(3) = Bw_start;


[sse fitvals] = gaussianfun(params);
% fprintf(1, 'Error before fit %f\n', sse);
options = optimset('fminsearch');
options = optimset(options,  'MaxFunEvals', 10000, 'TolX', 1.0000e-003, 'Display', 'off' );
params = fminsearch(@gaussianfun, params);
X0 = params(1);
A = params(2);
Bw = params(3);

[sse fitvals] = gaussianfun(params);
% fprintf(1, 'Error after fit %f\n', sse);

 function [sse, FittedCurve] = gaussianfun(params)
        X0 = params(1);
        A = params(2);
        Bw = params(3);

        if (X0 < xvals(1) )
            X0 = xvals(1);
        elseif (X0 > xvals(end) )
            X0 = xvals(end);
        end
        if Bw > 3*length(xvals)
            Bw = 3*length(xvals);
        end

        
        FittedCurve = A .* exp(-(0.5).*((xvals-X0)./Bw).^2);
      
        ErrorVector = FittedCurve - vals';
        sse = sum(ErrorVector .^ 2);
 end
end