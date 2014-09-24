function [A, X0, Sig, Bw, P, fitvals, sse] = fit_gabor(vals);
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
P_start = 0;
xvals=1:length(vals);

% obtain a guess for bandwidth and frequency.
if ( maxvals*minvals < 0 )  % Crossing point defines frequency - BW is set at twice that value
    Sig_start = 1/(2*abs(X0_start-Xf));
    Bw_start  = abs(X0_start-Xf)/sqrt(2*abs(log(abs(maxvals)/abs(minvals))));
else
    Bw_start  = abs(X0_start-Xf)/sqrt(2*abs(log(abs(maxvals)/abs(minvals))));
    Sig_start = 1/(4*abs(X0_start-Xf));
end
params(1)= X0_start;
params(2)= A_start;
params(3)= Sig_start;
params(4) = Bw_start;
params(5)= P_start;

[sse fitvals] = gaborfun(params);
% fprintf(1, 'Error before fit %f\n', sse);
options = optimset('fminsearch');
options = optimset(options,  'MaxFunEvals', 10000, 'TolX', 1.0000e-003, 'Display', 'off' );
params = fminsearch(@gaborfun, params);
X0 = params(1);
A = params(2);
Sig = params(3);
Bw = params(4);
P = params(5);
[sse fitvals] = gaborfun(params);
% fprintf(1, 'Error after fit %f\n', sse);

 function [sse, FittedCurve] = gaborfun(params)
        X0 = params(1);
        A = params(2);
        Sig = params(3);
        Bw = params(4);
        P = params(5);
        if (X0 < xvals(1) )
            X0 = xvals(1);
        elseif (X0 > xvals(end) )
            X0 = xvals(end);
        end
        if Bw > 2*length(xvals)
            Bw = 2*length(xvals);
        end
        if Sig < 1/(3*length(xvals))
            Sig = 1/(3*length(xvals));
        end
        
        FittedCurve = A .* exp(-(0.5).*((xvals-X0)./Bw).^2).*cos(2*pi()*Sig*(xvals-X0)+P);
      
        ErrorVector = FittedCurve - vals';
        sse = sum(ErrorVector .^ 2);
 end
end