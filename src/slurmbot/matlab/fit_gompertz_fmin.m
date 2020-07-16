function [b, c] = fit_gompertz_fmin(xvals, yvals, a, b0, c0)

    lb = [-Inf -Inf];
    ub = [-1e-3 -1e-3];
    
    if nargin < 4
        b0 = -abs(randn());
        c0 = -abs(randn());
    end
    
    x0 = [b0 c0];
    xvals = rv(xvals);
    yvals = rv(yvals);

    %opts = optimset('DerivativeCheck', 'on', 'Jacobian', 'on', 'Diagnostics', 'on');    
    %objFunc = @(x, xvals) gompertz_err_lsq(x, xvals, yvals, a);    
    %x = lsqcurvefit(objFunc, x0, xvals, yvals, lb, ub, opts);
        
    opts = optimset('DerivativeCheck', 'off', 'GradObj', 'on', 'Diagnostics', 'off');
    objFunc = @(x) gompertz_err_fmin(x, xvals, yvals, a);  
    x = fmincon(objFunc, x0, [], [], [], [], lb, ub, [], opts);
    
    b = x(1);
    c = x(2);
    
end


function [err, g] = gompertz_err_lsq(x, xvals, yvals, a)
    
    x = rv(x);
    xvals = rv(xvals);
    yvals = rv(yvals);
    
    b = x(1);
    c = x(2);
    y = gompertz(a, b, c, xvals);
    y = rv(y);
    
    rdiff = rv(y - yvals);
    err = rdiff.^2;
        
    if nargout > 1
        g = zeros(length(xvals), 2);
        
        gc = b*exp(c*xvals);
        %g(:, 1) = 2*rdiff .* exp(gc);
        g(:, 1) = 2*rdiff .* (a*exp(c*xvals + gc));
        g(:, 2) = 2*rdiff .* (a*b*xvals) .* (exp(c*xvals + gc));
    end
    
end

function [err, g] = gompertz_err_fmin(x, xvals, yvals, a)
    
    x = rv(x);
    xvals = rv(xvals);
    yvals = rv(yvals);
    
    b = x(1);
    c = x(2);
    y = gompertz(a, b, c, xvals);
    y = rv(y);
    
    rdiff = rv(y - yvals);
    err = sum(rdiff.^2);
        
    if nargout > 1
        g = zeros(1, 2);
        
        gc = b*exp(c*xvals);
        %g(1) = 2*rdiff .* exp(gc);
        g(1) = sum(2*rdiff .* (a*exp(c*xvals + gc)));
        g(2) = sum(2*rdiff .* (a*b*xvals) .* (exp(c*xvals + gc)));
    end
    
end
