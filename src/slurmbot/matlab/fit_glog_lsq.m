function [B, M] = fit_glog_lsq(xvals, yvals, B0, M0)

    yvals(yvals <= 0) = 1e-6;

    lb = [1e-1 -Inf];
    ub = [Inf Inf];
    
    if nargin < 3
        B0 = abs(randn()*2);
        M0 = randn();
    end
    
    x0 = [B0 M0];
    xvals = rv(xvals);
    yvals = rv(yvals);

    %opts = optimset('DerivativeCheck', 'on', 'Jacobian', 'on', 'Diagnostics', 'on');
    %objFunc = @(x, xvals) glog_err_lsq(x, xvals, yvals);    
    %x = lsqcurvefit(objFunc, x0, xvals, yvals, lb, ub, opts);
    
    opts = optimset('DerivativeCheck', 'on', 'GradObj', 'on', 'Diagnostics', 'on');
    objFunc = @(x) glog_err_fmin(x, xvals, yvals);  
    x = fmincon(objFunc, x0, [], [], [], [], lb, ub, [], opts);
    
    B = x(1);
    M = x(2);
    
end


function [err, g] = glog_err_lsq(x, xvals, yvals)
    
    xvals = rv(xvals);
    yvals = rv(yvals);
    
    B = x(1);
    M = x(2);
    y = glogistic(xvals, B, M);
    y = rv(y);
    
    rdiff = rv(yvals - y);
    err = rdiff.^2;
        
    if nargout > 1
        g = zeros(length(xvals), 2);        
        h = exp(-B*(xvals - M));
        gh = 1 + h;
        
        g(:, 1) = (2*rdiff) .* (gh.^-2) .* h .* (M-xvals);
        g(:, 2) = (2*rdiff) .* (gh.^-2) .* h * B;
    end
    
end

function [err, g] = glog_err_fmin(x, xvals, yvals)
    
    xvals = rv(xvals);
    yvals = rv(yvals);
    
    B = x(1);
    M = x(2);
    y = glogistic(xvals, B, M);
    y = rv(y);
    
    rdiff = rv(yvals - y);
    err = sum(rdiff.^2);
        
    if nargout > 1
        g = zeros(1, 2);        
        h = exp(-B*(xvals - M));
        gh = 1 + h;
        
        g(1) = sum( (2*rdiff) .* (gh.^-2) .* h .* (M-xvals) );
        g(2) = sum( (2*rdiff) .* (gh.^-2) .* h * B );
    end
    
end
