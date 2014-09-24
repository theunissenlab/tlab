function [cf, bw, forward_f_fit] = fit_cf_bw(f, forward);
% fits the spectral slice with largest amp with a guassian to get bandwidth
% The cf is found at the peak

% find maximum time and frequency
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));
if ( abs(maxforward) > abs(minforward) )
    [fpeak tpeak] = find(forward==maxforward);
else
    [fpeak tpeak] = find(forward==minforward);
end;

% find best guess for standard deviation
f_slice = forward(:,tpeak)';
fw_ind = length(f_slice);
if ( abs(maxforward) > abs(minforward) )
    f_width = exp(-0.5).*maxforward;
    for fp = fpeak+1:length(f_slice)
        if ( f_slice(fp-1) >= f_width & f_slice(fp) < f_width )
            fw_ind = fp;
            break;
        end
    end
else
    f_width = exp(-0.5).*minforward;
    for fp = fpeak+1:length(f_slice)
        if ( f_slice(fp-1) <= f_width & f_slice(fp) > f_width )
            fw_ind = fp;
            break;
        end
    end 
end;


beta(1)= forward(fpeak,tpeak);
beta(2)= f(fpeak);
beta(3)= f(fw_ind)-f(fpeak);

options = optimset('lsqcurvefit');
options = optimset(options,'Jacobian','on','DerivativeCheck','on','MaxFunEvals',6000);
[new_beta,resnorm,residual,exitflg] = lsqcurvefit('gaussfun',beta,f,f_slice,[],[],options);
if (exitflg < 0)
    new_beta = beta;
end
forward_f_fit = gaussfit(new_beta,f);
cf = f(fpeak);
bw = new_beta(3);