function [bmf, fpow, pow] = calc_tms(forward, samprate);

maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));
if ( abs(maxforward) > abs(minforward) )
    [fpeak tpeak] = find(forward==maxforward);
else
    [fpeak tpeak] = find(forward==minforward);
end;

n_forward = size(forward,2);
[pow fpow]=periodogram(forward(fpeak,:), ones(1,n_forward),n_forward, samprate);
powmax=max(pow);
bmf_index=find(pow==powmax);
bmf=fpow(bmf_index);