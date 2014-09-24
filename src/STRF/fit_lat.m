function lat = fit_lat(t, forward);
% Returns the latency of peak.

% find maximum time and frequency
maxforward = max(max(forward));
minforward = min(min(forward));
absforward = max(abs(minforward),abs(maxforward));
if ( abs(maxforward) > abs(minforward) )
    [fpeak tpeak] = find(forward==maxforward);
else
    [fpeak tpeak] = find(forward==minforward);
end;

lat = t(tpeak);