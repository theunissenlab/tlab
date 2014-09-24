function ev=calcExplainableVariance2(d)
% function ev=calcExplainableVariance2(d)
%
% Calculate explainable variance for repeated data
%
% d: a 3D matrix of size [nSamples x nRepeats x nUnits]
%
% (c.f. Pasupathy and Miller 2005; Mante et al. 2008)
%

repnum = size(d,2);
dm = nanmean(d,2);
derr = bsxfun(@minus, d, dm);

derr = reshape(derr, [], size(derr,3));
d = reshape(d, [], size(d,3));

ev = 1-nanvar(derr)./nanvar(d);

% correction for small repeats
ev = ev-(1-ev)/(repnum-1);
