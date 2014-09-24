function output = MatchMagnitudes(reference, input)
% MatchMagnitudes to take the magnitude from the reference and the phase from the input
%	output = MatchMagnitudes(reference, input)
% Take the magnitude from the reference and the phase from the input and return a new array.
referenceMag = abs(reference);
inputPhase = angle(input);
output = referenceMag.*exp(i.*inputPhase);
% inputMag = input .* conj(input);
% output = input .* sqrt(referenceMag ./ (inputMag+.00000001));
%ms = s .* conj(s);
%fprintf('reference mag is %g, original mag is %g\n', max(max(m2)), max(max(mi)));
%fprintf('new mag is %g.\n', max(max(ms)));

