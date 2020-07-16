function signal_out = cosramp(signal_in,ramp_nsamples)

% Apply half-raised-cosine ramp to a signal

ramp_nsamples = round(ramp_nsamples);

cos_mask = hanning(2*ramp_nsamples-1);
%mask = ones(1,length(signal_in));
mask = ones(size(signal_in));
mask(1) = 0;
mask(2:ramp_nsamples) = cos_mask(1:ramp_nsamples-1);
mask(end) = 0;
mask(end-ramp_nsamples+1:end-1) = cos_mask(ramp_nsamples+1:end);

signal_out = signal_in .* mask;