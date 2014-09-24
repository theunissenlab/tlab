function output = resample_lcm(input,fs_in,fs_out,rat_precision)

% Use Matlab's resample function with Patrick Gill's upsample-downsample LCM routine

if nargin < 4
	rat_precision = .0000001;
end

[p,q] = rat(fs_out/fs_in,rat_precision);
output = resample(input,p,q);