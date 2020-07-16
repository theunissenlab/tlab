function [d_SR, template_nspikes, test_nspikes] = SR_distanceB(spike_times1, spike_times2, stim_len1, stim_len2, winSize, num_trains, weights, templateSize, norm_flg, shift_ST, gauss_flg)
% Calculates the distance between a spike train and a template
% constructed from a collection of #num_trains spike trains
% The distance is calculated as a difference of spike count
% the winSize, weights, gauss_flg, norm_flg and shift_ST are not used here but are needed for the ensemble version of
% info_distanceB or for the Van Rossum version of info_distanceB

if nargin <= 7 
    template_npts = 1000;    % Use the first 1.0 seconds as a template
else
    template_npts = fix(templateSize);
end

min_stim_len = ceil(min(stim_len1, stim_len2)); % This line to use a full window


if (template_npts > min_stim_len)
    fprintf(1, 'Warning in SR_distanceB: the template length (%d) is longer than the smallest stimulus (%d ms)\n', ...
        template_npts, min_stim_len);
end


% Count Number of spikes for first spike train
time_ind1 = round(spike_times1);
template_nspikes = sum((time_ind1 > 0).*(time_ind1 <= template_npts));

% Count number of spikes for second spike train
time_ind2 = round(spike_times2);
test_nspikes = sum((time_ind2 > 0).*(time_ind2 <= template_npts));

%Divide spike count by the number of possible spikes (number of
%trains included)
test_nspikes = test_nspikes./num_trains;


% Calculate the difference between spike counts
d_SR = abs(template_nspikes - test_nspikes);


return