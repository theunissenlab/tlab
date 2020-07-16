function d_VR = VR_distance(spike_times1, spike_times2, stim_len1, stim_len2, winSize, templateSize)
% Calculates the Von Rossum distance between two spike trains

if nargin == 5 
    template_npts = 1000;    % Use the first 1.0 seconds as a template
else
    template_npts = fix(templateSize);
end

min_stim_len = ceil(min(stim_len1, stim_len2)); % This line to use a full window

debugfigs = 0;
tauVR = winSize;         % ms for the exponential smoothing

if (template_npts > min_stim_len)
    fprintf(1, 'Warning in indeal_VR: the template length (%d) is longer than the smallest stimulus (%d ms)\n', ...
        template_npts, min_stim_len);
end

expwav = 0:5*round(tauVR);
expwav = exp(-expwav./tauVR);
exp_npts = length(expwav);


% The test_beg variable allows for jitter between spike trains and
% returns optimal -- CURRENTLY NOT IN USE
test_beg = 0;
test_npts = template_npts+1;
diff_npts =  test_npts - template_npts;

% Make waveform for first spike train
time_ind = round(spike_times1);
nspikes = length(time_ind);
template_spikes = zeros(1, template_npts);
for is=1:nspikes
    if (time_ind(is) > 0 && time_ind(is) <= template_npts)
        template_spikes(1,time_ind(is)) = 1;
    end
end
template_wav = conv(template_spikes, expwav);
% template_wav has size template_npts+exp_npts-1

% Make waveform for second spike train
time_ind = round(spike_times2) + test_beg;
nspikes = length(time_ind);
test_spikes = zeros(1, test_npts);
for is=1:nspikes
    if (time_ind(is) > 0 && time_ind(is) <= test_npts)
        test_spikes(1,time_ind(is)) = 1;
    end
end
test_wav = conv(test_spikes, expwav);

% Calculate the euclidian distances with shifts
dvals = zeros(1, diff_npts);
for dt=1:diff_npts
    dvals(1,dt) = sum((test_wav(dt:template_npts+exp_npts+dt-2)-template_wav).^2);
end
d_VR = min(dvals);
offset_to_template = find(dvals == d_VR, 1);
d_VR = sqrt(d_VR);                              % The square root makes the distribution more Normal


% Debugging figure
if (debugfigs)
    figure_wav = zeros(2, test_npts+exp_npts-1);
    figure_wav(1,:) = test_wav;

    figure_wav(2, offset_to_template:offset_to_template+template_npts+exp_npts-2) = template_wav;
    figure(1);
    plot(figure_wav(1,:));
    hold on;
    plot(figure_wav(2,:),'r');
    title(sprintf('%f', d_VR));
    hold off;
    pause;
end


return