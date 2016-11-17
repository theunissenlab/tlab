function [d_VR, template_wav, test_wav] = VR_distanceB(spike_times1, spike_times2, stim_len1, stim_len2, winSize, num_trains, weights, templateSize, norm_flg, shift_ST, gauss_flg)
% Calculates the Von Rossum distance between a spike train and a template
% constructed from a collection of #num_trains spike trains
% the weights are not used here but are needed for the ensemble version of
% info_distanceB
if nargin<= 10
    gauss_flg=0;
end

if nargin<=9
    shift_ST=0; % Set to 1 to find the best allignment of spike train (the one that give the shorter distance)
end

if nargin <= 8
    norm_flg=0; % Set to 1 to nomarlize distance by rate
end

if nargin <= 7 
    template_npts = 1000;    % Use the first 1.0 seconds as a template
else
    template_npts = fix(templateSize);
end

min_stim_len = ceil(min(stim_len1, stim_len2)); % This line to use a full window


debugfigs = 0;
tauVR = round(winSize);         % ms for the exponential or gaussian smoothing
%tauVR2 = winSize/num_trains; % the same exponential is used for both
dVR=round(tauVR/2);

if (template_npts > min_stim_len)
    fprintf(1, 'Warning in indeal_VR: the template length (%d) is longer than the smallest stimulus (%d ms)\n', ...
        template_npts, min_stim_len);
end


if gauss_flg
    nStd =3;
    t_pts = 0:2*nStd*tauVR;
    expwav = exp(-0.5*(t_pts-nStd*tauVR).^2./tauVR^2);
    expwav2 = expwav;
else
    t_pts=0:4*tauVR;
    expwav = exp(-t_pts./tauVR);
    % expwav2 = exp(-t_pts./tauVR2);   % Use the same tau for both...
    expwav2 = expwav;
end
expwav = expwav./sum(expwav);
expwav2 = expwav2./sum(expwav2);

% Make waveform for first spike train
time_ind = round(spike_times1);
nspikes = length(time_ind);
template_spikes = zeros(1, template_npts);
for is=1:nspikes
    if (time_ind(is) > 0 && time_ind(is) <= template_npts)
        template_spikes(1,time_ind(is)) = template_spikes(1,time_ind(is)) + 1;
    end
end

template_wav = conv(template_spikes, expwav, 'same');
template_npts2=length(template_wav);


% Make waveform for second spike train
time_ind = round(spike_times2);
nspikes = length(time_ind);
test_spikes = zeros(1, template_npts);
for is=1:nspikes
    if (time_ind(is) > 0 && time_ind(is) <= template_npts)
        test_spikes(1,(time_ind(is))) = test_spikes(1,time_ind(is)) + 1;
    end
end




%Divide each spike window by the number of possible spikes (number of
%trains included)
test_spikes = test_spikes./num_trains;

test_wav = conv(test_spikes, expwav2, 'same');
test_npts = length(test_wav);
meanVal = mean(test_wav);

if shift_ST==1
    %extend the second spike train to allow for maximum jitter between spike trains and
    % returns optimal value.
    test_wav=[zeros(1,template_npts2) test_wav zeros(1,template_npts2)];
    test_wav(1,1:template_npts2)=meanVal;%the extended spike wa is padded with a model of the spike train=the mean rate per ms
    test_wav(1,(test_npts+template_npts2+1):end)=meanVal;%the extended spike wav is padded with a model of the spike train=the mean rate per ms
end





% Calculate the euclidian distances with shifts

if shift_ST==0 %no shift, a single value is calculated with alligned waveform
    test_wav_local=test_wav;
    if (norm_flg == 1)  % normalize distance
        norm_test = norm(test_wav_local);
        norm_template = norm(template_wav);
        if (norm_test == 0 )
            norm_test = 1;
        end
        if (norm_template == 0)
            norm_template = 1;
        end        
        dvals = norm(test_wav_local./norm_test - template_wav./norm_template);   
    else  % don't normalize
        dvals = norm(test_wav_local - template_wav);
    end
else
    Step = 1:dVR:(test_npts+template_npts2);
    dvals = zeros(1, length(Step));
    for tt=1:length(Step);
        dt=Step(tt);
        %select the area of the spike trains on which the distance is calculated
        test_wav_local=test_wav(dt:(dt+template_npts2-1));
        template_wav_local=template_wav;
        if (norm_flg == 1)  % normalize distance
            norm_test = norm(test_wav_local);
            norm_template = norm(template_wav_local);
            if (norm_test == 0 )
                norm_test = 1;
            end
            if (norm_template == 0)
                norm_template = 1;
            end        
            dvals(1,tt) = (norm(test_wav_local./norm_test - template_wav_local./norm_template)); 
        else  % don't normalize
            dvals(1,tt) = (norm(test_wav_local - template_wav_local)); 
        end
    end
end
d_VR = min(dvals);
if length(d_VR)>1
    fprintf(1,'Warning, there is more than one minimum distance here\n');
end
offset_to_template = find(dvals == d_VR, 1);
d_VR = sqrt(d_VR);                              % The square root makes the distribution more Normal


% Debugging figure
if (debugfigs)
    figure(2)
    figure_wav = zeros(2, (length(test_wav)));
    figure_wav(1,1:length(test_wav)) = test_wav;
    figure_wav(2,offset_to_template:offset_to_template+template_npts2-1) = template_wav;
    %maxX=max(length(test_wav), (-dt+length(template_wav)-1));
    %maxY=max(max(test_wav),max(template_wav));
    plot(figure_wav(1,:));
    hold on;
    plot(figure_wav(2,:),'r');
    vline(offset_to_template)
    vline(offset_to_template + template_npts2-1)
    vline(template_npts2, 'b:')
    vline((length(test_wav)-template_npts2), 'b:')
    %axis([0 maxX 0 maxY+0.1])
    title(sprintf('%f', d_VR));
    hold off;
    pause;
    
end


return