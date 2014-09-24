function d_VR = VR_multidistance(spike_times1, spike_times2, stim_len1, stim_len2, winSize, num_trains, cell_wts)
% Calculates the multineuron Von Rossum distance between a spike train and a template
% constructed from a collection of #num_trains spike trains
% as defined by Houghton and Sen (2008)

min_stim_len = ceil(min(stim_len1, stim_len2)); % This line to use a full window
template_npts = min_stim_len;  % Set at 1000 ms for songs
template_npts = 1000;

% theta is the 'angle' between cells.  theta = pi/2 means that each cell is independent, and that a spike from one cell
% should have no bearing on any other cell.  theta = 0 means that we don't care which cell a spike comes frome, so we 
% treat each cell as equivalent, effectively putting all the spikes in one 'big' cell.  Intermediate thetas yield inter-
% mediate results.
% theta = pi/2; 
theta = pi/2;

ncells = length(spike_times1);

thetas = zeros(1,ncells);
vectors = zeros(ncells);
vectors(1,1) = 1;
thetas(1) = theta;
for vi = 2:ncells
    rmat = eye(ncells);
    th = thetas(vi-1);
    rmat(vi-1:vi,vi-1:vi) = [cos(th), -sin(th); sin(th), cos(th)];
    %vectors(vi-1:vi,vi) = vectors(vi-1:vi,vi-1) .* [cos(thetas(vi-1)); sin(thetas(vi-1))];
    vectors(:,vi) = rmat * vectors(:,vi-1);
    thetas(vi) = acos( (cos(thetas(vi-1)) - cos(thetas(vi-1)).^2) / sin(thetas(vi-1)).^2 );
end

debugfigs = 0;
tauVR = winSize; % ms for the exponential smoothing

if (template_npts > min_stim_len)
    fprintf(1, 'Warning in indeal_VR: the template length (%d) is longer than the smallest stimulus (%d ms)\n', ...
        template_npts, min_stim_len);
end

expwav = 0:5*round(tauVR);
expwav = exp(-expwav./tauVR);
exp_npts = length(expwav);


% Make waveform for first spike train
spikes1 = zeros(ncells, template_npts);
for ci = 1:ncells
    good_spikes_ind = find(round(spike_times1{ci}) > 0 & round(spike_times1{ci}) <= template_npts);
    spikes1(ci,round(spike_times1{ci}(good_spikes_ind))) = spikes1(ci,round(spike_times1{ci}(good_spikes_ind))) + 1;
end

%spikes1_wav = filter(expwav,1,spikes1')'; %% replaced with faster fftfilt (0.11 vs. 0.16 seconds / call)
spikes1_wav = fftfilt(expwav,spikes1')';


% Make waveform for second spike train
nreps2 = size(spike_times2,1);
spikes2 = zeros(ncells, template_npts);
for ri = 1:nreps2
    for ci = 1:ncells
        good_spikes_ind = find(round(spike_times2{ri,ci}) > 0 & round(spike_times2{ri,ci}) <= template_npts);
        spikes2(ci,round(spike_times2{ri,ci}(good_spikes_ind))) = spikes2(ci,round(spike_times2{ri,ci}(good_spikes_ind))) + 1;
        
    end
end

%spikes2_wav = filter(expwav,1,spikes2')'; %% replaced with faster fftfilt (0.11 vs. 0.16 seconds / call)
spikes2_wav = fftfilt(expwav./num_trains,spikes2')';


% Rotate the filtered spike rasters through theta
spikes1_rot = vectors * spikes1_wav;
spikes2_rot = vectors * spikes2_wav;

% Calculate the euclidian distances with shifts -- there are a few methods to do this, to switch between methods uncomment 
%d_VR = norm(spikes1_rot(:)/norm(spikes1_rot(:)) - spikes2_rot(:)/norm(spikes2_rot(:))); %% normalized L2-norm
%d_VR = norm(spikes1_rot(:) - spikes2_rot(:)); %% not  - normalized L2-norm


%% Calculate distance per cell, then weight distances and add them together
cell_dists = zeros(1,ncells);
for ci = 1:ncells
    cell_dists(ci) = norm(spikes1_rot(ci,:) - spikes2_rot(ci,:));
    %cell_dists(ci) = norm(spikes1_rot(ci,:)/(norm(spikes1_rot(ci,:))+eps) - spikes2_rot(ci,:)/(norm(spikes2_rot(ci,:))+eps));
end

weighted_cell_dists = cell_dists .* cell_wts;
d_VR = sum(weighted_cell_dists);


%d_VR = sum(sum(abs(spikes1_rot/max(spikes1_rot(:)) - spikes2_rot/max(spikes2_rot(:))))); %% L1-norm

%cc = corrcoef(spikes1_rot(:),spikes2_rot(:)); %% Correlation coeffiicent
%d_VR = 1 - cc(1,2);

%d_VR = -1 * spikes1_rot(:)' * spikes2_rot(:); %% Dot product

return