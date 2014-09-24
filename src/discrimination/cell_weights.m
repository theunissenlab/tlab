function wts = cell_weights(spike_times,stim_lens,tau)

% Get various dimensions
nstims = numel(spike_times);
nreps = numel(spike_times{1});
if iscell(spike_times{1}{1})
    ncells = numel(spike_times{1}{1});
else
    ncells = 1;
end
if ncells == 1
    wts = 1;
    return;
end
tempdists = zeros(ncells,nstims,nstims);
wts = zeros(1,ncells);

% Get template length
template_npts = 1000;   % This matches the template length in VRdistance
if (min(stim_lens) < template_npts) 
        fprintf(1, 'Warning in cell_weights: the template length (%d) is longer than the smallest stimulus (%d ms)\n', ...
        template_npts, min(stim_lens));
end

% The exponential decay used for the templates
expwav = 0:5*round(tau);
expwav = exp(-expwav./tau);

% An index for the upper diagonals
utri_inds = find(triu(ones(nstims),1));

for ci = 1:ncells
    disp(sprintf('Checking template distances for cell %d / %d...',ci,ncells));
    %nspikes = sum(cellfun(@(y) sum(cellfun(@(x) sum(numel(x{ci})), y)), spike_times));
    for temp1 = 1:(nstims-1)
        for temp2 = 2:nstims
            %% Make waveform for first template, first cell
            spikes1 = zeros(1, template_npts);
            for ri = 1:nreps
                good_spikes_ind = find(round(spike_times{temp1}{ri}{ci}) > 0 & round(spike_times{temp1}{ri}{ci}) <= template_npts);
                spikes1(round(spike_times{temp1}{ri}{ci}(good_spikes_ind))) = spikes1(round(spike_times{temp1}{ri}{ci}(good_spikes_ind))) + 1;
            end
            
            %spikes1_wav = filter(expwav,1,spikes1')'; %% replaced with faster fftfilt (0.11 vs. 0.16 seconds / call)
            spikes1_wav = fftfilt(expwav./nreps,spikes1')';
            
            %% Make waveform for second spike train
            spikes2 = zeros(1, template_npts);
            for ri = 1:nreps
                good_spikes_ind = find(round(spike_times{temp2}{ri}{ci}) > 0 & round(spike_times{temp2}{ri}{ci}) <= template_npts);
                spikes2(round(spike_times{temp2}{ri}{ci}(good_spikes_ind))) = spikes2(round(spike_times{temp2}{ri}{ci}(good_spikes_ind))) + 1;
            end
            
            %spikes2_wav = filter(expwav,1,spikes2')'; %% replaced with faster fftfilt (0.11 vs. 0.16 seconds / call)
            spikes2_wav = fftfilt(expwav./nreps,spikes2')';
            
            %% Store distance
            tdiff = spikes1_wav-spikes2_wav;
            tempdist = norm(tdiff(:));
            %tempdist = norm(spikes1_wav/norm(spikes1_wav) - spikes2_wav/norm(spikes1_wav));
            tempdists(ci,temp1,temp2) = tempdist;
            tempdists(ci,temp2,temp1) = tempdist;
            
        end
    end
    tdcell = squeeze(tempdists(ci,:,:)); 
    wts(ci) = mean(tdcell(utri_inds));   % Averaging the upper triangular part of the distance matrix
end
