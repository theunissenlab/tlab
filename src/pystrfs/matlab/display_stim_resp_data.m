function display_stim_resp_data(stimRespData)

    groupIndex = stimRespData.groupIndex;
    wholeStim = stimRespData.wholeStim;
    wholeResp = stimRespData.wholeResp;

    %% display stim and responses
    nGroups = length(unique(groupIndex));
    for k = 1:nGroups

        tRng = find(groupIndex == k);

        stim = wholeStim(tRng, :);
        resp = wholeResp(tRng);

        tInc = 1 / stimRespData.sampleRate;
        t = 0:tInc:(size(stim, 1)-1)*tInc;

        figure; hold on;

        %plot spectrogram
        subplot(2, 1, 1);
        imagesc(t, stimRespData.f, stim'); axis tight;
        axis xy; colorbar;
        v_axis = axis;
        v_axis(1) = min(t); v_axis(2) = max(t);
        v_axis(3) = min(stimRespData.f); v_axis(4) = max(stimRespData.f);
        axis(v_axis);
        xlabel('Time (s)'); ylabel('Frequency (Hz)');

        %plot PSTH
        subplot(2, 1, 2);
        plot(t, resp, 'k-'); 
        xlabel('Time (s)'); ylabel('P[spike]');
        axis([0 max(t) 0 1]);

        title(sprintf('Pair %d', k));
    end
