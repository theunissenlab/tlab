function display_stim_resp(preprocData)

    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    stimInfo = preprocData.stimInfo;
    respInfo = preprocData.respInfo;

    %% display stim and responses
    nGroups = length(unique(groupIndex));
    for k = 1:nGroups

        tRng = find(groupIndex == k);

        stim = wholeStim(tRng, :);
        resp = wholeResponse(tRng);

        tInc = 1 / stimInfo.sampleRate;
        t = 0:tInc:(size(stim, 1)-1)*tInc;

        figure; hold on;

        %plot spectrogram
        subplot(2, 1, 1);
        imagesc(t, stimInfo.f, stim'); axis tight;
        axis xy;
        v_axis = axis;
        v_axis(1) = min(t); v_axis(2) = max(t);
        v_axis(3) = min(stimInfo.f); v_axis(4) = max(stimInfo.f);
        axis(v_axis);
        xlabel('Time (s)'); ylabel('Frequency (Hz)');

        %plot PSTH
        subplot(2, 1, 2);
        plot(t, resp, 'k-'); 
        xlabel('Time (s)'); ylabel('P[spike]');
        if ~respInfo.zscored
            axis([0 max(t) 0 1]);
        else
            axis tight;
        end

        title(sprintf('Pair %d', k));
    end
