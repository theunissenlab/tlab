function display_model_data(modelData, displayAllStims)
 
    if ~exist('displayAllStims', 'var')
        displayAllStims = 0;
    end

    m = modelData;
    grps = unique(m.data.groupIndex);    
    
    if displayAllStims
        %% display stimuli and responses
        for k = 1:length(grps)        
            indx = find(m.data.groupIndex == grps(k));

            t = ((1:length(indx)) - 1) / m.data.sampleRate;
            f = m.data.f;
            stim = m.data.wholeStim(indx, :)';
            realResp = m.data.wholeResp(indx);
            modelResp = m.model.response(indx);

            figure; hold on;

            %subplot(2, 1, 1); hold on;
            imagesc(t, f, stim); axis tight;
            title(sprintf('Stim %d', grps(k)));

            %{
            subplot(2, 1, 2); hold on;
            plot(t, realResp, 'k-', 'LineWidth', 2);
            plot(t, modelResp, 'r-');               
            axis tight;
            %}
        end
    end

    figure; hold on;
    
    %% plot STRF
    subplot(2, 2, 1); hold on;    
    imagesc(1:size(m.model.strf, 2), m.data.f, m.model.strf);
    cval = max(abs(m.model.strf(:)));
    caxis([-cval cval]);
    title(sprintf('STRF: bias=%0.4f', m.model.bias));
    axis tight; colorbar;    
    
    %% plot output NL
    subplot(2, 2, 3); hold on;    
    xmin = min(m.model.response);
    xmax = max(m.model.response);
    xinc = (xmax - xmin) / 200;
    x = xmin:xinc:xmax;
    y = m.model.outputNL(x);
    plot(x, y, 'k-'); axis tight;
    title('Output NL');
    
    %% plot real and model responses
    subplot(2, 2, 2); hold on;
    plot(m.data.wholeResp, 'k-', 'LineWidth', 2);
    plot(m.model.response, 'r-');
    axis tight;
    legend('PSTH', 'Model');
    
    %% plot coherences
    subplot(2, 2, 4); hold on;
    plot(m.coherence.bound.f, m.coherence.bound.mean, 'k-', 'LineWidth', 2);
    plot(m.coherence.validation.f, m.coherence.validation.mean, 'r-', 'LineWidth', 2);
    legend('Bound', 'Model');
    axis tight;
    
    %% model performance
    iBound = m.coherence.bound.info_mean;
    iTrain = m.coherence.training.info_mean;
    iValid = m.coherence.validation.info_mean;
    fprintf('Performance: %s, %s, %s\n', m.cellName, m.data.stimClass, m.responseFile);    
    fprintf('\tInfo Upper Bound: %0.2f bits\n', iBound);
    fprintf('\tTraining Info: %0.2f bits (%0.3f)\n', iTrain, iTrain/iBound);
    fprintf('\tValidation Info: %0.2f bits (%0.3f)\n', iValid, iValid/iBound);
    