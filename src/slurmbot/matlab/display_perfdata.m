function display_perfdata(preprocData, perfData, modelParamsTrained, trainingIndex, validationIndex)

    stimInfo = preprocData.stimInfo;
    respInfo = preprocData.respInfo;
    
    linModel = modelParamsTrained;
    if strcmp('lnl', modelParamsTrained.type)
        linModel = modelParamsTrained.linModel;
    end
    
    %% plot STRF
    figure; hold on;
    imagesc(linModel.delays, stimInfo.f, squeeze(linModel.w1)); axis tight;
    absmax = max(max(abs(linModel.w1)));
    caxis([-absmax absmax]);
    colorbar;
    title(sprintf('STRF | bias=%f', linModel.b1));

    respTrainReal = perfData.wholeResp(trainingIndex);
    respTrain = perfData.train.resp;
    respValidReal = perfData.wholeResp(validationIndex);
    respValid = perfData.valid.resp;
    
    %% plot training response
    tInc = 1 / stimInfo.sampleRate;
    t = 0:tInc:(length(respTrainReal)-1)*tInc;

    figure; hold on;
    plot(t, respTrainReal, 'k-');
    plot(t, respTrain, 'r-');
    if ~respInfo.zscored
        axis([min(t) max(t) 0 1]);
    else
        axis tight;
    end
    title('Training Data');
    legend('Real', 'ThreshGradDesc');

    %% plot early stop response
    tInc = 1 / stimInfo.sampleRate;
    t = 0:tInc:(length(respValidReal)-1)*tInc;
    
    figure; hold on;
    subplot(2, 1, 1); hold on;
    plot(t, respValidReal, 'k-');
    plot(t, respValid, 'r-');
    if ~respInfo.zscored
        axis([min(t) max(t) 0 1]);
    else
        axis tight;
    end
    title('Validation Data');
    
    subplot(2, 1, 2);
    rdiff = rv(respValid) - rv(respValidReal);
    plot(rdiff, 'r-');
    title('Response Differences');    

    legend('Real', 'ThreshGradDesc');
    
    

    %% plot coherences and information for training
    cBoundTrain = perfData.train.cBound;
    cTrain = perfData.train.c;
    cBoundValid = perfData.valid.cBound;
    cValid = perfData.valid.c;

    figure; hold on;
    plot(cBoundTrain.f, cBoundTrain.c, 'k-', 'LineWidth', 2);
    plot(cTrain.f, cTrain.c, 'b-');
    axis tight;
    title(sprintf('Coherences for Training Set, perf ratio=%f', perfData.train.perf));
    legend('Upper Bound', 'ThreshGradDesc');

    %% plot coherences and information for training
    figure; hold on;
    plot(cBoundValid.f, cBoundValid.c, 'k-', 'LineWidth', 2);
    plot(cValid.f, cValid.c, 'b-');
    axis tight;
    title(sprintf('Coherences for Early Stopping Set, perf ratio=%f', perfData.valid.perf));
    legend('Upper Bound', 'ThreshGradDesc');

