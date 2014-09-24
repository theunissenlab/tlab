function fig_display_modeldata(mdata)

    strf = mdata.model.strf;
    bias = mdata.model.bias;

    figName = sprintf('%s: %s', mdata.cellName, mdata.name);
    
    psth = mdata.data.wholeResp;    
    modelResp = mdata.model.response;
    modelRespRaw = mdata.model.rawResponse;
    vIndx = mdata.data.validationIndex;
    
    cBound = mdata.coherence.bound;
    cTrain = mdata.coherence.training;
    cValid = mdata.coherence.validation;
    
    perfTrain = cTrain.info_mean / cBound.info_mean;
    perfValid = cValid.info_mean / cBound.info_mean;
    
    figure('Name', ['STRF | ' figName]); hold on;
    cval = max(abs(strf(:)));
    imagesc(strf); axis tight;
    caxis([-cval cval]);
    title(sprintf('bias=%f', bias));
   
    figure('Name', ['Resp | ' figName]); hold on;
    plot(psth(vIndx), 'k-'); axis tight;
    plot(modelResp(vIndx), 'r-'); axis tight;
    legend('PSTH', 'Model');
    
    fprintf('%s\n', figName);
    fprintf('\tUpper Bound: %0.1f bits\n', cBound.info_mean);
    fprintf('\tTraining: %0.1f bits, perf=%0.3f\n', cTrain.info_mean, perfTrain);
    fprintf('\tValidation: %0.1f bits, perf=%0.3f\n', cValid.info_mean, perfValid);
    
    