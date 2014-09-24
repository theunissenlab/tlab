function preproc_model_response(responseFileTemplate, cellName, label, outputFile)

    [strf, bias, outputNL] = compute_avg_strf_and_nl(responseFileTemplate);

    responseFile = sprintf(responseFileTemplate, 0);
    h5 = h5utils();
    fid = h5.open(responseFile);
    
    preprocFile = h5.get_attr(fid, '/data', 'preproc_file');
    
    trStrs = h5.get_attr(fid, '/model', 'transforms');
    trStrs = strtrim(trStrs);
    if ~isempty(trStrs)
        transforms = regexp(trStrs,',','split');
    else
        transforms = '';
    end
    
    modelType = h5.get_attr(fid, '/model', 'type');
    
    isLyons = ~isempty(strfind(preprocFile, 'lyons'));
    
    %% create strflab model
    sz = size(strf);
    numChannels = sz(1);
    strfLength = sz(2);
    strfDelays = 0:(strfLength-1);
    extraParams = struct;
    extraParams.outputNL = outputNL;
    if strcmp('linear', modelType)
        modelParams = linInit(numChannels, strfDelays, modelType);
    elseif ismember(modelType, {'poisson', 'binomial', 'leglm'})        
        extraParams.numTrials = 10;
        modelParams = init_glm_model(numChannels, strfDelays, modelType, extraParams);
    elseif ismember(modelType, {'nl_dists', 'nl_spline'})        
        modelParams = init_glm_model(numChannels, strfDelays, 'gaussian', extraParams);
    end
        
    modelParams.w1 = strf;
    modelParams.b1 = bias;
    
    %% go through each stimulus in the preproc file and compute a response
    fid = h5.open(preprocFile);
    if ~exist(outputFile)
        fout = h5.create(outputFile);
    else
        fout = h5.open(outputFile);
    end
    stimMd5s = h5.get_subgroups(fid, '/');
    for k = 1:length(stimMd5s)
        sgrp = sprintf('/%s/spectrogram', stimMd5s{k});        
        spec = h5.get_ds(fid, sgrp);
        %% transform data
        for m = 1:length(transforms)           
            tfunc = sprintf('transform_%s(spec)', transforms{m});
            spec = eval(tfunc);        
        end
        
        sz = size(spec);
        nchans = sz(1);
        nt = sz(2);
        
        %% produce response
        gindx = ones(1, nt);
        strfData(spec', gindx, gindx);
        [modelParams, resp] = strfFwd(modelParams, 1:length(gindx));
        if isLyons        
            sz = size(resp);
            if sz(1) == 1
                resp = [0 resp];
            else
                resp = [0; resp];
            end
        end
        
        %{
        figure; hold on;
        subplot(4, 1, 1);
        imagesc(spec); axis tight;
        subplot(4, 1, 2);
        imagesc(strf); axis tight;
        subplot(4, 1, 3);
        x = linspace(-1, 1, 200);
        plot(x, outputNL(x), 'k-');        
        subplot(4, 1, 4);
        plot(resp, 'r-'); axis tight;        
        %}
        
        ogrp = sprintf('/%s/%s', stimMd5s{k}, cellName);
        h5.set_ds(fout, ogrp, label, resp);
        
    end
    
    h5.close(fout);
    h5.close(fid);
        
   
end
    
function [y, dy] = output_nl_wrapper(x, fnObj)    
    y = fnval(fnObj, x);
    if nargout > 1
        dy = zeros(size(y));
    end
end


function modelParams = init_glm_model(numChannels, strfDelays, modelType, extraParams)

    switch modelType
        case 'gaussian'
            family = init_glm_family_gaussian();
            outputNL = @(x) identity_outputnl(x);            
        case 'binomial'
            numTrials = extraParams.numTrials;
            family = init_glm_family_binomial(numTrials);
            outputNL = @(x) sigmoid_outputnl(x);
        case 'poisson'
            family = init_glm_family_poisson();
            outputNL = @(x) exp_outputnl(x);
        case 'leglm'
            family = init_glm_family_poisson();            
    end
    
    if isfield(extraParams, 'outputNL')
        outputNL = extraParams.outputNL;
    end
    
    dispersion = 1;
    if ~strcmp(modelType, 'leglm')
        modelParams = glmInit(numChannels, strfDelays, family, outputNL, dispersion);
    else
        modelParams = leglmInit(numChannels, strfDelays, family, dispersion);
    end
end

function [y, dy] = exp_outputnl(x)
    y = exp(x);
    if nargout > 1
        dy = exp(x);
    end
end

function [y, dy] = sigmoid_outputnl(x)
    y = 1 ./ (1 + exp(-x));
    if nargout > 1
        dy = y .* (1 - y);
    end
end

function [y, dy] = identity_outputnl(x)
    y = x;
    if nargout > 1
        dy = 1;
    end
end

