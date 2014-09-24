function test_leglm_real()

    unitDir = '/auto/k6/mschachter/pystrfs/units';
    preprocDir = '/auto/k6/mschachter/pystrfs/preproc';
    cellName = 'yy1617_4_A';
    
    unitFile = fullfile(unitDir, cellName, 'unit.h5');
    preprocFile = fullfile(preprocDir, 'stft.nstd_6.fband_125.h5');
    stimClass = 'Con';

    stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass);
    
    %% transform stimulus
    wholeStim = stimRespData.wholeStim;
    wholeStim = transform_log(wholeStim);
    wholeStim = transform_zscore(wholeStim);
    
    stimRespData.wholeStim = wholeStim;
    
    numTrials = round(mean(stimRespData.numTrials));
    
    %family = init_glm_family_gaussian();
    %family = init_glm_family_binomial(numTrials);
    family = init_glm_family_poisson();
    
    dispersion = 1;
    
    %% transform response    
    if strcmp(family.type, 'poisson')       
        wholeResp = stimRespData.wholeResp;
        gindx = stimRespData.groupIndex;
        grps = unique(gindx);
        for k = 1:length(grps)           
            gi = gindx == grps(k);
            wholeResp(gi) = wholeResp(gi)*stimRespData.numTrials(k);            
        end
        stimRespData.wholeResp = wholeResp;
    end
    
    %{
    figure; hold on;
    subplot(2, 1, 1); hold on;
    plot(stimRespData.wholeResp); axis tight;
    subplot(2, 1, 2); hold on;
    rindx = stimRespData.wholeResp > 0;
    hist(stimRespData.wholeResp(rindx), 20); axis tight;
    %}
    
    %% Initialize strflab global variables with our stim and responses
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

    
    %% Initialize a linear model
    strfLength = 60;
    strfDelays = 0:(strfLength-1);    
    modelParams = leglmInit(stimRespData.numChannels, strfDelays, family, dispersion);
    
    %% Initialize Threshold Gradient Descent options
    optMeth = 'tg';
    if strcmp(optMeth, 'tg')
        optOptions = trnThreshGradDescLS();        
        optOptions.threshold = [0.50 0];
        optOptions.earlyStop = 1;
        optOptions.display = 1;        
        optOptions.maxLSIter = 100;
        optOptions.maxIter = 500;
        optOptions.errLastN = 20;
        optOptions.errStartN = 2;
        optOptions.stepConvWindow = 8;
        optOptions.iterFunc = @(x, y, z) iterfunc(x, y, z);
        
        pg = ones(1, modelParams.nWts);
        pg(end-1) = 2;
        optOptions.paramGroup = pg;
        
        switch family.type
            case 'poisson'
                optOptions.maxStep = 5;        
                optOptions.minStep = 1e-6;
                optOptions.stepSize = 1e-2;
                optOptions.stepConv = 1e-5;
                
            otherwise
                optOptions.maxStep = 5;        
                optOptions.minStep = 1e-5;
                optOptions.stepSize = 1e-4;
                optOptions.stepConv = 1e-4;
        end        

    elseif strcmp(optMeth, 'scg')
        optOptions = trnSCG();
        optOptions.display = 1;
        optOptions.earlyStop = 0;
        optOptions.maxIter = 75;
        switch family.type
            case 'poisson'
                optOptions.minErrChange = 1;
            otherwise
                optOptions.minErrChange = .1;
        end            
    end
    
    trIndex = findIdx(1:18, stimRespData.groupIndex);
    esIndex = findIdx(19:20, stimRespData.groupIndex);        
    
    %% run direct fit
    if optOptions.earlyStop
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trIndex, optOptions, esIndex);
    else
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trIndex, optOptions);
    end
    
    save(sprintf('/auto/fhome/mschachter/leglm_strflab.%s.%s.mat', family.type, optMeth), 'modelParamsTrained');
    
end

function iterfunc(x, grad, bestStep)    
    gm = grad(end)
    m = x(end)
end
