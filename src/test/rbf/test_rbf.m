function test_rbf()

    N = 10000;
    ndim = 1;
    
    nb = 10;
    xmin = -5;
    xmax = 5;
    
    %% create model
    ctrs = linspace(xmin, xmax, nb)';
    modelParamsReal = rbfInit(@exp_bfunc, ndim, nb);
    
    modelParamsReal.w1 = randn(1, nb);
    modelParamsReal.centers = ctrs;

    %% generate data
    global globDat;
    datIdx = 1:N;
    x = randn(N, ndim);
    strfData(x, 1:N, datIdx);    
    [modelParamsReal, y] = rbfFwd(modelParamsReal, datIdx);
    strfData(x, y, datIdx);
    
     %% create model to fit
    modelParams = rbfInit(@exp_bfunc, ndim, nb);
    modelParams.centers = modelParamsReal.centers;
    modelParams.w1 = randn(1, nb);
            
    %% fit model
    %{
    optOptions = trnThreshGradDescLS();    
    optOptions.threshold = 0;
    optOptions.earlyStop = 0;
    optOptions.display = 1;
    optOptions.maxIter = 2000;    
    optOptions.maxStep = 5;
    optOptions.minStep = 1e-7;
    optOptions.maxLSIter = 50;
    optOptions.stepSize = 1e-4;
    optOptions.stepConvWindow = 8;
    optOptions.stepConv = 5e-6;
    %}
    
    optOptions = trnSCG();
    optOptions.display = 1;
    
    
    %% run training
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
    
    realWeights = modelParamsReal.w1
    trainedWeights = modelParamsTrained.w1
    wdiff = realWeights - trainedWeights
    
    %% plot functions
    npts = 120;
    xrng = linspace(xmin, xmax, npts)';
    strfData(xrng, 1:npts, ones(1, npts));
    [modelParams, yReal] = rbfFwd(modelParamsReal, 1:npts);
    [modelParamsTrained, yTrain] = rbfFwd(modelParamsTrained, 1:npts);
    
    figure; hold on;    
    plot(xrng, yReal, 'k-');
    plot(xrng, yTrain, 'r-');
    axis tight;
    legend('Real', 'Trained');
    
    
end



function [y, dy] = exp_bfunc(r)

    y = exp(-(r.^2));    
    dy = y; %wrong, needs to be fixed
end
    
    


