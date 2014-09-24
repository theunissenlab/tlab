function test_rbf_grad()

    N = 5000;
    ndim = 1;
    
    nb = 10;
    xmin = -5;
    xmax = 5;
    
    %% create model
    ctrs = linspace(xmin, xmax, nb)';
    modelParams = rbfInit(@exp_bfunc, ndim, nb);
    
    modelParams.w1 = randn(1, nb);
    modelParams.centers = ctrs;

    %% generate data
    global globDat;
    datIdx = 1:N;
    x = randn(N, ndim);
    strfData(x, 1:N, datIdx);    
    [modelParams, y] = rbfFwd(modelParams, datIdx);
    strfData(x, y, datIdx);
    %{
    figure; hold on;
    tmpw1 = modelParams.w1;
    modelParams.w1 = ones(size(modelParams.w1));
    xtmp = linspace(xmin, xmax, N)';
    strfData(xtmp, 1:N, datIdx);
    [modelParams, ytmp] = rbfFwd(modelParams, datIdx);
    strfData(xtmp, ytmp, datIdx);
    modelParams.w1 = tmpw1;
    subplot(2, 1, 1); hold on;
    plot(xtmp, ytmp, 'b-'); axis tight;
    title('basis');
    
    strfData(xtmp, y, datIdx);
    [modelParams, ytmp] = rbfFwd(modelParams, datIdx);
    subplot(2, 1, 2); hold on;
    plot(xtmp, ytmp, 'k-'); axis tight;
    strfData(x, y, datIdx);
    title('output nl');
    %}
    
    %% test gradient
    modelParams.w1 = randn(1, nb);
    
    %% compute gradient
    [modelParams, g] = rbfGrad(modelParams, datIdx);

    %% compute FD approximation to gradient
    [modelParams, gfd] = rbfGradFD(modelParams, datIdx);
    
    g
    gfd
    
end



function [y, dy] = exp_bfunc(r)

    y = exp(-(r.^2));    
    dy = y; %wrong, needs to be fixed
end
    
    


