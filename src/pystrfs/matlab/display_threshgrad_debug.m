function display_threshgrad_debug(debugFile)

    vars = load(debugFile);
    steps = vars.optOptions.diagnostics.stepSizes;
    clear vars;

    infoFile = sprintf('%s.info.mat', debugFile);
    vars = load(infoFile);
    
    errs = vars.errs;
    infos = vars.infos;
    gradMax = vars.gradMax;
    gradMean = vars.gradMean;
    gradStd = vars.gradStd;    
    gradMaxThresh = vars.gradMaxThresh;
    gradMeanThresh = vars.gradMeanThresh;
    gradStdThresh = vars.gradStdThresh;    
    
    maxStep = max(steps);
    
    
    nw = 20;
    
    %% normalize errs and infos
    nerrs = errs / max(abs(errs));
    [nerrsMin, nerrsMinIndx] = min(nerrs);
    
    ninfos = infos / max(abs(infos));
    [ninfosMax, ninfosMaxIndx] = max(ninfos);
    
    
    
    fprintf('Best iteration: err=%d, info=%d\n', nerrsMinIndx, ninfosMaxIndx);
    
    figure; hold on;
    subplot(4, 1, 1); hold on;
    plot(nerrs, 'r-');
    plot(nerrsMinIndx, nerrsMin, 'ko', 'MarkerSize', 14);    
    plot(ninfos, 'b-');
    plot(ninfosMaxIndx, ninfosMax, 'ko', 'MarkerSize', 14);
    
    subplot(4, 1, 2); hold on;
    plot(runmean(abs(gradMax), nw));
    title('Abs Max Gradient (running mean)');
    
    subplot(4, 1, 3); hold on;
    %errorbar(1:length(gradMean), runmean(gradMean, nw), runmean(gradStd, nw));
    plot(runmean(gradMean, nw));
    title('Abs Mean Gradient (running mean)');
    axis tight;
    
    subplot(4, 1, 4); hold on;
    plot(runmean(steps, 10) / maxStep);
    title('Normalized Step sizes (running mean)');
    axis tight;
    
    %{
    subplot(4, 1, 4); hold on;
    plot(runmean(gradStd, nw));
    title('Std Gradient (running mean)');
    axis tight;
    %}
    
end


function xnew = runmean(x, nw)

    xnew = zeros(size(x));
    
    for k = 1:length(x)
       
        sindx = max(1, k-nw+1);
        xnew(k) = mean(x(sindx:k));
    end
end
    
function xnew = cummax(x)
    xnew = zeros(size(x));
    x = abs(x);
    for k = 1:length(x)
        xnew(k) = max(x(1:k));        
    end
end
