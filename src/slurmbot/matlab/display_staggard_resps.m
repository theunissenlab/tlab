function display_staggard_resps(cellName, stimType, preprocType, toSvg)

    if nargin < 4
        toSvg = 0;
    end

    dataRootDir = '/auto/fdata/mschachter/data';
    cellDir = fullfile(dataRootDir, cellName);
    outputDir = fullfile(cellDir, stimType, 'output');
    
    if ~exist(outputDir, 'file')       
        fprintf('Output directory does not exist: %s', outputDir);        
        return;
    end
    
    outputFileName = sprintf('strflab.staggard.tfType_%s.%s.mat', preprocType, stimType);
    outputFilePath = fullfile(outputDir, outputFileName);

        
    if strcmp('stft', preprocType)
        linFileName = sprintf('strflab.tfType_stft.%s.mat', stimType);
        linFilePath = fullfile(outputDir, linFileName);
    else
        linFilePath = get_best_lyons(outputDir);        
    end
    
    %% get performance of normal linear model
    lvars = load(linFilePath);
    linPerfData = lvars.perfData;
    linPerfTrain = linPerfData.train.perf;    
    linPerfValid = linPerfData.valid.perf;
    clear lvars;
        
    %% get parameters of lnl model
    vars = load(outputFilePath);    
    metaOptions = vars.optOptions;
    modelParamsTrained = vars.modelParamsTrained;
    preprocData = vars.preprocData;
    perfData = vars.perfData;
    trainingIndex = vars.trainingIndex;
    validationIndex = vars.validationIndex;
    stimInfo = preprocData.stimInfo;
    
    clear vars;
    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    
    freqCutoff = 90;
    infoWindowSize = 0.250;
    
    strfData(wholeStim, wholeResponse, groupIndex);
    
    %% find minimum error
    [minErr, minErrIter] = min(metaOptions.diagnostics.nlErrs(2:end));
    minErrIter = minErrIter + 1;    
    fprintf('Min NL error at iteration %d, err=%f\n', minErrIter, minErr);
    
    bestIter = metaOptions.diagnostics.bestIteration;
    
    fprintf('Best Model: lsq=%d | info=%d\n', minErrIter, bestIter);
    
    
    %% get best model
    origStrf = metaOptions.diagnostics.linModels{1}.w1;
    origBias = metaOptions.diagnostics.linModels{1}.b1;
    
    bestLinModel = metaOptions.diagnostics.linModels{bestIter};
    bestStrf = bestLinModel.w1;
    bestBias = bestLinModel.b1;
    
    bestNLModel = metaOptions.diagnostics.nlModels{bestIter};
    bestb = bestNLModel.b;
    bestc = bestNLModel.c;
    
    strfDiff = origStrf - bestStrf;
    
    bestLNLModel = lnlInit(bestLinModel, bestNLModel);
    
    
    
    %compute range of output NL
    [bestLinModel, bestLinResp] = strfFwd(bestLinModel, 1:length(wholeResponse));
    minX = min(bestLinResp);
    maxX = max(bestLinResp);    
    
    xinc = (maxX - minX) / 100;
    gx = minX:xinc:maxX;
    gy = gompertz(1, bestb, bestc, gx);
    
    
    
    
    figName = sprintf('%s | %s | %s', cellName, stimType, preprocType);
    
    %% Make figures       
    display = 1;
    if display
        figure('Name', figName); hold on;
        subplot(3, 1, 1); hold on;
        imagesc(origStrf); axis tight;
        absmax = max(max(abs(origStrf)));
        caxis([-absmax absmax]);
        colorbar;
        title(sprintf('First STRF | Bias=%f', origBias));

        subplot(3, 1, 2); hold on;
        imagesc(bestStrf); axis tight;
        absmax = max(max(abs(bestStrf)));
        caxis([-absmax absmax]);
        colorbar;
        title(sprintf('Best STRF | Bias=%f', bestBias));

        subplot(3, 1, 3); hold on;
        imagesc(strfDiff); axis tight;
        absmax = max(max(abs(strfDiff)));
        caxis([-absmax absmax]);
        colorbar;
        title('STRF Difference');

        h2 = figure('Name', figName); hold on;
        plot(gx, gy, 'k-');
        axis tight;
        title(sprintf('Output NL | b=%f | c=%f', bestb, bestc));
        axis tight;
        
        h1 = figure('Name', figName); hold on;
        subplot(4, 1, 1); hold on;
        plot(linPerfData.wholeResp(validationIndex), 'k-');
        title('PSTH');
        axis tight;
        
        subplot(4, 1, 2); hold on;
        plot(linPerfData.valid.resp, 'r-');
        title('Linear Prediction');
        axis tight;
        
        subplot(4, 1, 3); hold on;
        plot(perfData.valid.resp, 'b-');
        title('LNL Prediction');
        axis tight;
        
        rdiff = linPerfData.valid.resp - perfData.valid.resp;
        
        subplot(4, 1, 4); hold on;
        plot(rdiff, 'k-');
        title(sprintf('Difference | mean=%f', mean(rdiff)));
        axis tight;
        
        if toSvg           
            svgFile = sprintf('respinfo_%s_%s_%s.svg', cellName, stimType, preprocType);
            svgPath = fullfile('~/berkeley/tlab/trunk/src/slurmbot/docs/outputNL', svgFile);
            plot2svg_2d(svgPath, h1);
            
            svgFile = sprintf('nl_%s_%s_%s.svg', cellName, stimType, preprocType);
            svgPath = fullfile('~/berkeley/tlab/trunk/src/slurmbot/docs/outputNL', svgFile);
            plot2svg_2d(svgPath, h2);
        end
        
    end
    
    
end

function outputFile = get_best_lyons(outputDir)   
    bfile = fullfile(outputDir, 'best.lyons.txt');    
    fid = fopen(bfile);
    a = textscan(fid, '%s');
    outputFile = a{1}{1};
    outputFile = strrep(outputFile, '.txt', '.mat');
end
    