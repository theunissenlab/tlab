function figure_compare_lyons_lnl(cellNames, stimType)


    dataRootDir = '/auto/fdata/mschachter/data';

    linTrainPerfs = [];
    lnlTrainPerfs = [];
    lyonsTrainPerfs = [];
    
    linValidPerfs = [];
    lnlValidPerfs = [];
    lyonsValidPerfs = [];

    for k = 1:length(cellNames)
       
        cname = cellNames{k};
        cellDir = fullfile(dataRootDir, cname);
        outputDir = fullfile(cellDir, stimType, 'output');
               
        linFileName = sprintf('strflab.tfType_stft.%s.mat', stimType);
        linFilePath = fullfile(outputDir, linFileName);
        
        lnlFileName = sprintf('strflab.staggard.tfType_stft.%s.mat', stimType);
        lnlFilePath = fullfile(outputDir, lnlFileName);

        lyonsFilePath = get_best_lyons(outputDir);
        
        linVars = load(linFilePath);
        linTrainPerf = linVars.perfData.train.perf;
        linValidPerf = linVars.perfData.valid.perf;
        clear linVars;
        
        lnlVars = load(lnlFilePath);
        lnlTrainPerf = lnlVars.perfData.train.perf;
        lnlValidPerf = lnlVars.perfData.valid.perf;
        clear lnlVars;
        
        lyonsVars = load(lyonsFilePath);
        lyonsTrainPerf = lyonsVars.perfData.train.perf;
        lyonsValidPerf = lyonsVars.perfData.valid.perf;
        clear lyonsVars;
        
        linTrainPerfs(end+1) = linTrainPerf;
        linValidPerfs(end+1) = linValidPerf;
        lnlTrainPerfs(end+1) = lnlTrainPerf;
        lnlValidPerfs(end+1) = lnlValidPerf;
        lyonsTrainPerfs(end+1) = lyonsTrainPerf;
        lyonsValidPerfs(end+1) = lyonsValidPerf;        
    end

    trainDiffs = lyonsTrainPerfs - lnlTrainPerfs;    
    validDiffs = lyonsValidPerfs - lnlValidPerfs;
    
    figure; hold on;
    
    subplot(3, 2, 1);
    hist(linTrainPerfs, 15);
    title(sprintf('Lin Train | mean=%0.3f +/- %0.4f', mean(linTrainPerfs), std(linTrainPerfs)));
    
    subplot(3, 2, 2);
    hist(linValidPerfs, 15);
    title(sprintf('Lin Validation | mean=%0.3f +/- %0.4f', mean(linValidPerfs), std(linValidPerfs)));
    
    subplot(3, 2, 3);
    hist(lnlTrainPerfs, 15);
    title(sprintf('LNL Train | mean=%0.3f +/- %0.4f', mean(lnlTrainPerfs), std(lnlTrainPerfs)));
    
    subplot(3, 2, 4);
    hist(lnlValidPerfs, 15);
    title(sprintf('LNL Validation | mean=%0.3f +/- %0.4f', mean(lnlValidPerfs), std(lnlValidPerfs)));
    
    subplot(3, 2, 5);
    hist(lyonsTrainPerfs, 15);
    title(sprintf('Lyons Train | mean=%0.3f +/- %0.4f', mean(lyonsTrainPerfs), std(lyonsTrainPerfs)));
    
    subplot(3, 2, 6);
    hist(lyonsValidPerfs, 15);
    title(sprintf('Lyons Validation | mean=%0.3f +/- %0.4f', mean(lyonsValidPerfs), std(lyonsValidPerfs)));
    
    
    
    h1 = figure; hold on;
    subplot(2, 1, 1);    
    hist(trainDiffs, 15);
    title(sprintf('Training | mean=%0.3f +/- %0.4f', mean(trainDiffs), std(trainDiffs)));
    
    [h, p, ci] = ttest(trainDiffs);
    fprintf('training | mean=%0.3f +/- %0.4f | sig diff=%d | p=%f | ci=(%f,%f)\n', mean(trainDiffs), std(trainDiffs), h, p, ci(1), ci(2));
    
    subplot(2, 1, 2);
    hist(validDiffs, 15);    
    title(sprintf('Validation | mean=%0.3f +/- %0.4f', mean(validDiffs), std(validDiffs)));
    
    [h, p, ci] = ttest(validDiffs);
    fprintf('validation | mean=%0.3f +/- %0.4f | sig diff=%d | p=%f | ci=(%f,%f)\n', mean(validDiffs), std(validDiffs), h, p, ci(1), ci(2));
        
    %{
    svgFile = 'lyons_vs_lnl.svg';
    svgPath = fullfile('~/berkeley/tlab/trunk/src/slurmbot/docs/outputNL', svgFile);
    plot2svg_2d(svgPath, h1);    
    %}
end


function outputFile = get_best_lyons(outputDir)   
    bfile = fullfile(outputDir, 'best.lyons.txt');    
    fid = fopen(bfile);
    a = textscan(fid, '%s');
    outputFile = a{1}{1};
    outputFile = strrep(outputFile, '.txt', '.mat');
end
    