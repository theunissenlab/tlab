function display_staggard_summary(cellNames, preprocType, toSvg)

    if nargin < 3
        toSvg = 0;
    end

    cinfo = get_cell_info();
    
    regions = {'mld', 'OV', 'L'};
    stimTypes = {'conspecific', 'flatrip', 'songrip'};
    
    allData = struct;
    for m = 1:length(regions)
        for k = 1:length(stimTypes)
            allData.(regions{m}).(stimTypes{k}) = struct;            
            allData.(regions{m}).(stimTypes{k}).lin.trainPerfs = [];
            allData.(regions{m}).(stimTypes{k}).lin.validPerfs = [];
            allData.(regions{m}).(stimTypes{k}).nl.b = [];
            allData.(regions{m}).(stimTypes{k}).nl.c = [];
            allData.(regions{m}).(stimTypes{k}).nl.trainPerfs = [];
            allData.(regions{m}).(stimTypes{k}).nl.validPerfs = [];
            allData.(regions{m}).(stimTypes{k}).nl.bestIters = [];
        end
    end
    
    for k = 1:length(cellNames)
        
        cname = cellNames{k}; 
        reg = cinfo.(cname).region;
        if reg(1) == 'L'
            reg = 'L';
        end
        stypes = cinfo.(cname).stimTypes;        
        nstypes = length(stypes);
        
        for m = 1:nstypes
            
            stimType = stypes{m};
            
            if ismember(stimType, stimTypes)

                outputDir = fullfile(cinfo.(cname).rootDir, stimType, 'output');
                sfile = sprintf('staggard.tfType_%s.%s.txt', preprocType, stimType);
                outputFilePath = fullfile(outputDir, sfile);

                if ~exist(outputFilePath, 'file')
                    fprintf('File does not exist: %s\n', outputFilePath);
                end

                sdata = csvread(outputFilePath);
                bestIter = sdata(2);
                trainPerf = sdata(3);
                validPerf = sdata(4);
                bestb = sdata(5);
                bestc = sdata(6);

                fprintf('%s: b=%0.2f | c=%0.2f | train=%0.3f | valid=%0.3f | iter=%d\n',...
                        cname, bestb, bestc, trainPerf, validPerf, bestIter);

                    
                [linTrainPerf, linValidPerf] = get_lin_perfs(outputDir, stimType, preprocType);
                
                allData.(reg).(stimType).lin.trainPerfs(end+1) = linTrainPerf;
                allData.(reg).(stimType).lin.validPerfs(end+1) = linValidPerf;

                allData.(reg).(stimType).nl.b(end+1) = bestb;
                allData.(reg).(stimType).nl.c(end+1) = bestc;
                allData.(reg).(stimType).nl.trainPerfs(end+1) = trainPerf;
                allData.(reg).(stimType).nl.validPerfs(end+1) = validPerf;
                allData.(reg).(stimType).nl.bestIters(end+1) = bestIter;
                
            end
            
        end
        
    end
    
    all.b = [];
    all.c = [];
    all.trainDiffs = [];
    all.validDiffs = [];
    
    for k = 1:length(regions)        
        reg = regions{k};
        allByReg.(reg).b = [];
        allByReg.(reg).c = [];
        allByReg.(reg).trainDiffs = [];
        allByReg.(reg).validDiffs = [];
    end
        
    
    for k = 1:length(regions)
    
        reg = regions{k};
        
        for m = 1:length(stimTypes)
           
            stimType = stimTypes{m};
            
            b = allData.(reg).(stimType).nl.b;
            c = allData.(reg).(stimType).nl.c;
            bmean = mean(b);
            bstd = std(b);
            cmean = mean(c);
            cstd = std(c);
            
            nlTrainPerfs = allData.(reg).(stimType).nl.trainPerfs;
            nlValidPerfs = allData.(reg).(stimType).nl.validPerfs;
            bestIters = allData.(reg).(stimType).nl.bestIters;
            
            linTrainPerfs = allData.(reg).(stimType).lin.trainPerfs;
            linValidPerfs = allData.(reg).(stimType).lin.validPerfs;
            
            trainDiffs = nlTrainPerfs - linTrainPerfs;
            validDiffs = nlValidPerfs - linValidPerfs;
            
            all.b = [all.b b];
            all.c = [all.c c];
            all.trainDiffs = [all.trainDiffs trainDiffs];
            all.validDiffs = [all.validDiffs validDiffs];
            
            allByReg.(reg).b = [allByReg.(reg).b b];
            allByReg.(reg).c = [allByReg.(reg).c c];
            allByReg.(reg).trainDiffs = [allByReg.(reg).trainDiffs trainDiffs];
            allByReg.(reg).validDiffs = [allByReg.(reg).validDiffs validDiffs];
            
            tdMean = mean(trainDiffs);
            tdStd = std(trainDiffs);
            vMean = mean(validDiffs);
            vStd = std(validDiffs);            
            
            %{
            figname = sprintf('%s | %s', reg, stimType);
            figure('Name', figname); hold on;
            
            subplot(2, 2, 1);
            hist(b);
            title(sprintf('b | mean=%f +/- %f', bmean, bstd));
            subplot(2, 2, 2);
            hist(c);
            title(sprintf('c | mean=%f +/- %f', cmean, cstd));
            
            subplot(2, 2, 3);
            hist(trainDiffs);
            title(sprintf('Training Diffs | mean=%f +/- %f', tdMean, tdStd));
            
            subplot(2, 2, 4);
            hist(validDiffs);
            title(sprintf('Validation Diffs | mean=%f +/- %f', vMean, vStd));
            
            fprintf('Best Iters: reg=%s | stimType=%s:\n', reg, stimType);
            bestIters
            %}
            
        end        
    end
    
    figure('Name', sprintf('All | %s', preprocType)); hold on;
    subplot(2, 2, 1);
    hist(all.b, 20);
    title(sprintf('b | mean=%f +/- %f', mean(all.b), std(all.b)));
    
    subplot(2, 2, 2);
    hist(all.c, 20);
    title(sprintf('c | mean=%f +/- %f', mean(all.c), std(all.c)));
    
    subplot(2, 2, 3);
    hist(all.trainDiffs, 20);
    title(sprintf('trainDiffs | mean=%f +/- %f', mean(all.trainDiffs), std(all.trainDiffs)));
    
    subplot(2, 2, 4);
    hist(all.validDiffs, 20);
    title(sprintf('validDiffs | mean=%f +/- %f', mean(all.validDiffs), std(all.validDiffs)));
    
    for k = 1:length(regions)
        
        reg = regions{k};
        figName = sprintf('%s | %s', reg, preprocType);
        
        h1 = figure('Name', figName); hold on;
        subplot(2, 2, 1);
        hist(allByReg.(reg).b, 15);
        title(sprintf('b | mean=%0.0f +/- %0.0f', mean(allByReg.(reg).b), std(allByReg.(reg).b)));

        cvals = allByReg.(reg).c;
        cvals(cvals < -150) = [];
        
        subplot(2, 2, 2);
        hist(cvals, 15);
        title(sprintf('c | mean=%0.0f +/- %0.0f', mean(cvals), std(cvals)));
        
        [h, p, ci] = ttest(allByReg.(reg).trainDiffs);
        fprintf('%s: training | mean=%f +/- %f | sig diff=%d | p=%f | ci=(%f,%f)\n', reg, mean(allByReg.(reg).trainDiffs), std(allByReg.(reg).trainDiffs), h, p, ci(1), ci(2));
        
        subplot(2, 2, 3);
        hist(allByReg.(reg).trainDiffs, 15);
        title(sprintf('Training | mean=%0.3f +/- %0.4f', mean(allByReg.(reg).trainDiffs), std(allByReg.(reg).trainDiffs)));

        [h, p, ci] = ttest(allByReg.(reg).validDiffs);
        fprintf('%s: validation | mean=%0.3f +/- %0.4f | sig diff=%d | p=%f | ci=(%f,%f)\n', reg, mean(allByReg.(reg).validDiffs), std(allByReg.(reg).validDiffs), h, p, ci(1), ci(2));
        
        subplot(2, 2, 4);
        hist(allByReg.(reg).validDiffs, 15);
        title(sprintf('Validation | mean=%0.3f +/- %0.4f', mean(allByReg.(reg).validDiffs), std(allByReg.(reg).validDiffs)));        
        
        if toSvg           
            svgName = sprintf('summary_data_%s_%s.svg', reg, preprocType);
            svnPath = fullfile('~/berkeley/tlab/trunk/src/slurmbot/docs/outputNL', svgName);
            plot2svg_2d(svnPath, h1);
        end
        
        
    end
    
end

function [trainPerf, validPerf] = get_lin_perfs(outputDir, stimType, preprocType)


    if strcmp('stft', preprocType)
        stftFile = sprintf('strflab.tfType_%s.%s.txt', preprocType, stimType);
        outputFilePath = fullfile(outputDir, stftFile);        
    else
       
        lyonsFile = get_best_lyons(outputDir);
        outputFilePath = strrep(lyonsFile, '.mat', '.txt');
    end
    linData = csvread(outputFilePath);
    trainPerf = linData(2);
    validPerf = linData(3);
end

function outputFile = get_best_lyons(outputDir)   
    bfile = fullfile(outputDir, 'best.lyons.txt');    
    fid = fopen(bfile);
    a = textscan(fid, '%s');
    outputFile = a{1}{1};
    outputFile = strrep(outputFile, '.txt', '.mat');
end
    
    
    
    
    