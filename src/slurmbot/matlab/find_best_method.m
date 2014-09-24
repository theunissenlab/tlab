function [linScores, expScores] = find_best_method(dataDir, cellName, display)

    if nargin < 3
        display = 0;
    end

    outputDir = fullfile(dataDir, cellName, 'conspecific', 'output');
    
    
    descs = {'DF', 'Thresh', 'LARS'};
    
    expOutputFiles = {'strflab.stft.directfit.nlType_exponential.mat', ...
                      'strflab.stft.threshgrad.best.nlType_exponential.mat', ...
                      'strflab.stft.LARS.best.nlType_exponential.mat'};
    
    linOutputFiles = {'strflab.stft.directfit.nlType_linear.mat', ...
                      'strflab.stft.threshgrad.best.nlType_linear.mat', ...
                      'strflab.stft.LARS.best.nlType_linear.mat'};
                  

    expScores = zeros(1, length(expOutputFiles));
    linScores = zeros(1, length(linOutputFiles));
    
    expStrfs = cell(1, length(expOutputFiles));    
    linStrfs = cell(1, length(linOutputFiles));
    expBias = zeros(1, length(expOutputFiles));
    linBias = zeros(1, length(linOutputFiles));
    
    for k = 1:length(expOutputFiles)       
        fname = fullfile(outputDir, expOutputFiles{k});
        evars = load(fname);
        
        eperf = evars.esPerfRatio;
        tperf = evars.trainingPerfRatio;
        
        score = 0.30*tperf + 0.70*eperf;
        expScores(k) = score;
        
        expStrfs{k} = evars.modelParamsTrained.w1;
        expBias(k) = evars.modelParamsTrained.b1;
        
        clear evars;
    end
    
    for k = 1:length(linOutputFiles)       
        fname = fullfile(outputDir, linOutputFiles{k});
        evars = load(fname);
        
        eperf = evars.esPerfRatio;
        tperf = evars.trainingPerfRatio;
        
        score = 0.30*tperf + 0.70*eperf;
        linScores(k) = score;
        
        linStrfs{k} = evars.modelParamsTrained.w1;        
        linBias(k) = evars.modelParamsTrained.b1;
        
        clear evars;
    end
    
    if display
        figure('Name', cellName); hold on;
        for k = 1:length(linOutputFiles)

            strf = linStrfs{k};
            sp = (k-1)*2 + 1;
            subplot(3, 2, sp); hold on;        
            imagesc(squeeze(strf)); axis tight;
            absmax = max(max(abs(strf)));
            caxis([-absmax absmax]);
            colorbar;
            title(sprintf('%s: score=%f, bias=%f', descs{k}, linScores(k), linBias(k)));                
        end

        for k = 1:length(expOutputFiles)

            strf = expStrfs{k};
            sp = k*2;
            subplot(3, 2, sp); hold on;        
            imagesc(squeeze(strf)); axis tight;
            absmax = max(max(abs(strf)));
            caxis([-absmax absmax]);
            colorbar;
            title(sprintf('%s: score=%f, bias=%f', descs{k}, expScores(k), expBias(k)));

        end
    end
       
