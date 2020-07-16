function mdata = get_method_data(dataRootDir, stimsDir, cellName)

    cellDir = fullfile(dataRootDir, cellName);
    outputDir = fullfile(cellDir, 'conspecific', 'output');
    
    srdata = get_sr_data(dataRootDir, stimsDir, cellName);
    
    descs = {'df', 'tg', 'lars'};
    
    fmts = {'strflab.stft.directfit.nlType_%s.mat', 'strflab.stft.threshgrad.best.nlType_%s.mat', 'strflab.stft.LARS.best.nlType_%s.mat'};
    sfx = {'ES', 'ES', 'Valid'};
        
    hparamNames = {'', 'threshold', 'lambda2'};
    nlTypes = {'linear', 'exponential'};
        
    mdata = struct;
    mdata.resp = srdata.resp;
    mdata.respInfo = srdata.respInfo;
    mdata.groupIndex = srdata.groupIndex;
    mdata.respHalf1 = srdata.splitResp(1, :);
    mdata.respHalf2 = srdata.splitResp(2, :);
    mdata.avgNumTrials = mean(srdata.respInfo.numTrials);
    mdata.stimInfo = srdata.stimInfo;
        
    for m = 1:length(descs)
        
        dname = descs{m};
        dstruct = struct;
        dFmt = fmts{m};
        
        rsfx = sprintf('resp%s', sfx{m});
        cbsfx = sprintf('cBound%s', sfx{m});
        csfx = sprintf('c%s', sfx{m});
        
        for k = 1:length(nlTypes)
            
            nstruct = struct;
        
            mfile = fullfile(outputDir, sprintf(dFmt, nlTypes{k}));
            mvars = load(mfile);            

            mparams = mvars.modelParamsTrained;
            nstruct.bias = mparams.b1;
            nstruct.strf = squeeze(mparams.w1);

            nstruct.respTrain = mvars.respTrain;
            nstruct.respValid = mvars.(rsfx);
            
            nstruct.hparam = -1;
            if ~isempty(hparamNames{m})
                nstruct.hparam = mvars.hyperParams.(hparamNames{m});
            end

            mdata.infobndTrain = mvars.cBoundTrain.info;
            mdata.infobndValid = mvars.(cbsfx).info;
            mdata.tGroups = mvars.trainingGroups;
            mdata.vGroups = mvars.validationGroups;
            
            nstruct.infoTrain = mvars.cTrain.info;
            nstruct.infoValid = mvars.(csfx).info;

            nstruct.perf = nstruct.infoValid / mdata.infobndValid;
            
            nstruct.etime = mvars.elapsedTime;
                        
            dstruct.(nlTypes{k}) = nstruct;
            
            clear mvars;
            
        end
        
        mdata.(dname) = dstruct;
        
    end
    
    