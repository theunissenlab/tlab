function onecell_data = lnl_comp_onecell(cellName, preprocTypes, stimTypes)

    cellDir = fullfile('/auto/fdata/mschachter/data', cellName);
        
    onecell_data = struct;
    onecell_data.cellName = cellName;
    onecell_data.preprocTypes = preprocTypes;
    onecell_data.stimTypes = stimTypes;
    for m = 1:length(stimTypes)
        
        stimType = stimTypes{m};
        outputDir = fullfile(cellDir, stimType, 'output');
        
        for k = 1:length(preprocTypes)

            preprocType = preprocTypes{k};

            outputFileLin = sprintf('strflab.tfType_%s.%s.mat', preprocType, stimType);
            if ~strcmp('surprise', preprocType)           
                outputFileNonlin = sprintf('strflab.staggard.tfType_%s.%s.mat', preprocType, stimType);                
            else
                outputFileNonlin = sprintf('strflab.staggard.tfType_stft.surprise.%s.mat', stimType);
            end

            varsLin = load(fullfile(outputDir, outputFileLin));
            varsNonlin = load(fullfile(outputDir, outputFileNonlin));
            
            vindx = varsNonlin.validationIndex;

            %% get strfs
            odata.lin.strf = varsLin.modelParamsTrained.w1;
            odata.lin.bias = varsLin.modelParamsTrained.b1;
            
            odata.nl.origStrf = varsNonlin.optOptions.diagnostics.linModels{1}.w1;
            odata.nl.origBias = varsNonlin.optOptions.diagnostics.linModels{1}.b1;
            
            odata.nl.strf = varsNonlin.modelParamsTrained.linModel.w1;
            odata.nl.bias = varsNonlin.modelParamsTrained.linModel.b1;
            
            odata.nl.bestIter = varsNonlin.optOptions.diagnostics.bestIteration;

            %% get output NL range
            stim = varsNonlin.preprocData.stim;
            resp = varsNonlin.preprocData.resp;
            model = varsNonlin.modelParamsTrained.linModel;
            groupIndex = varsNonlin.preprocData.groupIndex;

            %% convert to surprise stim if that's what we're dealing with
            if strcmp('surprise', preprocType)           
                wholeStim = varsNonlin.preprocData.stim;
                groupIndex = varsNonlin.preprocData.groupIndex;
                surpriseParams = struct;
                surpriseParams.outputPath = fullfile(cellDir, stimType, 'preproc');
                surpriseParams.outputDesc = 'default';
                [surpriseStimLouder, surpriseStimQuieter, groupIndex, surpriseStimInfo, surpriseParams] = preprocSoundSurprise(wholeStim, groupIndex, surpriseParams);
                stim = [surpriseStimLouder; surpriseStimQuieter]';            
            end

            %% get actual fit output nonlinearity
            strfData(stim, resp, groupIndex);        
            [model, mresp] = strfFwd(model, 1:length(groupIndex));
            mresp(isnan(mresp)) = 0;
            maxResp = max(mresp);        
            minResp = min(mresp);
            edist = 0.05*(maxResp - minResp);

            odata.nl.maxLinResp = maxResp + edist;
            odata.nl.minLinResp = minResp - edist;

            odata.nl.b = varsNonlin.modelParamsTrained.nlModel.b;
            odata.nl.c = varsNonlin.modelParamsTrained.nlModel.c;

            %% get binned static output nonlinearity
            [xvals, yvals, stds] = get_binned_outputnl(mresp, resp, 1);
            odata.nl.binned.x = xvals;
            odata.nl.binned.y = yvals;
            odata.nl.binned.stds = stds;

            odata.realResp = varsNonlin.perfData.wholeResp(vindx);
            odata.modelRespLin = varsLin.perfData.valid.resp;
            odata.modelRespNL = varsNonlin.perfData.valid.resp;
            
            odata.lin.perf = varsLin.perfData.valid.perf;
            odata.nl.perf = varsNonlin.perfData.valid.perf;
            
            odata.infoBound = varsNonlin.perfData.valid.cBound.info;
            
            onecell_data.(preprocType).(stimType) = odata;
            
            
            clear stim;
            clear resp;
            clear model;
            clear groupIndex;
            clear mresp;
            clear model;
            clear xvals;
            clear yvals;
            clear stds;
            clear varsLin;
            clear varsNonlin;

        end
        
    end
    