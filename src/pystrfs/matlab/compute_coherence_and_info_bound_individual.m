function compute_coherence_and_info_bound_individual(unitFile)

    unit = read_unit_file(unitFile);
    h5 = h5utils();
    fid = h5.open(unitFile);
    
    windowSize = 0.250;    
    freqCutoff = 500;
    
    for k = 1:length(unit.classes)
       
        %% get class response data
        stimClass = unit.classes{k};        
        resps = unit.class_responses.(stimClass);        
        numTrials = zeros(1, length(resps));
        allSpikeTrials = cell(1, length(resps));
        stimLengths = zeros(1, length(resps));
        
        for m = 1:length(resps)                       
            numTrials(m) = resps{m}.numTrials;
            allSpikeTrials{m} = resps{m}.spikeTrials;
            stimLengths(m) = resps{m}.trialDuration;
        end
        
        %% create PSTHs
        preprocRespParams = struct;      
        preprocRespParams.units = 's';   
        preprocRespParams.binSize = 0.001;
        preprocRespParams.stimLengths = stimLengths;

        [psth, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);
        
        preprocRespParams.split = 1;
        [splitPsth, groupIndex, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);
        
        avgNumTrials = mean(numTrials);
        
        uniqueGroups = unique(groupIndex);        
        spreads = zeros(length(uniqueGroups), 1);
        infoVals = zeros(length(uniqueGroups), 1);
        for m = 1:length(uniqueGroups)
            gindx = groupIndex == uniqueGroups(m);
            
            respPsth = psth(gindx);
            respSplitPsth = splitPsth(:, gindx);
            %% compute coherence            
            sampleRate = 1 / preprocRespParams.binSize;
            cStruct = compute_coherence_bound(respSplitPsth(1, :), respSplitPsth(2, :), avgNumTrials, sampleRate, freqCutoff, windowSize);
            meanSpikeRate = mean(respPsth) * sampleRate;

            spreads(m) = (cStruct.infoUpper - cStruct.infoLower) / 2;
            infoVals(m) = cStruct.info;
            
            %{
            if strcmp(stimClass, 'Con')
                figure; hold on;
                
                subplot(2, 1, 1); hold on;
                plot(respPsth, 'k-');
                axis([1 length(respPsth) 0 1]);                
                subplot(2, 1, 2); hold on;
                plot(cStruct.f, cStruct.c, 'k-');
                plot(cStruct.f, cStruct.cLower, 'b-');
                plot(cStruct.f, cStruct.cUpper, 'r-');
                plot(cStruct.freqCutoff, 0, 'gx', 'MarkerSize', 10);
                axis([min(cStruct.f) max(cStruct.f) 0 1]);
                title(sprintf('%s, resp %d: info=%f, up=%f, low=%f', stimClass, uniqueGroups(m), cStruct.info, cStruct.infoUpper, cStruct.infoLower));
            end
            %}

            %% write to unit file
            pathPrefix = sprintf('/class_info/%s/responses/%d', stimClass, uniqueGroups(m));

            h5.set_ds(fid, pathPrefix, 'mean_rate', meanSpikeRate);
            h5.set_ds(fid, [pathPrefix '/coherence'], 'cutoff_frequency', cStruct.freqCutoff);
            h5.set_ds(fid, [pathPrefix '/coherence'], 'frequencies', cStruct.f);
            h5.set_ds(fid, [pathPrefix '/coherence'], 'upper', cStruct.cUpper);
            h5.set_ds(fid, [pathPrefix '/coherence'], 'lower', cStruct.cLower);
            h5.set_ds(fid, [pathPrefix '/coherence'], 'mean', cStruct.c);

            h5.set_ds(fid, [pathPrefix '/coherence/info'], 'upper', cStruct.infoUpper);
            h5.set_ds(fid, [pathPrefix '/coherence/info'], 'lower', cStruct.infoLower);
            h5.set_ds(fid, [pathPrefix '/coherence/info'], 'mean', cStruct.info);        

        end
        
        ranks = (infoVals / max(infoVals));
        fprintf('Info Stats for %s: %0.0f +/ %0.1f\n', stimClass, mean(infoVals), std(infoVals));
        fprintf(' Ranks:\n');
        for m = 1:length(uniqueGroups)
            fprintf('    %d: %0.2f\n', uniqueGroups(m), ranks(m));
        end
        fprintf('   Spreads for %s: %0.0f +/ %0.1f\n', stimClass, mean(spreads), std(spreads));
        
        
    end
    
    h5.close(fid);
    
