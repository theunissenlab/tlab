function compute_coherence_for_unit(unitFile, makePlot)

    if nargin < 2
        makePlot = 0;
    end    

    unit = read_unit_h5file(unitFile, 'a');
    h5 = h5utils();
    fid = h5.open(unitFile, 'a');
    
    for k = 1:length(unit.classes)
       
        %% get class response data
        stimClass = unit.classes{k};        
        resps = unit.class_responses.(stimClass);
        numTrials = zeros(1, length(resps));
        allSpikeTrials = cell(1, length(resps));
        stimLengths = zeros(1, length(resps));
        
        for m = 1:length(resps)                       
            numTrials(m) = length(resps{m}.trials);
            respSpikeTrials = cell(numTrials(m), 1);
            for n = 1:numTrials(m)
               respSpikeTrials{n} = resps{m}.trials{n}.spikeTimes;
            end
            allSpikeTrials{m} = respSpikeTrials;
            stimLengths(m) = resps{m}.stim_duration;
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
                
        %% compute coherence                
        windowSize = 0.250;
        freqCutoff = 500;
        sampleRate = 1 / preprocRespParams.binSize;
        cStruct = compute_coherence_bound(splitPsth(1, :), splitPsth(2, :), avgNumTrials, sampleRate, freqCutoff, windowSize);
        meanSpikeRate = mean(psth) * sampleRate;
        
        if makePlot
            figure; hold on;
            plot(cStruct.f, cStruct.c, 'k-');
            plot(cStruct.f, cStruct.cLower, 'b-');
            plot(cStruct.f, cStruct.cUpper, 'r-');
            plot(cStruct.freqCutoff, 0, 'gx', 'MarkerSize', 10);
            axis([min(cStruct.f) max(cStruct.f) 0 1]);
            title(sprintf('%s: info=%f, up=%f, low=%f', stimClass, cStruct.info, cStruct.infoUpper, cStruct.infoLower));
        end
        
        %% write to unit file
        pathPrefix = sprintf('/class_info/%s', stimClass);
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
    
    h5.close(fid);
    
    