function unit = read_unit_h5file(unitPath,perm)
% perm should be 'r' for reading only, 'a' to add fields

    h5 = h5utils();

    unit = struct;
    
    fid = h5.open(unitPath, perm);
    
    unit.electrode = h5.get_attr(fid, '/', 'electrode');
    try
        unit.ldepth = h5.get_attr(fid, '/', 'ldepth');
    catch err
        unit.ldepth = [];
    end
    try
        unit.rdepth = h5.get_attr(fid, '/', 'rdepth');    
    catch err
        unit.rdepth = [];
    end
	 try
		  unit.sortType = h5.get_attr(fid, '/', 'sortType');
	 catch err
		 unit.rdepth = [];
    end
    try
	    unit.sortid = h5.get_attr(fid, '/', 'sortid');
	 catch err
		 unit.sortid = [];
	 end

    unit.site = h5.get_attr(fid, '/', 'site');
    unit.source_directory = h5.get_attr(fid, '/', 'source_directory');
    
    hasCoherenceInfo = 0;
    hasExtraInfo = 0;
    
    crList = h5.get_subgroups(fid, '/');
    realClasses = {};
    for k = 1:length(crList)        
        if ~ismember(crList{k}, {'class_info', 'extra_info'})
            realClasses{end+1} = crList{k};
        elseif strcmp(crList{k}, 'class_info')
            hasCoherenceInfo = 1;
        elseif strcmp(crList{k}, 'extra_info')
            hasExtraInfo = 1;
        end
    end
    crList = realClasses;
    unit.classes = realClasses;
    unit.class_responses = struct;
    
    for k = 1:length(crList)

        className = crList{k};
        crPath = sprintf('/%s', className);
        stimNumbers = h5.get_subgroups(fid, crPath);
        
        allStims = {};
        
        for m = 1:length(stimNumbers)
           
            stimPath = sprintf('%s/%s', crPath, stimNumbers{m});
            
            stim = h5.get_attrs(fid, stimPath);
            stim.number = stimNumbers{m};
            trialNames = h5.get_subgroups(fid, stimPath);
            stim.trials = {};
            for j = 1:length(trialNames)
                
               try
                   respPath = sprintf('%s/%s', stimPath, trialNames{j});
                   resp = h5.get_attrs(fid, respPath);
               catch err
                   fprintf('Problem finding responses at path %s\n', respPath);
                   rethrow(err);
               end
               try
                   spikeIdsPath = sprintf('%s/%s', respPath, 'spike_ids');
                   resp.spikeIds = h5.get_ds(fid, spikeIdsPath);
               catch err
                   fprintf('Problem finding spike ids at path %s\n', spikeIdsPath);
                   rethrow(err);
               end
               try
                   spikeTimesPath = sprintf('%s/%s', respPath, 'spike_times');
                   resp.spikeTimes = h5.get_ds(fid, spikeTimesPath);    
               catch err
                   fprintf('Problem finding spike times at path %s\n', spikeTimesPath);
                   rethrow(err);
               end
               if resp.spikeIds == -999
                   resp.spikeIds = [];
                   resp.spikeTimes = [];
               end
               stim.trials{end+1} = resp;                
            end
            
            allStims{end+1} = stim;
        end        
        
        unit.class_responses.(className) = allStims;
        
        if hasCoherenceInfo
            
            cstruct = struct;
            cpath = sprintf('class_info/%s/coherence', className);
            
            cstruct.meanSpikeRate = h5.get_ds(fid, sprintf('class_info/%s/mean_rate', className));
            cstruct.cutoffFreq = h5.get_ds(fid, [cpath '/cutoff_frequency']);
            cstruct.f = h5.get_ds(fid, [cpath '/frequencies']);
            cstruct.c = h5.get_ds(fid, [cpath '/mean']);
            cstruct.cUpper = h5.get_ds(fid, [cpath '/upper']);
            cstruct.cLower = h5.get_ds(fid, [cpath '/lower']);
            cstruct.info = h5.get_ds(fid, [cpath '/info/mean']);
            cstruct.infoUpper = h5.get_ds(fid, [cpath '/info/upper']);
            cstruct.infoLower = h5.get_ds(fid, [cpath '/info/lower']);
            
            unit.class_info.(className) = cstruct;        
        end
        if hasExtraInfo
            if ismember(className, h5.get_subgroups(fid, '/extra_info/'))
                eiList = h5.get_subgroups(fid, sprintf('/extra_info/%s/',className));
                if ismember('sections_data', eiList)
                    sstruct = struct;
                    spath = sprintf('extra_info/%s/sections_data', className);
                
                    sstruct.nb_signif_sections_Bonf = h5.get_ds(fid, [spath '/num_signif_sections_Bonferroni']);
                    sstruct.nb_signif_sections_BonfDeg = h5.get_ds(fid, [spath '/num_signif_sections_BonferroniDeg']);
                    sstruct.sections_pval = h5.get_ds(fid, [spath '/sections_pvalues']);
                
                    unit.extra_info.(className).sections_data = sstruct;
                end
                if ismember('inter_spike_interval', eiList)
                    istruct = struct;
                    ipath = sprintf('extra_info/%s/inter_spike_interval', className);
                
                    istruct.Stim_rate = h5.get_ds(fid, [ipath '/Stim_rate']);
                    istruct.Spike_proba = h5.get_ds(fid,[ipath '/Spike_probability']);
                    istruct.Interval_Inter_Spike = h5.get_ds(fid, [ipath '/Interval_Inter_spike']);
                    istruct.Spike_proba_1ms = h5.get_ds(fid, [ipath '/Spike_proba_1ms']);
                    istruct.Spike_proba_1ms_theo = h5.get_ds(fid, [ipath '/Spike_proba_1ms_theo']);
                
                    unit.extra_info.(className).inter_spike_interval = istruct;
                end
            end
        end
    end
        
    h5.close(fid);
    
    
