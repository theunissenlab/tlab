function zUnit = calc_zscore_unit(unit, preTime)
%Calculates the z score and the pvalue for each stimulus in the h5 data
%structure unit.  This function can be called after calling:
%   unit=read_unit_h5file(h5Path)
% preTime is the duration of the background preceeding the stimulus in s
% and is used in the calculation of the background rate.  If preTime is [],
% the function will use the maximum amount of time available according to
% the pre_time field for each trial.

% Loop through all classes and all stims.
nclasses = length(unit.classes);


for ic=1:nclasses
    className =  unit.classes{ic};
    responses = unit.class_responses.(className);
    nstims = length(responses);
    spike_zscore = zeros(nstims, 1);
    tvalue = zeros(nstims,1);
    pvalue = zeros(nstims,1);
    
    for ns=1:nstims
        stimLen = responses{ns}.stim_duration;   %here stimLen is in s        
        Trials = responses{ns}.trials;
        nt = length(Trials);

        spike_count = zeros(1,nt);
        spike_count_bg = zeros(1,nt);
        rate_diff = zeros(1,nt);
        
        totStimLen = 0;
        totBgLen = 0;
        for it = 1:nt
            trial = Trials{it};
            
            if isempty(preTime)
                if isfield(trial, 'pre_time')
                    bgStart = -trial.pre_time;
                else
                    fprintf(1, 'Warning no pre_time in h5 data file: 1 s will be used\n');
                    bgStart = -1.0;
                end
            else
                bgStart = -preTime;
                if isfield(trial, 'pre_time')
                    if preTime > trial.pre_time
                        fprintf(1, 'Warning: trial %d in stim %s has shorter background time (%.3f) than asked (%.3f)\n', it, ns, trial.pre_time, preTime);
                        bgStart = -trial.pre_time;
                    end
                end
               
            end
            % The definition of time intervals matches the one in
            % convert.py
            spike_count(it) = length(find(trial.spikeTimes >= 0.0 & trial.spikeTimes <= stimLen));
            totStimLen = totStimLen + stimLen;
            
            spike_count_bg(it) = length(find(trial.spikeTimes < 0.0 & trial.spikeTimes >= bgStart));
            totBgLen = totBgLen - bgStart;
            
            % There is a minus sign because bgStart is negative...
            rate_stim = spike_count(it)/stimLen;
            rate_back = -spike_count_bg(it)/bgStart;
            
%              if (abs(rate_stim - trial.peri_rate) > 0.001)
%                 fprintf(1,'Rate difference for trial %d: Old StimR=%.4f New StimR=%.4f Old BgR=%.4f New BgR=%.4f\n', it, trial.peri_rate, rate_stim, trial.bg_rate, rate_back);
%              end
            
            unit.class_responses.(className){ns}.trials{it}.peri_rate = rate_stim;
            unit.class_responses.(className){ns}.trials{it}.bg_rate = rate_back;
            
            
            rate_diff(it) = rate_stim - rate_back;
        end
        
        % The spike rates are in spikes per second - this is not saved for
        % now...
        stim_mean_rate = sum(spike_count)/totStimLen;
        back_mean_rate = sum(spike_count_bg)/totBgLen;

        mdrate = mean(rate_diff);
        if (mdrate == 0 )
            spike_zscore(ns) = 0;
        else
            spike_zscore(ns) = mdrate./std(rate_diff);
        end
        tvalue(ns) = spike_zscore(ns)*sqrt(nt);
        
        % Two tailed t-test pvalue
        pvalue(ns) = 2.0.*(1-tcdf(abs(tvalue(ns)),nt-1));
        
%         if (abs(pvalue(ns)-responses{ns}.pvalue) > 0.001)
%             fprintf(1,'Stim %s: Old t= %.4f p= %.4f New t= %.4f p= %.4f\n', responses{ns}.number, responses{ns}.tstat, responses{ns}.pvalue, tvalue(ns), pvalue(ns));
%             changeFound = 1;
% %             pause();
%         end
        unit.class_responses.(className){ns}.pvalue = pvalue(ns);
        unit.class_responses.(className){ns}.zscore = spike_zscore(ns);
        unit.class_responses.(className){ns}.tstat = tvalue(ns);
    end        
end
zUnit = unit;
return






