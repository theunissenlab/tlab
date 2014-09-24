function [stim_mean_rate, spike_std_rate, spike_rate_stim] = meanrate_cal(begInd, endInd, samprate,response)
% here Section_Trials is now transformed from s to ms
%% This function calculates the PSTH, isolate arrival times of spikes...
...in reference to the begining of the section (begInd).

%[Section_Trials,psth, stim_mean_rate, spike_std_rate, spike_zscore,
%pvalue, tvalue, pvalue_rand, spike_rate_stim] = meanrate_cal(begInd,
%endInd, samprate,response,spike_rate_bg, Win) put this back if the other
%values are calculated

stimdur=(endInd-begInd)/samprate; %Sould we correct for the 70ms silence periods we added before (10) and after the stim (60)
%ntimebins = round(1000*(endInd-begInd)/samprate); % Number of time bins in ms
%psth = zeros(1, ntimebins);
%t=(1:ntimebins)./1000.0;
ntrials=length(response.trials);
%Section_Trials=cell(ntrials,1);
spike_count = zeros(1,ntrials);
for it=1:ntrials
    trial = response.trials{it};
    spike_time_trials = trial.spikeTimes;
    spike_time1=(spike_time_trials>(begInd/samprate));
    spike_time2=(spike_time_trials<(endInd/samprate));
    SpikesIndices=find(spike_time1 .* spike_time2);
    spike_count(it)=length(SpikesIndices);
%    spike_time_trials_section=spike_time_trials(SpikesIndices);
%    Section_Trials{it}=(spike_time_trials_section-begInd/samprate).*1000; % here spike times arrival are referred to the begining of the section instead of the begining of the stim and in ms instead of s 
%    ns = length(spike_time_trials_section);
%    spike_array = zeros(1, ntimebins);
%    for is=1:ns
%        time_ind = ceil(spike_time_trials_section(is)*1000-begInd/samprate*1000);
%        if (time_ind < 1 || time_ind > ntimebins)
%            fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', response.number, it, time_ind, ntimebins);
%            continue;
%        end
%        spike_array(time_ind) = spike_array(time_ind) +1;
%    end
%    psth = psth + spike_array;
end
%psth = psth./ntrials;


stim_mean_rate = mean(spike_count)./(stimdur);
spike_std_rate = std(spike_count)./(stimdur);
spike_rate_stim = spike_count./stimdur;
%drate = spike_rate_stim - spike_rate_bg;
%rate_rand=[spike_rate_stim spike_rate_bg];
%rate_rand=rate_rand(randperm(2*ntrials));
%drate_rand=rate_rand(1:ntrials)-rate_rand(ntrials+1:2*ntrials);
%mdrate = mean(drate);
%mdrate_rand = mean(drate_rand);
%if (mdrate == 0 )
%    spike_zscore = 0;
%    spike_zscore_rand = 0;
%else
%    spike_zscore = mdrate./std(drate);
%    spike_zscore_rand = mdrate_rand./std(drate_rand);
%end
%tvalue = spike_zscore*sqrt(ntrials);
%tvalue_rand = spike_zscore_rand*sqrt(ntrials);
        
%% Two tailed t-test pvalue
%pvalue = 2.0.*(1-tcdf(abs(tvalue),ntrials-1));
%pvalue_rand = 2.0.*(1-tcdf(abs(tvalue_rand),ntrials-1));
end