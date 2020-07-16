function [psths,stims,stim_dir]=GetUnitPsth(unit_class,tmin,tmax,tbin)
% Return all psths and stims of a given unit
% 
% INPUT:
%           unit_class : unit should contain list of cells with stims
%
%           tmin: time in ms start of psth recording (before start of stims)
%
%           tmax : end in ms of the psth recording (after the start of stims)
%
%           tbin: bin time in ms (1=1ms 2=2ms etc ...) for the binning 
%
% OUTPUT:
%           psths = m by n matrix ; one rows for each stimuli m = length(unit_class)
%           n= number of bins ceil(tmax+tmin)/tbin
%
%           stims: the m by 1 vector of original stim name associated with the psths
%            vector
%
%

nl=length('/auto/fdata/solveig/ELECTROPHY/');
nstims=length(unit_class);
nbin=ceil((tmax+tmin)/tbin); 
psths=zeros(nstims,nbin);
stims={};
for ns=1:nstims;
    stim=unit_class{ns};
    ntrials=length(stim.trials);
    psth = zeros(1, nbin);   
    owave=stim.original_wavfile;
    [path,name,ext]=fileparts(owave);
    for it=1:ntrials
        trial =stim.trials{it};
        spike_time_trials = trial.spikeTimes;
        nspike = length(spike_time_trials);
        spike_array = zeros(1, nbin);    
        for is=1:nspike
            time_ind = round(spike_time_trials(is)*1000/tbin)+tmin;
            if (time_ind < 1 || time_ind > nbin)
                % spike ignored
            else
                spike_array(time_ind) = spike_array(time_ind) +1;
            end;
        end
    psth = psth + spike_array;
    end
    psths(ns,:)=psth/ntrials;
    stims{ns}=name;
    stim_dir=path((nl+1):end);
end;

