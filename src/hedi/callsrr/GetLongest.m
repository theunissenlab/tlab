function nl=GetLongest(stims,stim_dir)
% Compute the longest call in ms in the list sent
% INPUT: 
%        stims: array of nstims wavnames
%        stim_dir: directory of all calls 
% OUTPUT: nl longest call in ms

nl=-1;
nstims=length(stims);
for ns=1:nstims;
    wave_name=strcat(stim_dir,'/',stims(ns),'.wav');
    %wave_name{1}
    [wave,sampleRate]=wavread(wave_name{1});
    ns=ceil(length(wave)/sampleRate*1000)+1; % add one to be sure :) 
    if nl <ns;
        nl=ns;
    end;
end;


