function sps=GetAllSpectrum(stims,stim_dir);
% Return the spectrums of all stims as output by timefreq
% INPUT:
%       stims: the vector of all stims used nstims=length(stims)
%       stim_dir: where the wave files are 
%       
% OUTPUT:
%       sps:
%          .specs: a nstims*nspec matrix nspec=ntime*nfrec of the abs(spectrum) obtained via timefreq
%          .nms: length of the longest call in ms 
%          .nfreq: nb frequencies 
%          .f: frequency vectors 


% first compute the longest call
nms=GetLongest(stims,stim_dir);
nfreq=60; % from timefreq default 

% prepare storage 
nstims=length(stims);
specs=zeros(nstims,nfreq*nms); % we pad with zeros at the end 
sps={};
sps.nms=nms;
sps.nfreq=nfreq;

% loop on all stims 

for ns=1:nstims;
    % load wave
    wave_name=strcat(stim_dir,'/',stims(ns),'.wav');
    [wave,sampleRate]=wavread(wave_name{1});
    % compute spectrum 
    tz=timefreq(wave,sampleRate,'stft');
    nt=length(tz.t);
    % reshape into 1d vector 
    spec=reshape(tz.spec,1,nfreq*nt);
    % store it 
    specs(ns,1:(nfreq*nt))=spec;
    
    % store frequency vector 
    if ns==1;
        sps.f=tz.f;
    end;
        
end;
sps.specs=specs;