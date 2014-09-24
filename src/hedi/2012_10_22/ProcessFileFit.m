function res=ProcessFileFit(fname,a,b);
% this is the main function to retrieve the mechanical fit 
% fname = name of the wave file 
% return res as a structure 
% res.alpha, res.beta, res.fc and res.fw are the fitted parameters
% during the time course res.t
% res.sfinal is the final signal
% res.s is the signal unfiltered 
% res.spec is the final spectrum (timefreq's) 
% res.z is the complex spectrum 
% load the wave
addpath('../../data');
x=wavread(fname);
samplerate=44100.0; % always assumed 
z=timefreq(x,samplerate,'stft');

spec=z.spec;
X0=[a,b, 100,100,3500,2000,7000,500];
res={};

hashname=sprintf('pA.%s.mat',fname);

if ~exist(hashname, 'file')
    pax=GetFitCall(hashname,spec,z.t,z.f,X0); % very long unless hashname exist
else
    load(hashname);
    pax=sxpa;
end;

[final_signal,final_spectrum,raw_spectrum,raw_signal]=GetSignalSpectrum(pax);


