function res=ProcessFileFitL(fname,alpha,a,L);
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
%addpath('../data');
x=wavread(fname);
samplerate=44100.0; % always assumed 
z=timefreq(x,samplerate,'stft');
ev=sum(z.spec);
z.t=z.t(ev>max(ev)/10);
z.spec=z.spec(:,ev>max(ev)/10);

slope=-0.6025;
intercept=-0.5979;

[pathstr,name,ext]=fileparts(fname);


spec=z.spec;

X0=[a L];
res={};

hashname=sprintf('../output/pL.%s.mat',name);

if ~exist(hashname, 'file')
    pax=GetFitCallL(hashname,spec,z.t,z.f,X0,slope,intercept); % very long unless hashname exist
else
    load(hashname);
    pax=sxpa;
end;
wavname=sprintf('../output/pL.%s.wav',name);
[final_signal,final_spectrum,raw_spectrum,raw_signal]=GetSignalSpectrum(pax,slope,intercept);
wavwrite(final_signal/max(abs(final_signal)),samplerate,wavname);
resname=sprintf('../output/res.pL.%s.mat',name);
res.data=pax;
res.final=final_signal;
res.spectrum=final_spectrum;
res.raw_signal=raw_signal;
save(resname,'res');



