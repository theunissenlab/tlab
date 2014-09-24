function res=ProcessFileFitv3(fname,a,b);
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

spec=z.spec;
%-0.541404 -0.493281 413.751981 -1.001849 3386.176542 4.774814 4447.570730 -0.826540
X0=[a b 413.751981 -1.001849 3386.176542 4.774814 4447.570730 -0.826540];
%X0=[-0.576718 -0.539373 1158.551445 -0.488003 3180.517998 4.824780 4649.278908 -1.191427];
res={};

hashname=sprintf('../output/pAv3.%s.mat',fname);

if ~exist(hashname, 'file')
    pax=GetFitCallv3(hashname,spec,z.t,z.f,X0); % very long unless hashname exist
else
    load(hashname);
    pax=sxpa;
end;
wavname=sprintf('../output/pAv3.%s.mat.wav',fname);
[final_signal,final_spectrum,raw_spectrum,raw_signal]=GetSignalSpectrumv3(pax);
wavwrite(final_signal/max(abs(final_signal)),samplerate,wavname);
resname=sprintf('../output/res.pAv3.%s.mat',fname);
res.data=pax;
res.final=final_signal;
res.spectrum=final_spectrum;
res.raw_signal=raw_signal;
save(resname,'res');


