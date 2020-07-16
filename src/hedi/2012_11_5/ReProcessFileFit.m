function res=ReProcessFileFit(fname,a,thr);
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
X0=[a 2173.878558 1.182190 3903.084563 14.587081 6873.080217 0.631122];
res={};

hashname=sprintf('../output/pAv3.%s.mat',name);

if ~exist(hashname, 'file')
   fprintf(1,'skipping\n');
else
    load(hashname);
    
    pax=sxpa;
    pax(:,2)=min(-0.40,pax(:,2));
    wavname=sprintf('../output/pAv4.%s.wav',name);
    [final_signal,final_spectrum,raw_spectrum,raw_signal]=GetSignalSpectrumV2(spec,pax,slope,intercept);
    
    wavwrite(final_signal/max(abs(final_signal)),samplerate,wavname);
    resname=sprintf('../output/res.pAv4.%s.mat',name);
    res.data=pax;
    res.final=final_signal;
    res.spectrum=final_spectrum;
    res.raw_signal=raw_signal;
    save(resname,'res');
end;