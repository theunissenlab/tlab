function res=ProcessFileStressFit(fname,a);
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
[x,samplerate,nbits]=wavread(fname);
n=size(x);
if n(2)>1;
    x=(x(:,1)+x(:,2));
end;
%samplerate=44100.0; % always assumed 
z=timefreq(x,samplerate,'stft');
ev=sum(z.spec);
z.t=z.t(ev>max(ev)/10);
z.spec=z.spec(:,ev>max(ev)/10);

slope=-0.6025;
intercept=-0.5979;

[pathstr,name,ext]=fileparts(fname);


spec=z.spec;
%-0.541404 -0.493281 413.751981 -1.001849 3386.176542 4.774814 4447.570730 -0.826540
X0=[a 1173.878558 1.182190 3903.084563 14.587081 6873.080217 0.631122];
%X0=[-0.576718 -0.539373 1158.551445 -0.488003 3180.517998 4.824780 4649.278908 -1.191427];
res={};

hashname=sprintf('%s/output/st.%s.mat',pathstr,name);

if ~exist(hashname, 'file')
    pax=GetFitCall(hashname,spec,z.t,z.f,X0,slope,intercept,0); % very long unless hashname exist
else
    load(hashname);
    pax=sxpa;
end;
wavname=sprintf('%s/output/st.%s.wav',pathstr,name);
[final_signal,final_spectrum,raw_spectrum,raw_signal]=GetSignalSpectrum(pax,slope,intercept);
wavwrite(final_signal/max(abs(final_signal)),samplerate,wavname);
resname=sprintf('%s/output/res.st.%s.mat',pathstr,name);
res.data=pax;
res.final=final_signal;
res.spectrum=final_spectrum;
res.raw_signal=raw_signal;
save(resname,'res');


