function [final_signal,final_spectrum,raw_spectrum,raw_signal,final_filter]=GetSignalSpectrumv3(sxpa);


fband=125;
nstd=6;
sampleRate=44100.0;
twindow = nstd/(fband*2.0*pi);   % Window length
winLength = fix(twindow*sampleRate);  % Window length in number of points
winLength = fix(winLength/2)*2; % Enforce even window length
increment = fix(0.001*sampleRate);
iter = 0;                               % First time through
cInitialPhase = 2;                      % Linearly weighted cross-correlation
cNoCorrelation = 0;                     % Don't do any cross-correlation
cFFTShift = 0;  

alpha=myunique(sxpa(:,2));
beta=myunique(sxpa(:,3));
ttime=myunique(sxpa(:,1))/1000.;
tfinal=ttime(end);
ttime=[0;ttime];
alpha=[alpha(1);alpha];
beta=[beta(1);beta];
vfc1=myunique(sxpa(:,4));
vfw1=myunique(sxpa(:,5));
vfc2=myunique(sxpa(:,6));
vfw2=myunique(sxpa(:,7));
vfc3=myunique(sxpa(:,8));
vfw3=myunique(sxpa(:,9));
    

z=smBGAsl(tfinal,24000,alpha,ttime,beta,ttime);
z=z/max(abs(z));
[s, t0, f0, pg] = GaussianSpectrum(z, increment, winLength, sampleRate); 
ss=0*s;
final_filter=[];
for i=1:length(t0);
    t=t0(i);
    if t>sxpa(1,1)/1000.0;
        x1=interp1(ttime(2:end),vfc1,t);
        y1=interp1(ttime(2:end),vfw1,t); 
         x2=interp1(ttime(2:end),vfc2,t);
        y2=interp1(ttime(2:end),vfw2,t); 
         x3=interp1(ttime(2:end),vfc3,t);
        y3=interp1(ttime(2:end),vfw3,t); 
        %fw=1000.0;
        
        filter=abs(interp1([-10 x1 x2 x3 10000 50000],[0 y1 y2 y3 0 0],f0,'cubic'));
        final_filter=[final_filter,filter];
        ss(:,i)= s(:,i).*(filter/sum(filter));
    else
        ss(:,i)= 0*f0;
    end;
end;
raw_signal=z;
raw_spectrum=s;
final_spectrum=ss;
final_signal = InvertAndAdd(ss, increment, winLength);
final_signal=final_signal/max(abs(final_signal));
