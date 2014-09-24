function fs=MgramWave(Mgram);
addpath('../2012_11_5');
addpath('../sound');
intercept=-0.5979;
slope=-0.6025;

t=Mgram(:,1);
a=Mgram(:,2);
me=Mgram(:,3);
sm=Mgram(:,4);
pm=Mgram(:,5);

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

alpha = a;
beta=slope*alpha+intercept;
ttime=(t-t(1))/1000.;
tfinal=ttime(end);
z=smBGAsl(tfinal,24000,alpha,ttime,beta,ttime);
[s, t0, f0, pg] = GaussianSpectrum(z, increment, winLength, sampleRate); 
ss=0*s;
%smp=max(spec);
mes=interp1(ttime,me,t0);
sms=interp1(ttime,sm,t0);
pms=interp1(ttime,pm,t0);

final_filter=[];
for i=1:length(t0);
    t=t0(i);
    filter=exp(-(f0-mes(i)).^2/(2*sms(i).^2));
    ss(:,i)= pms(i)*s(:,i).*(filter)/(sum(s(:,i))*sum(filter));
    %u=sum(ss(:,i));    
    %ss(:,i)=ss(:,i);
end;
final_signal = InvertAndAdd(ss, increment, winLength);
fs=final_signal/max(abs(final_signal))*.8;
