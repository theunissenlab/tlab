function wvf=ComputeWaveFromMoto(moto,spec);
% version wit alpha; mean gaussian + si1 + si2 + amplitude spec

slope=-0.6025;
intercept=-0.5979;

time=moto(:,1);
alpha=moto(:,2);
beta=slope*alpha+intercept;

fband=125;
nstd=6;
sampleRate=44100.0;
twindow = nstd/(fband*2.0*pi);   % Window length
winLength = fix(twindow*sampleRate);  % Window length in number of points
winLength = fix(winLength/2)*2; % Enforce even window length
increment = fix(0.001*sampleRate);

ttime=(time-1)/1000.;
tfinal=ttime(end);
z=smBGAsl(tfinal,24000,alpha,ttime,beta,ttime);
[s, t0, f0, pg] = GaussianSpectrum(z, increment, winLength, sampleRate); 
ss=0*s;

final_filter=[];
for i=1:length(t0);
    t=t0(i);
    spi=sum(spec.spec(:,i));
    mu=moto(i,3);    
    sigma1=moto(i,4);
    sigma2=moto(i,5);
    mui=max(mu-f0,0).^2/(2*sigma2^2)+max(f0-mu,0).^2/(2*sigma1^2);
    filter=exp(-mui);
    ssi=sum(abs(s(:,i)).*(filter));
    ss(:,i)= s(:,i).*(filter)/ssi*spi;
end;
wvf = InvertAndAdd(ss, increment, winLength);