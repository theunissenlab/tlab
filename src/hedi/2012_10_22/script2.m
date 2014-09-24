
fband=125;
nstd=6;
sampleRate=44100.0;

twindow = nstd/(fband*2.0*pi);   % Window length
winLength = fix(twindow*sampleRate);  % Window length in number of points
winLength = fix(winLength/2)*2 % Enforce even window length
increment = fix(0.001*sampleRate)
iter = 0;                               % First time through
cInitialPhase = 2;                      % Linearly weighted cross-correlation
%cInitialPhase = 0;                     % Zero Phase at start for each frame.
cNoCorrelation = 0;                     % Don't do any cross-correlation
cFFTShift = 0;  
signal=wavread('01M05A_281106_001.wav');
[s, t0, f0, pg] = GaussianSpectrum(signal, increment, winLength, sampleRate); 
ytrunc = InvertAndAdd(s, increment, winLength);
plot(1:length(ytrunc),ytrunc,1:length(signal),signal);
z=smBGAsl(0.14,24000,alpha,ttime2,beta,ttime2);
[s, t0, f0, pg] = GaussianSpectrum(z, increment, winLength, sampleRate); 
ss=0*s;

for i=1:length(t0);
    t=t0(i);
    if t>xpas(1,1)/1000.0;
        fc=interp1(xpas(:,1)/1000.0,xpas(:,4),t);
        fw=interp1(xpas(:,1)/1000.0,xpas(:,5),t);    
        %fw=1000.0;
        ss(:,i)= s(:,i).*ZBFilter(f0,fc,fw);
    else
        ss(:,i)= 0*f0;
    end;
end;

ytrunc = InvertAndAdd(ss, increment, winLength);
wavwrite(ytrunc,44100,'yop.wav')
%signal=real(InvertAndAdd(spectrum, increment, winLength, ...
                                        %iter, cFFTShift, cNoCorrelation));
%sign05=SpectrumInversion(t05.spec,increment,winLength,sampleRate);
