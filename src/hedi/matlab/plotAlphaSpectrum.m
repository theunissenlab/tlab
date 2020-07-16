function plotAlphaSpectrum(be,ga)
t=0.1

spec=[];
for al=-0.2:0.001:-0.1;
    y=smBGA(t,al,be,ga,0,1);
    tfrep=timefreq(y,44100,'stft');
    spec=[spec,sum(tfrep.spec,2)];
end;

imagesc(spec);