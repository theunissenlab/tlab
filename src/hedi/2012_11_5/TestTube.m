function TestTube(tau,alpha)
sm=2*44100;
g=HarmSignals(800,sm);
z=timefreq(g,sm,'stft');
f=OneTube(g,tau,alpha);

zf=timefreq(f,sm,'stft');
plot(zf.f,sum(zf.spec,2)/max(sum(zf.spec,2)),z.f,sum(z.spec,2)/max(sum(z.spec,2)))