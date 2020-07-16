function TestTwoTubes(tau1,alpha1,tau2,alpha2)
sm=2*44100;
g=HarmSignals(800,sm);
z=timefreq(g,sm,'stft');
f=TwoTube(g,tau1,alpha1,tau2,alpha2);

zf=timefreq(f,sm,'stft');
subplot(2,1,1)
plot(zf.f,sum(zf.spec,2))
subplot(2,1,2)
plot(zf.f,sum(zf.spec,2)/max(sum(zf.spec,2)),z.f,sum(z.spec,2)/max(sum(z.spec,2)))