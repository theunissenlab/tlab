function ch2=CHI2spec(x,spec)
alpha0=x(1);
beta0=x(2);
fcenter=x(3);
fwidth=x(4);
gamma=24000;
ttime=linspace(0,0.20,10);
alpha=alpha0+0*ttime;
beta=beta0+0*ttime;           
y=smBGAs(0.20,gamma,alpha',ttime',beta',ttime');
y=y(end/2:end);
fy=timefreq(y,44100,'stft');
ffy=mean(fy.spec,2);
ffzy=ffy.*ZBfilter(fy.f,fcenter,fwidth)';
ffzy=ffzy/max(ffzy);
ch2=sum((ffzy-spec).^2);
plot(fy.f,ffzy,fy.f,spec);
drawnow;
