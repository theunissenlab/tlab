function ch2=CHI2specV2(x,spec)
alpha0=x(1);
beta0=x(2);
fcenter1=x(3);
fwidth1=x(4);
fcenter2=x(5);
fwidth2=x(6);
fcenter3=x(7);
fwidth3=x(8);

gamma=24000;
ttime=linspace(0,0.20,10);
alpha=alpha0+0*ttime;
beta=beta0+0*ttime;           
y=smBGAs(0.20,gamma,alpha',ttime',beta',ttime');
y=y(end/2:end);
fy=timefreq(y,44100,'stft');
ffy=mean(fy.spec,2);
ffzy=ffy.*(ZBfilter(fy.f,fcenter1,fwidth1)+ZBfilter(fy.f,fcenter2,fwidth2)+ZBfilter(fy.f,fcenter3,fwidth3))'/3;
ffzy=ffzy/max(ffzy);
ch2=sum((ffzy-spec).^2);
plot(fy.f,ffzy,fy.f,spec);
drawnow;
