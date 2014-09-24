function ch2=CHI2specV3(x,spec)
alpha0=x(1);
beta0=x(2);
x1=x(3);
y1=x(4);
x2=x(5);
y2=x(6);
x3=x(7);
y3=x(8);

gamma=24000;
ttime=linspace(0,0.20,10);
alpha=alpha0+0*ttime;
beta=beta0+0*ttime;           
y=smBGAs(0.20,gamma,alpha',ttime',beta',ttime');
y=y(end/2:end);
fy=timefreq(y,44100,'stft');
ffy=mean(fy.spec,2);
filter=abs(interp1([-10 x1 x2 x3 10000],[0 y1 y2 y3 0],fy.f,'cubic'));
ffzy=ffy.*filter';
%ffzy=ffy.*(ZBfilter(fy.f,fcenter1,fwidth1)+ZBfilter(fy.f,fcenter2,fwidth2)+ZBfilter(fy.f,fcenter3,fwidth3))'/3;
ffzy=ffzy/max(ffzy);
ch2=sum((ffzy-spec).^2);
plot(fy.f,ffzy,fy.f,spec);
drawnow;
