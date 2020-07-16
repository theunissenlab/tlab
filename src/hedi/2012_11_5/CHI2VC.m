function ch2=CHI2VC(x,spec,slope,intercept,alpha0)
%alpha0=x(1);
beta0=alpha0*slope+intercept;
a=x(1);
L=x(2);

c=340.0; % speed of sound

tmax=0.5;
gamma=24000;
ttime=linspace(0,tmax,10);
alpha=alpha0+0*ttime;
beta=beta0+0*ttime;           
y=smBGAs(tmax,gamma,alpha',ttime',beta',ttime');
y=y(end/2:end);

fy=timefreq(y,44100,'stft');
ffy=mean(fy.spec,2);
fw=c/(4*L)*1000;
tau=2*L/c;
%filter=exp(-(fy.f-fw).^2/(2*a^2));
filter=abs(1./(sqrt(2+a*2*cos(2*pi*tau*fy.f/1000.0))));

%filter=abs(interp1([-10 x1 x2 x3 10000],[0 y1 y2 y3 0],fy.f,'cubic'));
ffzy=ffy.*filter';
%ffzy=ffy.*(ZBfilter(fy.f,fcenter1,fwidth1)+ZBfilter(fy.f,fcenter2,fwidth2)+ZBfilter(fy.f,fcenter3,fwidth3))'/3;
ffzy=ffzy/max(ffzy);

ch2=sum((ffzy-spec).^2);

plot(fy.f,ffzy,fy.f,spec);
drawnow;
