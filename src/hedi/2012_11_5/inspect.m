function inspect(sex,bird,call,t,a,L)
fname=sprintf('../output/res.pAv4.%cbird%dcall%d.mat',sex,bird,call);
load(fname);
slope=-0.6025;
intercept=-0.5979;

x=res.data(t,:);
alpha0=x(2);
beta0=slope*alpha0 + intercept;
x1=x(3)
y1=x(4)
x2=x(5)
y2=x(6)
x3=x(7)
y3=x(8)



tmax=0.5;
gamma=24000;
ttime=linspace(0,tmax,10);
alpha=alpha0+0*ttime;
beta=beta0+0*ttime;           
y=smBGAs(tmax,gamma,alpha',ttime',beta',ttime');
y=y(end/2:end);

fy=timefreq(y,44100,'stft');
ffy=mean(fy.spec,2);

c=340;

fil=abs(interp1([-10 x1 x2 x3 10000],[0 y1 y2 y3 0],fy.f,'cubic'));
fw=c/(4*L)*1000
tau=2*L/c
1/tau

ffz=ffy.*fil';
ffz=ffz/max(ffz);

filL=abs(1./(sqrt(2+a*2*cos(2*pi*tau*fy.f/1000.0)))).^2;


[af,fm]=max(fil);
w=1500.0;

%fil2=exp(-(fy.f-fy.f(fm)).^2/(2*w^2));
fil2=exp(-(fy.f-fw).^2/(2*a^2));

ffzL=ffy.*filL';
ffzL=ffzL/max(ffzL);

ffz2=ffy.*fil2';
ffz2=ffz2/max(ffz2);


subplot(2,2,1)
plot(fy.f,ffz,fy.f,ffzL);
subplot(2,2,2)
plot(fy.f,fil/max(fil),fy.f,filL/max(filL));
subplot(2,2,3)
plot(fy.f,ffz,fy.f,ffz2);
subplot(2,2,4)
plot(fy.f,fil/max(fil),fy.f,filL/max(filL));
drawnow;