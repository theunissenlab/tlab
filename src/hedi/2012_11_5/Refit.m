function Refit(sex,bird,call,a,L)
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


% recomute the wav
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

fw=c/(4*L)*1000;

filter=exp(-(fy.f-fw).^2/(2*a^2));

ffzL=ffy.*filL';
ffzL=ffzL/max(ffzL);

ffz2=ffy.*fil2';
ffz2=ffz2/max(ffz2);

