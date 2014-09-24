function xhi2=Chi2(spec,alpha,mu,sigma1,sigma2,pwr,slope,intercept);

addpath('../2012_11_5');

beta=slope*alpha+intercept;
tmax=0.5;
gamma=24000;

ttime=linspace(0,tmax,10);

alpha=alpha+0*ttime;
beta=beta+0*ttime;           


y=smBGAs(tmax,gamma,alpha',ttime',beta',ttime');
y=y(end/2:end);


fy=htimefreq(y,44100,'stft');
mui=max(mu-fy.f,0).^2/(2*sigma2^2)+max(fy.f-mu,0).^2/(2*sigma1^2);
filter=abs(pwr)*exp(-mui);
fys=mean(fy.spec,2);
%sfys=fys;
fys=fys.*filter'; 
fys=fys/max(fys)*max(spec);

w=1;%(fy.f-fy.f(1)).*(fy.f(end)-fy.f);w=w/max(w);w=w.^0.5;
%plot(fy.f,spec/max(spec),fy.f,fys/max(fys),fy.f,filter/pwr);
%plot(fy.f,filter/pwr);
drawnow
alpha_min=-0.4115;
lambda=0;
if alpha>alpha_min;
    lambda=10;
end;
xhi2=sum((fys-spec).^2.*w')+lambda;