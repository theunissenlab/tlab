function sxpa=GetFitCallL(sex,bird,call,x0);
% try to fit a whole call spec with time vector t and frequency f

fname=sprintf('../data/%cbird%dcall%d.wav',sex,bird,call);
x=wavread(fname);
samplerate=44100.0; % always assumed 
z=timefreq(x,samplerate,'stft');
ev=sum(z.spec);
z.t=z.t(ev>max(ev)/10);
z.spec=z.spec(:,ev>max(ev)/10);

spec=z.spec;
fname=sprintf('../output/res.pAv4.%cbird%dcall%d.mat',sex,bird,call);
load(fname);

slope=-0.6025;
intercept=-0.5979;
tmax=z.t(end);


n=length(z.t);

xna=[];
for tt=fix(n/2):1:n;
    x=res.data(tt,:);
    alpha0=x(2);
    sp=spec(:,tt);
    xpa=fitProcL(sp,x0,xna,z.f,tt,slope,intercept,alpha0);
    xna=[xna;tt xpa];
    
end;
for tt=fix(n/2)-1:-1:1;
    sp=spec(:,tt);
    x=res.data(tt,:);
    alpha0=x(2);
    xpa=fitProcL(sp,x0,xna,f,tt,slope,intercept,alpha0);
    xna=[xna;tt xpa];
end;

sxpa=sortrows(xna,1);

save(hashname,'sxpa');
