function mz=findAlphaBeta(f,spec,tab,divisor);
% give best alpha and beta using brute force methods



[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax)
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 

nd=floor(fmax/fmin);
if nd<1; 
    nd=1;
end;

ftarget=fmax/divisor;
fm=ftarget-fband;
fM=ftarget+fband;
cho=find(tab(:,3)>fm & tab(:,3)<fM);
atab=tab(cho,1);
btab=tab(cho,2);
n=length(cho)
mz=[];

for i=1:n;
    alpha0=atab(i);
    beta0=btab(i);
    fcenter=3500;
    fwidth=1000;
    x0=[alpha0,beta0,fcenter,fwidth];
    g=@(x)CHI2spec(x,spec);
    options = optimset('MaxIter',100);
    xf=fminsearch(g,x0,options);
    mz=[mz;xf g(xf)];
    xf(1)
    xf(2)
    xf(3)
    xf(4)
    g(xf)
end;



% rz=[];
% ndz=size(dz,1);
% for k=1:ndz;
%     alpha0=dz(k,2);
%     beta0=dz(k,3);
%     gamma=24000;
%     ttime=linspace(0,0.15,10);
%     alpha=alpha0+0*ttime;
%     beta=beta0+0*ttime;
%     y=smBGAs(0.15,24000,alpha',ttime',beta',ttime');
%     fy=timefreq(y,44100,'stft');
%     ffy=mean(fy.spec,2);
%     fcenter=3000;
%     fwidth=1000;
%     ffzy=ffy.*ZBfilter(f,fcenter,fwidth)';
%     rz=[rz,ffzy/max(ffzy)];
% end;
% plot(f,spec,f,rz);