function [a,b,mm]=diplayFreqs(f,spec,tab);
% compute the correlation of a given spectrum with the one found with the
% model


[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax);
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 

nd=floor(fmax/fmin);
if nd<1; 
    nd=1;
end;
dz=[];

for divisor=2:nd;
    divisor
    ftarget=fmax/divisor;
    fm=ftarget-fband;
    fM=ftarget+fband;
    
    cho=find(tab(:,3)>fm & tab(:,3)<fM);
    atab=tab(cho,1);
    btab=tab(cho,2);
    n=length(cho);
    mz=[];
    MZ=1e6;
    for i=1:n;
        alpha0=atab(i);
        beta0=btab(i);
        gamma=24000;
        ttime=linspace(0,0.15,10);
        alpha=alpha0+0*ttime;
        beta=beta0+0*ttime;           
        y=smBGAs(0.15,gamma,alpha',ttime',beta',ttime');
        fy=timefreq(y,44100,'stft');
        ffy=mean(fy.spec,2);
        fcenter=3500;
        fwidth=1000;
        ffzy=ffy.*ZBfilter(f,fcenter,fwidth)';
        ffzy=ffzy/max(ffzy);
        %z=((ffzy-mean(ffzy))/std(ffzy)).*((spec-mean(spec))/std(spec));
        z=(ffzy-spec).^2;
        mz=[mz;sum(z)];
        if sum(z)<MZ;
            MZ=sum(z);
            plot(fy.f,ffzy,fy.f,spec);
            drawnow
        end;
    end;    
    [cmax,cimax]=min(mz);
    if n==0;
        dz=[dz; divisor 0 0 1e5];
    else
        dz=[dz; divisor atab(cimax) btab(cimax) cmax];
    end
end;
[mm,cm]=min(dz(:,4));
a=dz(cm,2);
b=dz(cm,3);


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