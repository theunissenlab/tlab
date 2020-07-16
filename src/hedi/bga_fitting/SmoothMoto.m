function [new_moto,new_nw,new_sp]=SmoothMoto(moto,spec);


time=moto(:,1);
n=length(time)
nfix=20;
freq=200:10000;
nbad=10;
slope=-0.6025;
intercept=-0.5979;
min_L2=1e9;
xchose=1:n;
for nparam=5:min(35,n-nfix-nbad);
    for ncase=1:30;
        % choose param
        [v,ns]=sort(moto(:,end));
        nchose=ns(1:nfix);
        nperm=randperm(n-nfix-nbad)+nfix;       
        nchose=[nchose;ns(nperm(1:nparam))];
        all_nd=GetSampleCall(nchose,moto,spec);
        if all_nd<min_L2;
            min_L2=all_nd;
            xchose=nchose;
        end;
    end;
end;
%xchose
[all_nd,new_nw,new_sp,new_moto]=GetSampleCall(xchose,moto,spec);

% all_nd
subplot(2,2,1)
plot(moto(:,1),new_moto(:,2),moto(:,1),moto(:,2),'--')
drawnow
subplot(2,2,2)
plot(moto(:,1),new_moto(:,3),'k',moto(:,1),new_moto(:,3)+new_moto(:,4),'--r',new_moto(:,1),new_moto(:,3)-new_moto(:,5),'--r');
drawnow
subplot(2,2,3)
imagesc(spec.t,spec.f,spec.spec)
axis xy
drawnow
subplot(2,2,4)
imagesc(new_sp.t,new_sp.f,new_sp.spec)
axis xy
drawnow
    